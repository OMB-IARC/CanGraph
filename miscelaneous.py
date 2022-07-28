#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

# This is just a collection of functions used by the "main" script

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from urllib.request import urlopen   # Extensible library for opening URLs
from zipfile import ZipFile          # Work with ZIP files
import os                            # Integration with the system
import xml.etree.ElementTree as ET   # To parse and split XML files
import re                            # To split XML files with a regex pattern
from time import sleep               # Cute go slowly

from MeSHandMetaNetX import build_database as MeSHandMetaNetXDataBases

# ********* Manage the Neo4J Database Connection and Transactions ********* #

def get_neo4j_path(tx):
    """ Runs an autocommit transaction to find Neo4J's import path """
    result = tx.run("Call dbms.listConfig() YIELD name, value WHERE name='dbms.directories.import' RETURN value")
    return [record["value"] for record in result]

def get_import_path(current_driver):
    """ Runs the autocommit transaction and returns the Import Path """
    with current_driver.session() as session:
        Neo4JImportPath = session.read_transaction(get_neo4j_path)[0]
    return Neo4JImportPath

def repeat_transaction(tx, num_retries, session, bar, number=None, query = "WikiData"):
    """
    A function that repeats transactions whenever an error is found.
    This may make an incorrect script unnecessarily repear; however, since the error is printed,
    one can discriminate those out, and the function remains helpful to prevent WikiData Read Time-Outs.
    """
    for attempt in range(num_retries):
        try:
            session.write_transaction(tx, number, query = "WikiData")
            if attempt > 0: print(f"Error solved on attempt #{attempt}")
            break
        except Exception as error:
            if attempt < (num_retries - 1):
                print(f"The following error was found while processing Nodes ending in #{number}")
                print(error)
                print(f"Retrying... ({attempt + 1}/{num_retries})")
            else:
                print(f"{num_retries} consecutive attempts were made at processing Nodes ending in #{number}. Aborting...")
                raise error

# ********* Interact with the Neo4J Database ********* #

def call_db_schema_visualization(tx):
    """
    Shows the DB Schema. TODO: Make it download the image
    """
    return tx.run("""
        CALL db.schema.visualization()
        """)

def clean_database():
    """
    Gets all the nodes in a Neo4J database and removes them
    NOTE: If running directly from neo4j browser, you will need to add ```:auto ```
    at the beginning of the query in order to force it to autocommit
    NOTE: To run from a driver, you must force the autocommit using: ```session.run( clean_database() )```
    """
    return f"""
        MATCH (n)
        CALL {{ WITH n
        DETACH DELETE n
        }} IN TRANSACTIONS OF 1000 ROWS;
        """

def create_n10s_graphconfig(tx):
    """
    Creates a neosemantic constraint to hold all the RDF we will import.
    Taken from: https://neo4j.com/labs/neosemantics/how-to-guide/
    """
    return tx.run("""
        CALL n10s.graphconfig.init({
            handleVocabUris: 'MAP',
            handleMultival: 'ARRAY',
            keepLangTag: true,
            keepCustomDataTypes: true,
            applyNeo4jNaming: true
        })
        """)

def remove_n10s_graphconfig(tx):
    """
    Removes the "_GraphConfig" node, which is necessary for querying SPARQL endpoints
    but not at all useful in our final export
    """
    return tx.run("""
        MATCH (n:`_GraphConfig`) DETACH DELETE n
        """)

def remove_duplicate_relationships(tx):
    """
    Removes duplicated relationships between ANY existing pair of nodes.
    NOTE: Only deletes DIRECTED relationships between THE SAME nodes, combining their properties
    Taken from: https://stackoverflow.com/questions/18724939/neo4j-cypher-merge-duplicate-relationships
    """
    return tx.run("""
            MATCH (s)-[r]->(e)
            WITH s, e, type(r) as typ, collect(r) as rels
            CALL apoc.refactor.mergeRelationships(rels, {properties:"combine"})
            YIELD rel
            RETURN rel
        """)

def remove_duplicate_nodes(tx, node_type, condition, optional_where=""):
    """
    Removes any two nodes of any given ```node_type``` with the same ```condition```.
    This condition can be expressed as: n.{node_property} as a, n.{other_property} as b
    NOTE: Optionally, one can specify mandatory conditions with a ```optional_where``` clause
    """
    if len(str(node_type))>0:
        return tx.run(f"""
                MATCH (n:{node_type})
                {optional_where}
                WITH {condition}, COLLECT(n) AS ns
                WHERE size(ns) > 1
                        CALL apoc.refactor.mergeNodes(ns, {{properties:"combine"}}) YIELD node
                RETURN node;
            """)
    else:
        return tx.run(f"""
                MATCH (n)
                {optional_where}
                WITH {condition}, COLLECT(n) AS ns
                WHERE size(ns) > 1
                        CALL apoc.refactor.mergeNodes(ns, {{properties:"combine"}}) YIELD node
                RETURN node;
            """)

def purge_database(driver):
    """
    A series of commands that purge a database, removing unnecessary, duplicated or empty nodes and merging those without necessary properties
    This has been converted into a common function to standarize the ways the nodes are merged.
    # WARNING: When modifying, take good care on how the keys names are written: if a key is not present, all nodes will be merged!
    """
    with driver.session() as session:
        # Fist, we purge Publications by PubMed_ID, using the abstract to merge those that have no PubMed_ID
        session.write_transaction(remove_duplicate_nodes, "Publication", "n.Pubmed_ID as id", "WHERE n.Pubmed_ID IS NOT null")
        session.write_transaction(remove_duplicate_nodes, "Publication", "n.Abstract as abs", "WHERE n.Abstract IS NOT null AND n.Pubmed_ID IS null")

        # Now, we work on Proteins/Metabolites/Drugs:
        # We merge those that have the same InChI or InChIKey:
        session.write_transaction(remove_duplicate_nodes, "", "n.InChI as inchi", "WHERE n:Protein OR n:Metabolite OR n:Drug AND n.InChI IS NOT null")
        session.write_transaction(remove_duplicate_nodes, "", "n.InChIKey as key", "WHERE n:Protein OR n:Metabolite OR n:Drug AND n.InChIKey IS NOT null")
        # We merge Proteins by UniProt_ID, and, when there is none, by Name:
        session.write_transaction(remove_duplicate_nodes, "Protein", "n.UniProt_ID as id", "WHERE n.UniProt_ID IS NOT null")
        session.write_transaction(remove_duplicate_nodes, "Protein", "n.Name as id", "WHERE n.UniProt_ID IS null AND n.Name IS NOT null")
        # We merge Metabolites by HMDB_ID (normally, it should be unique):
        session.write_transaction(remove_duplicate_nodes, "", "n.HMDB_ID as hmdb_id", "WHERE n:Protein OR n:Metabolite AND n.HMDB_ID IS NOT null")
        # We can also remove all OriginalMetabolites (or other) that are non-unique
        session.write_transaction(remove_duplicate_nodes, "", """n.InChI as inchi, n.InChIKey as inchikey, n.Name as name, n.SMILES as smiles,
                                                                      n.Identifier as hmdb_id, n.ChEBI as chebi, n.Monisotopic_Molecular_Weight as mass""",
                                                                   "WHERE n:Metabolite OR n:Protein OR n:OriginalMetabolite OR n:Drug")
        # And all Metabolites that do not match our Schema
        session.run("MATCH (m:Metabolite) WHERE n.ChEBI_ID IS NULL OR n.ChEBI_ID = "" OR n.CAS_Number IS NULL OR n.CAS_Number = "" DETACH DELETE m")

        # We also remove all non-unique Subjects. We do this by passing on all three parameters this nodes may have to apoc.mergeNodes
        # NOTE: This concerns only those nodes that DO NOT COME from Exposome_Explorer
        session.write_transaction(remove_duplicate_nodes, "Subject", "n.Age_Mean as age, n.Gender as gender, n.Information as inf", "WHERE n.Exposome_Explorer_ID IS null")

        # We can do the same for the different Dosages:
        session.write_transaction(remove_duplicate_nodes, "Dosage", "n.Form as frm, n.Stength as str, n.Route as rt")

        # For products, we merge all those with the same EMA_MA_Number or FDA_Application_Number to try to minimize duplicates, although this is not the best approach
        session.write_transaction(remove_duplicate_nodes, "Product", "n.EMA_MA_Number as ema_nb", "WHERE n.EMA_MA_Number IS NOT null")
        session.write_transaction(remove_duplicate_nodes, "Product", "n.FDA_Application_Number as fda_nb", "WHERE n.FDA_Application_Number IS NOT null")

        # For CelularLocations and BioSpecimens, we merge those with the same Name:
        session.write_transaction(remove_duplicate_nodes, "CelularLocation", "n.Name as name")
        session.write_transaction(remove_duplicate_nodes, "BioSpecimen", "n.Name as name")

        # Finally, we delete all empty nodes. This shouldn't be created on the first place, but, in case anyone escapes, this makes the DB cleaner.
        # NOTE: In the case of Taxonomys, these "empty nodes" are actually created on purpose. This, they are here removed.
        session.run("MATCH (n) WHERE size(keys(properties(n))) < 1 CALL { WITH n DETACH DELETE n } IN TRANSACTIONS OF 1000 ROWS")
        # For Measurements and Sequences, 2 properties are the minimum, since they always have some boolean values
        session.run("MATCH (m:Measurement) WHERE size(keys(properties(m))) < 2 DETACH DELETE m")
        session.run("MATCH (s:Sequence) WHERE size(keys(properties(s))) < 2 DETACH DELETE s")

        #At last, we may remove any duplicate relationships, which, since we have merged nodes, will surely be there:
        session.write_transaction(remove_duplicate_relationships)

def find_synonyms(driver, hmdb_ids, chebi_ids, inchis, names, query_type, query):
    with driver.session() as session:
        graph_response = session.read_transaction(MeSHandMetaNetXDataBases.read_synonyms_in_metanetx, query_type, query)
        for element in graph_response:
            if element["databasename"].lower() == "hmdb":
                if element["databaseid"] not in hmdb_ids: hmdb_ids.append(element["databaseid"])
            if element["databasename"].lower() == "chebi":
                if element["databaseid"] not in chebi_ids: chebi_ids.append(element["databaseid"])
            if element["InChI"] not in inchis: inchis.append(graph_response[0]["InChI"])
            if element["Name"] not in names: names.append(graph_response[0]["Name"])

# ********* Work with Files ********* #

def export_graphml(tx, exportname):
    """
    Exports a Neo4J graph to GML format. The graph will be exported to Neo4JImportPath
    NOTE: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.export.graphml.all("{exportname}", {{useTypes:true, storeNodeIds:false}})
        """)


def import_graphml(tx, importname):
    """
    Imports a GML file into a Neo4J graph. The file has to be located in Neo4JImportPath
    NOTE: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.import.graphml("{importname}", {{useTypes:true, storeNodeIds:false, readLabels:True}})
        """)

def download_and_unzip(url, folder):
    """
    https://svaderia.github.io/articles/downloading-and-unzipping-a-zipfile/
    https://stackoverflow.com/questions/32123394/workflow-to-create-a-folder-if-it-doesnt-exist-already
    """
    zipresp = urlopen(url)
    tempzip = open("/tmp/tempfile.zip", "wb")
    tempzip.write(zipresp.read())
    tempzip.close()
    zf = ZipFile("/tmp/tempfile.zip")
    if not os.path.exists(f"./{folder}"):
        os.makedirs(f"./{folder}")
    zf.extractall(path = f"./{folder}/")
    zf.close()

def split_xml(filename, filetype, bigtag):
    """
    Splits a given .xml file in n smaller XML files, one for each first-level record (the name of which has to be specified)
    on the original file, so that it does not crash when processing it.
    """
    # Set counters and text variables to null
    current_text = ""; current_line = 0; num_files = 0
    with open(f'{filename}', "r") as f:
        for line in f:
            current_text += line; current_line += 1
            # Whenever we find a closing tag, we know the info is over
            if line.startswith(f"</{filetype}>"):
                newfile = filename.split(".")[0] + "_" + str(num_files) + ".xml"
                with open(newfile, "w+") as f:
                    # Skip already-present headers for the first file
                    if num_files > 1:
                        f.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
                        f.write(f'<{bigtag}>\n')
                    f.write(current_text)
                    # This is the same since it will just stop processing at the last tag
                    f.write(f'</{bigtag}>\n')
                f.close()
                num_files += 1
                current_text = ""
    return num_files - 1
