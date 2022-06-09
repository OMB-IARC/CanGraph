#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: MIT

# This is just a collection of functions used by the "main" script

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from urllib.request import urlopen   # Extensible library for opening URLs
from zipfile import ZipFile          # Work with ZIP files
import os                            # Integration with the system
import xml.etree.ElementTree as ET   # To parse and split XML files
import re                            # To split XML files with a regex pattern

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

def repeat_transaction(tx, num_retries, session, bar, number=None):
    """
    A function that repeats transactions whenever an error is found.
    This may make an incorrect script unnecessarily repear; however, since the error is printed,
    one can discriminate those out, and the function remains helpful to prevent WikiData Read Time-Outs.
    """
    for attempt in range(num_retries):
        try:
            session.write_transaction(tx, number)
            bar()
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

# ********* Edit the Neo4J Database ********* #

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
    NOTE: To run from a driver, you must force the autocommit using: ```session.run( misc.clean_database() )```
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

# ********* Work with Files ********* #

def export_graphml(tx, exportname):
    """
    Exports a Neo4J graph to GML format. The graph will be exported to Neo4JImportPath
    NOTE: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.export.graphml.all("{exportname}", {{useTypes:true, storeNodeIds:false}})
        """)


def import_graphml(tx, exportname):
    """
    Imports a GML file into a Neo4J graph. The file has to be located in Neo4JImportPath
    NOTE: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.import.graphml("{exportname}", {{useTypes:true, storeNodeIds:false, readLabels:True}})
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
    Splits a given .xml file in n smaller XML files, one for each first-level
    record on the original file. This allows for slower processing.
    """
    splitted = re.split(f"<{filetype}>", open(filename).read())
    matches = re.findall(f"<{filetype}>", open(filename).read())
    for index, subunit in enumerate(splitted):
        newfile = filename.split(".")[0] + "_" + str(index) + ".xml"
        with open(newfile, "w+") as f:
            if index > 0:
                f.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
                f.write(f'<{bigtag}>\n')
                f.write(f'{matches[index-1]}')
                f.write(splitted[index])
                if index < len(splitted) -1: f.write(f'</{bigtag}>')
        f.close()
    return len(splitted) - 1
