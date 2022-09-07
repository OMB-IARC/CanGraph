#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that provides a collection of functions to be used across the different
scripts present in the CanGraph package, with various, useful functionalities
"""

# Import external modules necessary for the script
from neo4j import GraphDatabase     # The Neo4J python driver
import urllib.request as request    # Extensible library for opening URLs
from zipfile import ZipFile         # Work with ZIP files
import tarfile                      # Work with tar.gz files
import os                           # Integration with the system
import xml.etree.ElementTree as ET  # To parse and split XML files
import re                           # To split XML files with a regex pattern
from time import sleep              # Cute go slowly
import pandas as pd                 # Analysis of tabular data

from MeSHandMetaNetX import build_database as MeSHandMetaNetXDataBases

# ********* Manage the Neo4J Database Connection and Transactions ********* #

def neo4j_import_path_query(tx):
    """
    A CYPHER query able to find Neo4J's Import Path. When used in :obj:`~CanGraph.miscelaneous.get_import_path`,
    it returns a one-item list with the import path itself, but it **needs** to be executed in
    :obj:`~CanGraph.miscelaneous.get_import_path` for it to do so.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running

    Returns:
        str: A one-item list containing the Neo4J Import Path
    """
    result = tx.run("Call dbms.listConfig() YIELD name, value WHERE name='dbms.directories.import' RETURN value")
    return [record["value"] for record in result]

def get_import_path(current_driver):
    """
    A function that runs :obj:`~CanGraph.miscelaneous.neo4j_import_path_query` to get Neo4J's Import Path

    .. NOTE:: By doing the Neo4JImportPath search this way (in two functions), we are able to run the query as
        a :obj: read_transaction, which, unlike autocommit transactions, allows the query to be better controlled,
        and repeated in case it fails.

    Args:
        current_driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        str: Neo4J's Import Path, i.e., where Neo4J will pick up files to be imported using the ```file:///``` schema
    """
    with current_driver.session() as session:
        Neo4JImportPath = session.read_transaction(neo4j_import_path_query)[0]
    return Neo4JImportPath

def repeat_transaction(tx, num_retries, driver, **kwargs):
    """
    A function that repeats transactions whenever an error is found.
    This may make an incorrect script unnecessarily repeat; however, since the error is printed,
    one can discriminate those out, and the function remains helpful to prevent SPARQL Read Time-Outs.

    Args:
        tx (neo4j.work.simple.Session): The session under which the driver is running
        num_retries (int): The number of times that we wish the transaction to be retried
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        **kwargs: Any number of arbitrary keyword arguments

    Raises:
        Exception:
            An exception telling the user that the maximum number of retries
            has been exceded, if such a thing happens
    """
    for attempt in range(num_retries):
        try:
            with driver.session() as session:
                session.write_transaction(tx, **kwargs)
            if attempt > 0: print(f"Error solved on attempt #{attempt}")
            break
        except Exception as error:
            if attempt < (num_retries - 1):
                print(f"The following error was found while processing Nodes ending in #{number}")
                print(error)
                print(f"Retrying... ({attempt + 1}/{num_retries})")
            else:
                raise Exception(f"{num_retries} consecutive attempts were made at processing Nodes ",
                                f"ending in #{number}. Aborting...")

# ********* Interact with the Neo4J Database ********* #

def call_db_schema_visualization(tx):
    """
    Shows the DB Schema. This function is intended to be run only in Neo4J's console, since it produces no output when called from the driver.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running

    .. TODO:: Make it download the image
    """
    return tx.run("""
        CALL db.schema.visualization()
        """)

def clean_database():
    """
    A CYPHER query that gets all the nodes in a Neo4J database and removes them, in transactions of 1000 rows to alleviate memory load

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`

    .. NOTE:: This is an **autocommit transaction**. This means that, in order to not keep data in memory
        (and make running it with a huge amount of data) more efficient, you will need to add ```:auto ```
        when calling it from the Neo4J browser, or call it as using :obj:`neo4j.Session.run` from the driver.
    """
    return """
        MATCH (n)
        CALL { WITH n
        DETACH DELETE n
        } IN TRANSACTIONS OF 1000 ROWS
        """

def create_n10s_graphconfig(tx):
    """
    A CYPHER query that creates a *neosemantics* (n10s) constraint to hold all the RDF we will import.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. seealso:: More information on this approach can be found in `Neosemantics' 101 Guide <https://neo4j.com/labs/neosemantics/>`_
        and in `Neo4J's guide on how to import data from Wikidata <https://neo4j.com/labs/neosemantics/how-to-guide/>`_ ,
        where this approach was taken from

    .. deprecated:: 0.9
        Since we are importing based on apoc.load.jsonParams, this is not needed anymore
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

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. deprecated:: 0.9
        Since we are importing based on apoc.load.jsonParams, this is not needed anymore
    """
    return tx.run("""
        MATCH (n:`_GraphConfig`) DETACH DELETE n
        """)

def remove_ExternalEquivalent(tx):
    """
    Removes all nodes of type: ExternalEquivalent from he DataBase; since this do not add
    new info, one might consider them not useful.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.
    """
    return tx.run("""
        MATCH (e:ExternalEquivalent)
        DETACH DELETE e
        """)

def remove_duplicate_relationships(tx):
    """
    Removes duplicated relationships between ANY existing pair of nodes.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Only deletes DIRECTED relationships between THE SAME nodes, combining their properties

    .. seealso:: This way of working has been taken from
        `StackOverflow #18724939 <https://stackoverflow.com/questions/18724939/neo4j-cypher-merge-duplicate-relationships>`_
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

    Args:
        tx              (neo4j.work.simple.Session): The session under which the driver is running
        node_type       (str): The label of the nodes that will be selected for merging. If more than one, leave null and set
            optional_where = ``WHERE (n:Label1 ÒR n:Label2) AND ...``
        condition       (str): The node properties used for collecting, if not using all properties. This can be expressed as:
            ``n.{property_1} as a, n.{property_n} as b``
        optional_where  (str): An optional ``WHERE`` Neo4J Statement to constrain the nodes that are picked up; default is "" (no statement)

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. WARNING:: When using, take good care on how the keys names are written: sometimes, if a key is not present, all nodes will be merged!
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

     Args:
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.

    .. WARNING:: When modifying, take good care on how the keys names are written: with :obj:`~CanGraph.miscelaneous.remove_duplicate_nodes`,
        sometimes, if a key is not present, all nodes will be merged!
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
        # .. NOTE:: This concerns only those nodes that DO NOT COME from Exposome_Explorer
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
        # .. NOTE:: In the case of Taxonomys, these "empty nodes" are actually created on purpose. This, they are here removed.
        session.run("MATCH (n) WHERE size(keys(properties(n))) < 1 CALL { WITH n DETACH DELETE n } IN TRANSACTIONS OF 1000 ROWS")
        # For Measurements and Sequences, 2 properties are the minimum, since they always have some boolean values
        session.run("MATCH (m:Measurement) WHERE size(keys(properties(m))) < 2 DETACH DELETE m")
        session.run("MATCH (s:Sequence) WHERE size(keys(properties(s))) < 2 DETACH DELETE s")

        #At last, we may remove any duplicate relationships, which, since we have merged nodes, will surely be there:
        session.write_transaction(remove_duplicate_relationships)

def find_synonyms_in_metanetx(driver, query_type, query):
    """
    A function that automates finding synonyms in MetaNetX for a given metabolite,
    using :obj:`CanGraph.MeSHandMetaNetX.build_database.read_synonyms_in_metanetx`.

    Args:
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        query_type (str): one of ["HMDB_ID","ChEBI_ID","InChI","Name"]; the type of query we will be searching synonyms in
        query (str): The value of an identifier based on which we want to find synonyms for a metabolite

    Returns:
        list: a list of the synonyms that MetaNetx can find for ``query``, including the query itself
    """
    # First, we declare a list to store the synonyms and we store the original IDs there
    items_list = []; items_list.append(query)

    # Then, we find the synonyms based on the original query
    with driver.session() as session:
        graph_response = session.read_transaction(MeSHandMetaNetXDataBases.read_synonyms_in_metanetx, query_type, query)

    # And we process it to add it to the appropriate list
    for element in graph_response:

        if element["databasename"].lower() == "hmdb" and query_type == "HMDB_ID":
            if element["databaseid"] not in items_list: items_list.append(element["databaseid"])

        if element["databasename"].lower() == "chebi" and query_type == "ChEBI_ID":
            if element["databaseid"] not in items_list: items_list.append(element["databaseid"])

        if element["InChI"] not in items_list and query_type == "InChI": items_list.append(graph_response[0]["InChI"])

        if element["Name"] not in items_list and query_type == "Name": items_list.append(graph_response[0]["Name"])

    return items_list

# ********* Work with Files ********* #

def export_graphml(tx, exportname):
    """
    Exports a Neo4J graph to GraphML format. The graph will be exported to Neo4JImportPath

    Args:
        tx  (neo4j.work.simple.Session): The session under which the driver is running
        exportname (str): The name for the exported file, which will be saved under ./Neo4JImportPath/

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that exports the file, using batch optimizations and
            smaller batch sizes to try to keep the impact on memory use low

    .. NOTE:: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.export.graphml.all("{exportname}", {{batchSize: 5, useTypes:true, storeNodeIds:false,
                                                       useOptimizations:
                                                           {{type: "UNWIND_BATCH", unwindBatchSize: 5}} }})
        """)


def import_graphml(tx, importname):
    """
    Imports a GraphML file into a Neo4J graph. The file has to be located in Neo4JImportPath

    Args:
        tx  (neo4j.work.simple.Session): The session under which the driver is running
        importname (str): The name for the file to be imported, which must be under ./Neo4JImportPath/

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that imports the file, using batch optimizations and
            smaller batch sizes to try to keep the impact on memory use low

    .. NOTE:: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.import.graphml("{importname}", {{batchSize: 5, useTypes:true, storeNodeIds:false,
                                                   readLabels:True, useOptimizations:
                                                       {{type: "UNWIND_BATCH", unwindBatchSize: 5}} }})
        """)

def download(url, folder):
    """
    Downloads a file from the internet into a given folder

    Args:
        url (str): The Uniform Resource Locator for the Zipfile to be downloaded and unzipped
        folder (str): The folder under which the file will be stored.

    Returns:
        str: The path where the file we just downloaded has been stored
    """

    folder = os.path.abspath(folder) # Set the folder to be an absolute path

    # If the folder does not exist, we create it
    if not os.path.exists(f"{folder}"):
        os.makedirs(f"{folder}")

    # Create some naming variables and request the file from the URL
    filename = url.split('/')[-1].split('.')[0]
    file_ext = url.split('.')[-1]
    file_path = f"{folder}/{filename}.{file_ext}"
    request.urlretrieve(url, file_path)

    return file_path

def unzip(file_path, folder):
    """
    Unizps a file present at a given ``file_path`` into a given ``folder``

    Args:
        url (str): The Uniform Resource Locator for the Zipfile to be unzipped
        folder (str): The folder under which the file will be stored.

    Returns:
        str: The path where the file we just unzipped has been stored
    """

    folder = os.path.abspath(folder) # Set the folder to be an absolute path
    filename = file_path.split('/')[-1].split('.')[0] # Create some naming variables

    # If the folder does not exist, we create it
    if not os.path.exists(f"{folder}"):
        os.makedirs(f"{folder}")

    # And we unzip the file
    zf = ZipFile(file_path)
    zf.extractall(path = f"{folder}/")
    zf.close()

    return filename

def untargz(file_path, folder):
    """
    Untargzs a file present at a given ``file_path`` into a given ``folder``

    Args:
        url (str): The Uniform Resource Locator for the Tarfile to be untargz
        folder (str): The folder under which the file will be stored.

    Returns:
        str: The path where the file we just untargz has been stored
    """

    folder = os.path.abspath(folder) # Set the folder to be an absolute path
    filename = file_path.split('/')[-1].split('.')[0] # Create some naming variables

    # If the folder does not exist, we create it
    if not os.path.exists(f"{folder}"):
        os.makedirs(f"{folder}")

    # And we untargz the file
    tf = tarfile.open(file_path)
    tf.extractall(path = f"{folder}/")
    tf.close()

    return filename

def download_and_unzip(url, folder):
    """
    Downloads and unzips a given Zipfile from the internet; useful for databases which provide zip access.

    Args:
        url (str): The Uniform Resource Locator for the Zipfile to be downloaded and unzipped
        folder (str): The folder under which the file will be stored.

    Returns:
        This function downloads and unzips the file in the desired folder, but does not produce any particular return.

    .. seealso:: Code snippets for this function were taken from `Shyamal Vaderia's Github <https://svaderia.github.io/articles/downloading-and-unzipping-a-zipfile/>`_
        and from `StackOverflow #32123394 <https://stackoverflow.com/questions/32123394/workflow-to-create-a-folder-if-it-doesnt-exist-already>`_
    """
    file_path = download(url, "/tmp")
    unzip(file_path, folder)
    os.remove(file_path)

def download_and_untargz(url, folder):
    """
    Downloads and unzips a given ``tar.gz`` from the internet

    Args:
        url (str): The Uniform Resource Locator for the ``tar.gz`` to be downloaded and unzipped
        folder (str): The folder under which the file will be stored.

    Returns:
        This function downloads and unzips the file in the desired folder, but does not produce any particular return.
    """
    file_path = download(url, "/tmp")
    untargz(file_path, folder)
    os.remove(file_path)

def split_xml(filepath, splittag, bigtag):
    """
    Splits a given .xml file in n smaller XML files, one for each ``splittag`` section that is pressent in the original
    file, which should be of type ``bigtag``. For example, we might have an ``<hmdb>`` file which we want to slit based on
    the ``<metabolite>`` items therein contained. Ths is so that Neo4J does not crash when processing it.

    Args:
        filepath (str): The path to the file that needs to be xplitted
        splittag (str): The tag based on which the file will be split
        bigtag (str): The main tag of the file, which needs to be re-added.

    Returns:
        int: The number of files that have been produced from the original

    .. WARNING:: The original file will be removed
    """

    # Set the filepath to be an absolute path
    filepath = os.path.abspath(filename)

    # Set counters and text variables to null
    current_text = ""; current_line = 0; num_files = 0

    # Open the file that we wish to split (but without reading! neat!)
    with open(f'{filepath}', "r") as f:
        for line in f: # For each line, scan and recount
            current_text += line; current_line += 1

            # Whenever we find a closing tag, we know the info is over
            # and we generate a new file
            if line.startswith(f"</{splittag}>"):
                newfile = filepath.split(".")[0] + "_" + str(num_files) + ".xml"
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

    # Remove the original file
    os.remove(f"{filepath}")

    return num_files - 1

def split_csv(filename, folder, sep=",", sep_out=",", startFrom=0, withStepsOf=1):
    """
    Splits a given .csv/tsv file in n smaller csv files, one for each row on the original file,
    so that it does not crash when processing it. It also allows to start reading from ```startFrom``` lines

    Args:
        filepath (str): The path to the file that needs to be xplitted
        splittag (str): The tag based on which the file will be split
        bigtag (str): The main tag of the file, which needs to be re-added.

    Returns:
        int: The number of files that have been produced from the original

    .. WARNING:: The original file will be removed
    """

    # Set the filepath to be an absolute path
    filepath = os.path.abspath(f"{folder}/{filename}")

    # Read the file using pandas (and hope it does not crash"
    bigfile = pd.read_csv(filepath, sep=sep, skiprows=startFrom)
    filenumber = 0
    for index in range(0, len(bigfile), withStepsOf):
        new_df = bigfile.iloc[index:(index+withStepsOf)]
        new_df.to_csv(f"{folder}/{os.path.splitext(filename)[0]}_{filenumber}.csv", index = False, sep=sep_out)
        filenumber += 1

    os.remove(filepath) # And remove the file after finishing

    return filenumber
