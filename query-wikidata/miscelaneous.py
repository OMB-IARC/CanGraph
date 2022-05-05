#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This is just a collection of functions used by the "main" script

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver

def get_neo4j_path(tx):
    """ Runs an autocommit transaction to find Neo4J's import path """
    result = tx.run("Call dbms.listConfig() YIELD name, value WHERE name='dbms.directories.neo4j_home' RETURN value")
    return [record["value"] for record in result]

def get_import_path(current_driver):
    """ Runs the autocommit transaction and returns the Import Path """
    with current_driver.session() as session:
        Neo4JImportPath = session.read_transaction(get_neo4j_path)[0]+'/import'
    return Neo4JImportPath

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

def remove_n10s_graphconfig(tx):
    """
    Removes the "_GraphConfig" node, which is necessary for querying SPARQL endpoints
    but not at all useful in our final export
    """
    return tx.run("""
        MATCH (n:`_GraphConfig`) DETACH DELETE n
        """)

def call_db_schema_visualization(tx):
    """
    Shows the DB Schema. TODO: Make it download the image
    """
    return tx.run("""
        CALL db.schema.visualization()
        """)

def export_function(tx, exportname):
    """
    Exports a Neo4J graph to GML format. The graph will be exported to Neo4JImportPath
    Note: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.export.graphml.all("{exportname}", {{useTypes:true, storeNodeIds:false}})
        """)

