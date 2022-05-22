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
    result = tx.run("Call dbms.listConfig() YIELD name, value WHERE name='dbms.directories.import' RETURN value")
    return [record["value"] for record in result]

def get_import_path(current_driver):
    """ Runs the autocommit transaction and returns the Import Path """
    with current_driver.session() as session:
        Neo4JImportPath = session.read_transaction(get_neo4j_path)[0]
    return Neo4JImportPath

def export_function(tx, exportname):
    """
    Exports a Neo4J graph to GML format. The graph will be exported to Neo4JImportPath
    Note: for this to work, you HAVE TO have APOC availaible on your Neo4J installation
    """
    return tx.run(f"""
        CALL apoc.export.graphml.all("{exportname}", {{useTypes:true, storeNodeIds:false}})
        """)

