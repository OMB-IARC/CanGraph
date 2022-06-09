#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: MIT

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem

# Import subscripts for the program
import sparql_queries
# A hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# NOTE: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

with alive_bar(65) as bar:

    instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
    driver = GraphDatabase.driver(instance, auth=(user, passwd))

    Neo4JImportPath = misc.get_import_path(driver)

    print("Connected to Neo4J")

    with driver.session() as session:
        session.run( misc.clean_database() )
        session.write_transaction(misc.create_n10s_graphconfig)
        bar()
    print("Cleaned DataBase")

    with driver.session() as session:
        session.write_transaction(sparql_queries.initial_cancer_discovery)
        session.write_transaction(sparql_queries.remove_duplicate_nodes)
        bar()

    for number in range(3):
        with driver.session() as session:
            misc.repeat_transaction(sparql_queries.find_subclass_of_cancer, 10, session, bar)
            misc.repeat_transaction(sparql_queries.remove_duplicate_nodes, 10, session, bar)

    with driver.session() as session:
        misc.repeat_transaction(sparql_queries.find_instance_of_cancer, 10, session, bar)
        misc.repeat_transaction(sparql_queries.remove_duplicate_nodes, 10, session, bar)

    for number in range(10):
        with driver.session() as session:
            misc.repeat_transaction(sparql_queries.add_cancer_info, 10, session, bar, number)

    for number in range(10):
        with driver.session() as session:
            misc.repeat_transaction(sparql_queries.add_drugs, 10, session, bar, number)
            misc.repeat_transaction(sparql_queries.add_causes, 10, session, bar, number)
            misc.repeat_transaction(sparql_queries.add_genes, 10, session, bar, number)

    with driver.session() as session:
        session.write_transaction(sparql_queries.remove_duplicate_nodes)
        bar()

    with driver.session() as session:
        misc.repeat_transaction(sparql_queries.add_drug_external_ids, 10, session, bar)
        misc.repeat_transaction(sparql_queries.add_more_drug_info, 10, session, bar)

    with driver.session() as session:
        misc.repeat_transaction(sparql_queries.add_gene_info, 10, session, bar)
        bar()

    for number in range(4):
        with driver.session() as session:
            misc.repeat_transaction(sparql_queries.add_metabolite_info, 10, session, bar)
            bar()

    with driver.session() as session:
        #NOTE: We purge by merging all products with the same EMA_MA_Number or FDA_Application_Number
        session.write_transaction(misc.remove_duplicate_nodes, "", "n.WikiData_ID as wdt")
        # And, then, we delete duplicate relationships
        #NOTE: We will just do this once to decrease processing time
        session.write_transaction(misc.remove_duplicate_relationships)
        bar()

    with driver.session() as session:
        # We might want to remove ExternalEquivalent nodes
        #session.write_transaction(sparql_queries.remove_ExternalEquivalent)
        session.write_transaction(misc.remove_n10s_graphconfig)
        session.write_transaction(misc.export_graphml, "graph.graphml")
        bar()

print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
print(f"A copy of the file has been saved in this project's directory")
