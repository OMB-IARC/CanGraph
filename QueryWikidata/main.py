#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that leverages the functions present in the :obj:`~CanGraph.QueryWikidata.build\_database`
module to recreate selected parts of the `the Wikidata database <https://www.wikidata.org/wiki/Wikidata:Main_Page>`_
using a graph format and Neo4J, and then provides an GraphML export file.

Please note that, to work, the functions here pre-suppose you have internet access, which will be used to access
Wikidata's SPAQL endpoint and write info to the Neo4J database

For more details on how to run this script, please consult the package's README
"""

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem

# Import subscripts for the program
import build_database
# A hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# NOTE: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

def main():
    """
    The function that executes the code
    """
    with alive_bar(65) as bar:

        instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
        driver = GraphDatabase.driver(instance, auth=(user, passwd))

        Neo4JImportPath = misc.get_import_path(driver)

        print("Connected to Neo4J")

        with driver.session() as session:
            misc.manage_transaction( misc.clean_database(), driver )
            bar()
        print("Cleaned DataBase")

        with driver.session() as session:
            misc.manage_transaction(build_database.initial_cancer_discovery(), driver)
            bar()

        for number in range(3):
            with driver.session() as session:
                misc.manage_transaction(build_database.find_subclass_of_disease(), driver)
                bar(); bar()

        with driver.session() as session:
            misc.manage_transaction(build_database.find_instance_of_disease(), driver)
            bar(); bar()

        for number in range(10):
            with driver.session() as session:
                misc.manage_transaction(build_database.add_disease_info(), driver, number)
                bar()

        for number in range(10):
            with driver.session() as session:
                misc.manage_transaction(build_database.add_drugs(), driver, number)
                misc.manage_transaction(build_database.add_causes(), driver, number)
                misc.manage_transaction(build_database.add_genes(), driver, number)
                bar(); bar(); bar();

        with driver.session() as session:
            misc.manage_transaction(build_database.add_drug_external_ids(), driver)
            misc.manage_transaction(build_database.add_more_drug_info(), driver)
            bar(); bar()

        with driver.session() as session:
            misc.manage_transaction(build_database.add_gene_info(), driver)
            bar()

        for number in range(4):
            with driver.session() as session:
                misc.manage_transaction(build_database.add_metabolite_info(), driver)
                bar()

        with driver.session() as session:
            misc.manage_transaction(misc.purge_database(), driver)
            misc.manage_transaction(misc.merge_duplicate_nodes("", "n.WikiData_ID as wdt"), driver)
            # And, then, we delete duplicate relationships
            #NOTE: We will just do this once to decrease processing time
            misc.manage_transaction(misc.remove_duplicate_relationships(), driver)
            bar()

        # # At the end, purge the database
        misc.purge_database(driver)

        # And export it:
        with driver.session() as session:
            # We might want to remove ExternalEquivalent nodes
            #misc.manage_transaction(misc.remove_ExternalEquivalent)
            misc.manage_transaction(misc.export_graphml("graph.graphml"), driver)
            bar()

    print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
    shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
    print(f"A copy of the file has been saved in this project's work directory")

if __name__ == '__main__':

    main()
