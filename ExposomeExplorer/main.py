#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that leverages the functions present in the :obj:`~CanGraph.ExposomeExplorer.build\_database`
module to recreate `the exposome-explorer database <http://exposome-explorer.iarc.fr>`_ using a graph format
and Neo4J, and then provides an GraphML export file.

Please note that, to work, the functions here pre-suppose you have access to Exposome-Explorer internal CSVs,
and that you have placed them under a folder provided as ```sys.argv[4]```. These CSVs are confidential,
and can only be accessed under request to the *International Agency for Research on Cancer*.

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

def main(args):
    """
    The function that executes the code
    """
    with alive_bar(25) as bar:

        instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
        driver = GraphDatabase.driver(instance, auth=(user, passwd))

        Neo4JImportPath = misc.get_import_path(driver)
        print("Connected to Neo4J")
        bar()

        with driver.session() as session:
            session.run( misc.clean_database() )
        print("Cleaned DataBase")
        bar()

        shutil.copyfile(f"{os.path.abspath(sys.argv[4])}/components.csv", f"{Neo4JImportPath}/components.csv")

        with driver.session() as session:
            misc.manage_transaction(ExposomeExplorerDatabase.add_components, "components.csv")
        os.remove(f"{Neo4JImportPath}/components.csv")

        with driver.session() as session:
            build_database.build_from_file(sys.argv[4], Neo4JImportPath, driver, bar,
                                           do_all = True, keep_counts_and_displayeds = False)

        # At the end, purge the database
        misc.purge_database(driver)

        # And export it:
        with driver.session() as session:
            misc.manage_transaction(misc.export_graphml, "graph.graphml")

    print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
    shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
    print(f"A copy of the file has been saved in this project's work directory")

if __name__ == '__main__':

    main()
