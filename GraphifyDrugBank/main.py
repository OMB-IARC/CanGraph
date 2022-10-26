#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that leverages the functions present in the :obj:`~CanGraph.GraphifyDrugBank.build\_database`
module to recreate `the DrugBank database <https://go.drugbank.com/>`_ using a graph format
and Neo4J, and then provides an GraphML export file.

Please note that, to work, the functions here pre-suppose you have a DrugBank account with access to the whole Database.
You have to prevously apply for it at: https://go.drugbank.com/public_users/sign_up, and place the ```full_database.xml```
file at ```./xmlfolder```.

For more details on how to run this script, please consult the package's README
"""

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem
from time import sleep               # A hack to avoid starving the system resources

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

    instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
    driver = GraphDatabase.driver(instance, auth=(user, passwd))

    Neo4JImportPath = misc.get_import_path(driver)

    print("Connected to Neo4J")

    with driver.session() as session:
        session.run( misc.clean_database() )

    print("Cleaned DataBase")

    filename = sys.argv[4]
    filepath = os.path.abspath(f"./xmlfolder/{filename}")
    shutil.copyfile(filepath, f"{Neo4JImportPath}/{filename}")

    print("File downloaded, now splitting its contents...")
    total_subfiles = misc.split_xml(f"{Neo4JImportPath}/{filename}", "drug type=.*", "drugbank")

    print("File splitted, now on to import its contents to Neo4J...")
    sleep(5)

    with alive_bar((total_subfiles-1)*13 + 3) as bar:
        for i in range(1, total_subfiles):
            newfile = f"{filename.split('.')[0]}_{i}.xml"
            # For each file, add its contents to the DB
            build_database.build_from_file(newfile, driver)
            # And remove it
            os.remove(f"{Neo4JImportPath}/{filename.split('.')[0]}_{i}.xml")
            sleep(1) # Cool-off time

        # At the end, purge the database
        misc.purge_database(driver)

        # And export it:
        with driver.session() as session:
            session.execute_write(misc.export_graphml, "graph.graphml")
            bar()

    print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
    shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
    print(f"A copy of the file has been saved in this project's work directory")

if __name__ == '__main__':

    main()
