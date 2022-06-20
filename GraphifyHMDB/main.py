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
from time import sleep               # A hack to avoid starving the system resources

# Import subscripts for the program
import build_database
# A hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# NOTE: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

hmdb_urls = ["https://hmdb.ca/system/downloads/current/hmdb_proteins.zip",
             "https://hmdb.ca/system/downloads/current/urine_metabolites.zip",
             "https://hmdb.ca/system/downloads/current/serum_metabolites.zip",
             "https://hmdb.ca/system/downloads/current/csf_metabolites.zip",
             "https://hmdb.ca/system/downloads/current/saliva_metabolites.zip",
             "https://hmdb.ca/system/downloads/current/feces_metabolites.zip",
             "https://hmdb.ca/system/downloads/current/sweat_metabolites.zip"
            ]

instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
driver = GraphDatabase.driver(instance, auth=(user, passwd))

Neo4JImportPath = misc.get_import_path(driver)

print("Connected to Neo4J")

with driver.session() as session:
    session.run( misc.clean_database() )

print("Cleaned DataBase")

for url in hmdb_urls:

    filename = f"{url.split('/')[-1].split('.')[0]}.xml"
    print(f"Downloading and Unzipping: {filename}...")
    misc.download_and_unzip(url, "xmlfolder")

    filepath = os.path.abspath(f"./xmlfolder/{filename}")
    shutil.copyfile(filepath, f"{Neo4JImportPath}/{filename}")

    if filename == "hmdb_proteins.xml":

        print("File downloaded, now splitting its contents...")
        total_subfiles = misc.split_xml(f"{Neo4JImportPath}/{filename}", "protein", "hmdb")

        print("File splitted, now on to import its contents to Neo4J...")
        sleep(5)

        with alive_bar(total_subfiles*7) as bar:
            for i in range(1, total_subfiles):
                newfile = f"{filename.split('.')[0]}_{i}.xml"
                # For each file, add its contents to the DB
                build_database.build_from_protein_file(newfile, driver)
                # And remove it
                os.remove(f"{Neo4JImportPath}/{filename.split('.')[0]}_{i}.xml")

    else:

        print("File downloaded, now splitting its contents...")
        total_subfiles = misc.split_xml(f"{Neo4JImportPath}/{filename}", "metabolite", "hmdb")

        print("File splitted, now on to import its contents to Neo4J...")
        sleep(5)

        with alive_bar(total_subfiles*9) as bar:
            for i in range(1, total_subfiles):
                newfile = f"{filename.split('.')[0]}_{i}.xml"
                # For each file, add its contents to the DB
                build_database.build_from_metabolite_file(newfile, driver)
                # And remove it
                os.remove(f"{Neo4JImportPath}/{filename.split('.')[0]}_{i}.xml")

# At the end, purge the database
build_database.purge_database(driver)

with driver.session() as session:
    session.write_transaction(misc.export_graphml, "graph.graphml")
    bar()

print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
print(f"A copy of the file has been saved in this project's directory")
