#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem
from time import sleep               # A hack to avoid starving the system resources

# Import subscripts for the program
import create_nodes
import create_relations
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
    session.write_transaction(create_nodes.clean_database)

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
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_proteins, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_go_classifications, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_gene_properties, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_protein_properties, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_relations.add_metabolite_associations, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_relations.add_metabolite_references, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_general_references, newfile, "protein")
                    bar()
                sleep(1)
                os.remove(f"{Neo4JImportPath}/{filename.split('.')[0]}_{i}.xml")

    else:

        print("File downloaded, now splitting its contents...")
        total_subfiles = misc.split_xml(f"{Neo4JImportPath}/{filename}", "metabolite")

        print("File splitted, now on to import its contents to Neo4J...")
        sleep(5)

        with alive_bar(total_subfiles*9) as bar:
            for i in range(1, total_subfiles):
                newfile = f"{filename.split('.')[0]}_{i}.xml"
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_metabolites, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_relations.add_protein_associations, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_diseases, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_concentrations_normal, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_concentrations_abnormal, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_taxonomy, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_biological_properties, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_experimental_properties, newfile)
                    session.write_transaction(create_nodes.add_predicted_properties, newfile)
                    bar()
                with driver.session() as session:
                    session.write_transaction(create_nodes.add_general_references, newfile, "metabolite")
                    bar()
                sleep(1)
                os.remove(f"{Neo4JImportPath}/{filename.split('.')[0]}_{i}.xml")

with driver.session() as session:
    session.write_transaction(create_nodes.remove_duplicate_accessions)
    session.write_transaction(create_nodes.remove_duplicate_pubmed_ids)
    session.write_transaction(create_nodes.remove_duplicate_abstracts)
    session.write_transaction(create_nodes.remove_duplicate_relationships)
    session.write_transaction(create_nodes.purge_database)
    bar()
with driver.session() as session:
    session.write_transaction(misc.export_graphml, "graph.graphml")
    bar()

print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
print(f"A copy of the file has been saved in this project's directory")
