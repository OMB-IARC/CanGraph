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
import create_nodes
# A hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# NOTE: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

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
        with driver.session() as session:
            session.write_transaction(create_nodes.add_drugs, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_general_references, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_taxonomy, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_products, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_mixtures, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_categories, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_manufacturers, newfile)
            session.write_transaction(create_nodes.add_packagers, newfile)
            session.write_transaction(create_nodes.add_dosages, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_atc_codes, newfile)
            session.write_transaction(create_nodes.add_drug_interactions, newfile)
            session.write_transaction(create_nodes.add_sequences, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_experimental_properties, newfile)
            session.write_transaction(create_nodes.add_external_identifiers, newfile)
            session.write_transaction(create_nodes.add_external_equivalents, newfile)
            bar()
        with driver.session() as session:
            session.write_transaction(create_nodes.add_pathways_and_relations, newfile)
            bar()
        with driver.session() as session:
            for element in ["enzymes", "carriers", "transporters"]:
                session.write_transaction(create_nodes.add_targets_enzymes_carriers_and_transporters, newfile, element)
                bar()
        sleep(1)
        os.remove(f"{Neo4JImportPath}/{filename.split('.')[0]}_{i}.xml")

    with driver.session() as session:
        #NOTE: We are not doing: MATCH (n) WHERE size(keys(properties(n))) < 2 DETACH DELETE n
        #Because we assume the CREATE clase is only invoked if there are non-empty rows
        #NOTE: We purge by merging all products with the same EMA_MA_Number or FDA_Application_Number
        session.write_transaction(misc.remove_duplicate_nodes, "Product", "n.EMA_MA_Number as ema_nb")
        session.write_transaction(misc.remove_duplicate_nodes, "Product", "n.FDA_Application_Number as fda_nb")
        #Also remove duplicates by PubMed ID
        session.write_transaction(misc.remove_duplicate_nodes, "Publication", "n.Pubmed_ID as id", "WHERE n.Pubmed_ID IS NOT null")
        bar()
        #NOTE: We also merge all dosage nodes that are *exactly equal* in all three values
        session.write_transaction(misc.remove_duplicate_nodes, "Product", "n.Form as frm, n.Stength as str, n.Route as rt")
        #And in two (route *should* be duplicated)
        session.write_transaction(misc.remove_duplicate_nodes, "Product", "n.Stength as str, n.Route as rt", "WHERE n.Form IS null")
        session.write_transaction(misc.remove_duplicate_nodes, "Product", "n.Form as frm, n.Route as rt", "WHERE n.Stength IS null")
        #And we delete duplicate relationships
        session.write_transaction(misc.remove_duplicate_relationships)
        #No need for a dedicated "purge database" command this time!
        #(This shouldn't be necessary but just in case)
        session.run("MATCH (d:Dosage) WHERE size(keys(properties(d))) < 1 DETACH DELETE d")
        session.run("MATCH (p:Product) WHERE size(keys(properties(p))) < 1 DETACH DELETE p")
        bar()
    with driver.session() as session:
        session.write_transaction(misc.export_graphml, "graph.graphml")
        bar()

print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
print(f"A copy of the file has been saved in this project's directory")
