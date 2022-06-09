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
import create_nodes
import create_relations
# A hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# NOTE: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

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

    all_csvs = os.listdir(sys.argv[4])
    relation_tables = ["cancer_associations.csv", "metabolomic_associations.csv", "correlations.csv"]
    node_tables = [x for x in all_csvs if x not in relation_tables]
    bar()

    for filename in node_tables:
        filepath = os.path.abspath(f"{sys.argv[4]}/{filename}")
        shutil.copyfile(filepath, f"{Neo4JImportPath}/{filename}")

        with driver.session() as session:
            session.write_transaction(create_nodes.import_csv, filename, filename.split(".")[0])
            bar()

    with driver.session() as session:
        shutil.copyfile(f"{sys.argv[4]}/cancer_associations.csv", f"{Neo4JImportPath}/cancer_associations.csv")
        session.write_transaction(create_relations.cancer_associations)
        bar()
    with driver.session() as session:
        shutil.copyfile(f"{sys.argv[4]}/cancer_associations.csv", f"{Neo4JImportPath}/metabolomic_associations.csv")
        session.write_transaction(create_relations.metabolomic_associations)
        bar()
    with driver.session() as session:
        shutil.copyfile(f"{sys.argv[4]}/cancer_associations.csv", f"{Neo4JImportPath}/correlations.csv")
        session.write_transaction(create_relations.correlations)
        bar()
    with driver.session() as session:
        session.write_transaction(create_relations.auto_units)
        bar()
    with driver.session() as session:
        session.write_transaction(create_relations.measurements_stuff)
        bar()
    with driver.session() as session:
        session.write_transaction(create_relations.reproducibilities)
        bar()
    with driver.session() as session:
        session.write_transaction(create_relations.subjects)
        bar()
    with driver.session() as session:
        session.write_transaction(create_relations.samples)
        bar()
    with driver.session() as session:
        session.write_transaction(create_relations.microbial_metabolite_identifications)
        bar()
    with driver.session() as session:
        session.write_transaction(misc.export_graphml, "graph.graphml")
        bar()

print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
print(f"A copy of the file has been saved in this project's directory")
