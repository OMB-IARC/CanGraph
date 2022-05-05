#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem

# Import subscripts for the program
import create_nodes
import create_relations
import miscelaneous as misc

hmdb_urls = ["http://smpdb.ca/downloads/smpdb_pathways.csv.zip",
             "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip",
             "http://smpdb.ca/downloads/smpdb_proteins.csv.zip",
             "http://smpdb.ca/downloads/sequences/smpdb_protein.fasta.zip",
             "http://smpdb.ca/downloads/sequences/smpdb_gene.fasta.zip"
            ]

instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
driver = GraphDatabase.driver(instance, auth=(user, passwd))

Neo4JImportPath = misc.get_import_path(driver)

with driver.session() as session:
    session.write_transaction(create_nodes.clean_database)
