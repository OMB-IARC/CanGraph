##!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Import modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
import pandas as pd                  # WARNING not used, should be removed.

import create_nodes
import create_relations
import miscelaneous as misc

# WARNING : Hardcoded connection info: MUST REMOVE FROM PRODUCTION
driver = GraphDatabase.driver("neo4j://localhost:7687", auth=("neo4j", "MJwBfFD$a2@Qa7"))

all_urls = ["http://exposome-explorer.iarc.fr/system/downloads/current/biomarkers.csv.zip",
            "http://exposome-explorer.iarc.fr/system/downloads/current/microbial_metabolites.csv.zip",
            "http://exposome-explorer.iarc.fr/system/downloads/current/concentrations.csv.zip",
            "http://exposome-explorer.iarc.fr/system/downloads/current/reproducibilities.csv.zip",
            "http://exposome-explorer.iarc.fr/system/downloads/current/correlations.csv.zip",
            "http://exposome-explorer.iarc.fr/system/downloads/current/metabolomic_associations.csv.zip",
            "http://exposome-explorer.iarc.fr/system/downloads/current/microbial_metabolite_identifications.csv.zip",
            "http://exposome-explorer.iarc.fr/system/downloads/current/cancer_associations.csv.zip",
            "http://exposome-explorer.iarc.fr/system/downloads/current/publications.csv.zip"]

Neo4JImportPath = misc.get_import_path(driver)

for url in all_urls:
    misc.download_and_unzip(url, Neo4JImportPath)

with driver.session() as session:
    session.write_transaction(create_nodes.clean_database)
with driver.session() as session:
    session.write_transaction(create_nodes.biomarkers)
with driver.session() as session:
    session.write_transaction(create_nodes.microbial_metabolites)
with driver.session() as session:
    session.write_transaction(create_nodes.concentrations)
with driver.session() as session:
    session.write_transaction(create_nodes.publications)


with driver.session() as session:
    session.write_transaction(create_relations.correlation)

df = pd.read_csv("../correlations.csv", header = 0)
print(df.columns)
