#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem
from Bio import SeqIO                # To parse FASTA files
from time import sleep               # A hack to avoid starving the system resources

# Import subscripts for the program
import create_nodes
import miscelaneous as misc

smpdb_urls = ["http://smpdb.ca/downloads/smpdb_pathways.csv.zip",
             "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip",
             "http://smpdb.ca/downloads/smpdb_proteins.csv.zip",
             "http://smpdb.ca/downloads/sequences/smpdb_protein.fasta.zip",
             "http://smpdb.ca/downloads/sequences/smpdb_gene.fasta.zip",
             ]

instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
driver = GraphDatabase.driver(instance, auth=(user, passwd))

Neo4JImportPath = misc.get_import_path(driver)

print("Connected to Neo4J")

with driver.session() as session:
    session.run( misc.clean_database() )

print("Cleaned DataBase")

for url in smpdb_urls:

    split = url.split('/')[-1].split('.')[0]
    print(f"Downloading and Unzipping: {split}...")

    if split in ["smpdb_metabolites", "smpdb_proteins"]:
        misc.download_and_unzip(url, f"csvfolder/{split}")
    else:
        misc.download_and_unzip(url, "csvfolder/")

total_files = (len(os.listdir("csvfolder/smpdb_metabolites")) +
               len(os.listdir("csvfolder/smpdb_proteins")) +
               1 +
               len(list(SeqIO.parse(open("csvfolder/smpdb_protein.fasta"),'fasta'))) +
               len(list(SeqIO.parse(open("csvfolder/smpdb_gene.fasta"),'fasta'))) +
               1)

with alive_bar(total_files) as bar:

    print(f"Importing Metabolites...")
    for filename os.listdir("csvfolder/smpdb_metabolites"):
            filepath = os.path.abspath(f"./csvfolder/smpdb_metabolites/{filename}")
            shutil.copyfile(filepath, f"{Neo4JImportPath}/{filename}")
            with driver.session() as session:
                session.write_transaction(create_nodes.add_metabolites, filename)
                bar()
            os.remove(f"{Neo4JImportPath}/{filename}")

    print(f"Importing Proteins...")
    for filename in os.listdir("csvfolder/smpdb_proteins"):
            filepath = os.path.abspath(f"./csvfolder/smpdb_proteins/{filename}")
            shutil.copyfile(filepath, f"{Neo4JImportPath}/{filename}")
            with driver.session() as session:
                session.write_transaction(create_nodes.add_proteins, filename)
                bar()
            os.remove(f"{Neo4JImportPath}/{filename}")

    print(f"Importing Pathways...")
    with driver.session() as session:
        filepath = os.path.abspath(f"./csvfolder/smpdb_pathways.csv")
        shutil.copyfile(filepath, f"{Neo4JImportPath}/smpdb_pathways.csv")
        session.write_transaction(create_nodes.add_pathways, "smpdb_pathways.csv")
        bar()
        os.remove(f"{Neo4JImportPath}/smpdb_pathways.csv")

    print(f"Importing Genomic Sequences...")
    with driver.session() as session:
        genomic_sequences = SeqIO.parse(open("csvfolder/smpdb_gene.fasta"),'fasta')
        for fasta in genomic_sequences:
            name, sequence = fasta.id,9
            #NOTE: The [-1] here is ESSENTIAL, since there are some duplicated parenthesis sometimes
            UniProt_ID = fasta.description.split("(")[-1].split(")")[0]
            Name =  " ".join(fasta.description.split(" ")[1:]).replace("("+UniProt_ID+")", "")
            #NOTE: Could type also be: NUC?
            Sequence = str(fasta.description) + "\n" + str(fasta.seq); Format = "FASTA"; Type = "DNA"
            session.write_transaction(create_nodes.add_sequence, UniProt_ID, Name, Type, Sequence, Format)
            bar()

    print(f"Importing Protein Sequences...")
    with driver.session() as session:
        genomic_sequences = SeqIO.parse(open("csvfolder/smpdb_protein.fasta"),'fasta')
        for fasta in genomic_sequences:
            name, sequence = fasta.id,9
            #NOTE: The [-1] here is ESSENTIAL, since there are some duplicated parenthesis sometimes
            UniProt_ID = fasta.description.split("(")[-1].split(")")[0]
            Name =  " ".join(fasta.description.split(" ")[1:]).replace("("+UniProt_ID+")", "")
            #NOTE: Could type also be: NUC?
            Sequence = str(fasta.description) + "\n" + str(fasta.seq); Format = "FASTA"; Type = "PROT"
            session.write_transaction(create_nodes.add_sequence, UniProt_ID, Name, Type, Sequence, Format)
            bar();

    with driver.session() as session:
        # We purge the protein since there were some null UniProt_IDs
        session.write_transaction(misc.remove_duplicate_nodes, "Protein", "n.UniProt_ID as id", "WHERE n.UniProt_ID IS NOT null")
        session.write_transaction(misc.remove_duplicate_nodes, "Protein", "n.Name as id", "WHERE n.UniProt_ID IS null")
        session.write_transaction(misc.export_graphml, "graph.graphml")
        bar()

print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
print(f"A copy of the file has been saved in this project's directory")
