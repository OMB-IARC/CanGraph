#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that leverages the functions present in the :obj:`~CanGraph.GraphifySMPDB.build\_database`
module to recreate `the SMPDB database <http://smpdb.ca/>`_ using a graph format and Neo4J,
and then provides an GraphML export file.

Please note that, to work, the functions here pre-suppose you have internet access, which will be used to download
HMDB's CSVs under ```./csvfolder/``` (please ensure you have read-write access there).

For more details on how to run this script, please consult the package's README
"""

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem
from Bio import SeqIO                # To parse FASTA files
from time import sleep               # A hack to avoid starving the system resources

# Import subscripts for the program
import build_database
# A hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# NOTE: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc


def import_metabolites(filename, Neo4JImportPath):
    """
    Imports "Metabolite" files from the SMPDB Database.
    The function assumes ```filename``` is availaible at ```./csvfolder/smpdb_metabolites/```

    Args:
        filename        (str): The name of the CSV file that is being imported
        Neo4JImportPath (str): The path from which Neo4J is importing data

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    print(f"Importing Metabolites...")
    for filename in os.listdir("csvfolder/smpdb_metabolites"):
            filepath = os.path.abspath(f"./csvfolder/smpdb_metabolites/{filename}")
            shutil.copyfile(filepath, f"{Neo4JImportPath}/{filename}")
            with driver.session() as session:
                misc.manage_transaction(build_database.add_metabolites, filename)
                bar()
            os.remove(f"{Neo4JImportPath}/{filename}")

def import_proteins(filename, Neo4JImportPath):
    """
    Imports "Protein" files from the SMPDB Database.
    The function assumes ```filename``` is availaible at ```./csvfolder/smpdb_proteins/```

    Args:
        filename        (str): The name of the CSV file that is being imported
        Neo4JImportPath (str): The path from which Neo4J is importing data

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    print(f"Importing Proteins...")
    for filename in os.listdir("csvfolder/smpdb_proteins"):
            filepath = os.path.abspath(f"./csvfolder/smpdb_proteins/{filename}")
            shutil.copyfile(filepath, f"{Neo4JImportPath}/{filename}")
            with driver.session() as session:
                misc.manage_transaction(build_database.add_proteins, filename)
                bar()
            os.remove(f"{Neo4JImportPath}/{filename}")

def import_pathways(Neo4JImportPath):
    """
    Imports the ```smpdb_pathways.csv``` file from the SMPDB Database.
    The function assumes the file is availaible at ```./csvfolder/smpdb_gene.fasta```

    Args:
        Neo4JImportPath (str): The path from which Neo4J is importing data

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    print(f"Importing Pathways...")
    with driver.session() as session:
        filepath = os.path.abspath(f"./csvfolder/smpdb_pathways.csv")
        shutil.copyfile(filepath, f"{Neo4JImportPath}/smpdb_pathways.csv")
        misc.manage_transaction(build_database.add_pathways, "smpdb_pathways.csv")
        bar()
        os.remove(f"{Neo4JImportPath}/smpdb_pathways.csv")

def import_genomic_seqs():
    """
    Imports the ```smpdb_gene.fasta``` file from the SMPDB Database.
    The function assumes the file is availaible at ```./csvfolder/smpdb_gene.fasta```

    Args:
        Neo4JImportPath (str): The path from which Neo4J is importing data

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    print(f"Importing Genomic Sequences...")
    with driver.session() as session:
        genomic_sequences = SeqIO.parse(open("./csvfolder/smpdb_gene.fasta"),'fasta')
        for fasta in genomic_sequences:
            #NOTE: The [-1] here is ESSENTIAL, since there are some duplicated parenthesis sometimes
            UniProt_ID = fasta.description.split("(")[-1].split(")")[0]
            Name =  " ".join(fasta.description.split(" ")[1:]).replace("("+UniProt_ID+")", "")
            #NOTE: Could type also be: NUC?
            Sequence = str(fasta.description) + "\n" + str(fasta.seq); Format = "FASTA"; Type = "DNA"
            misc.manage_transaction(build_database.add_sequence, UniProt_ID, Name, Type, Sequence, Format)
            bar()

def import_proteic_seqs():
    """
    Imports the ```smpdb_protein.fasta``` file from the SMPDB Database.
    The function assumes the file is availaible at ```./csvfolder/smpdb_protein.fasta```

    Args:
        Neo4JImportPath (str): The path from which Neo4J is importing data

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    print(f"Importing Protein Sequences...")
    with driver.session() as session:
        proteic_sequences = SeqIO.parse(open("./csvfolder/smpdb_protein.fasta"),'fasta')
        for fasta in proteic_sequences:
            #NOTE: The [-1] here is ESSENTIAL, since there are some duplicated parenthesis sometimes
            UniProt_ID = fasta.description.split("(")[-1].split(")")[0]
            Name =  " ".join(fasta.description.split(" ")[1:]).replace("("+UniProt_ID+")", "")
            #NOTE: Could type also be: NUC?
            Sequence = str(fasta.description) + "\n" + str(fasta.seq); Format = "FASTA"; Type = "PROT"
            misc.manage_transaction(build_database.add_sequence, UniProt_ID, Name, Type, Sequence, Format)
            bar();

def main():
    """
    The function that executes the code
    """

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
            misc.download_and_unzip(url, f"./csvfolder/{split}")
        else:
            misc.download_and_unzip(url, "./csvfolder/")

    total_files = (len(os.listdir("csvfolder/smpdb_metabolites")) +
                len(os.listdir("csvfolder/smpdb_proteins")) +
                1 +
                len(list(SeqIO.parse(open("csvfolder/smpdb_protein.fasta"),'fasta'))) +
                len(list(SeqIO.parse(open("csvfolder/smpdb_gene.fasta"),'fasta'))) +
                1)

    with alive_bar(total_files) as bar:
        import_metabolites()
        import_proteins()
        import_pathways()
        import_genomic_seqs()
        import_proteic_seqs()

    # At the end, we purge the database
    misc.purge_database(driver)

    # And export it:
    with driver.session() as session:
        misc.manage_transaction(misc.export_graphml, "graph.graphml")
        bar()

    print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
    shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
    print(f"A copy of the file has been saved in this project's work directory")


if __name__ == '__main__':

    main()
