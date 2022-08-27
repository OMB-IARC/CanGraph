#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that provides the necessary functions to transition the SMPDB database to graph format,
either from scratch importing all the nodes (as showcased in :obj:`CanGraph.GraphifySMPDB.main`) or in a case-by-case basis,
to annotate existing metabolites (as showcased in :obj:`CanGraph.main`).
"""

# Import external modules necessary for the script
import os, sys, shutil               # Vital modules to interact with the filesystem
from Bio import SeqIO                # To parse FASTA files
import pandas as pd                  # Process tabular data

# Import subscripts for the program
# This hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# .. NOTE:: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

def add_metabolites(tx, filename):
    """
    Adds "Metabolite" nodes to the database, according to individual CSVs present in the SMPDB website

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Some of the node's properties might be set to "null" (important in order to work with it)

    .. NOTE:: This database clearly differentiates Metabolites and Proteins, so no overlap is accounted for
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
        MERGE (m:Metabolite {{ Metabolite_ID:line["Metabolite ID"] }})
        SET m.Name = line["Metabolite Name"], m.HMDB_ID = line["HMDB ID"], m.KEGG_ID = line["KEGG ID"], m.ChEBI_ID = line["ChEBI ID"],
            m.DrugBank_ID = line["DrugBank ID"], m.CAS_Number = line.CAS, m.Formula = line.Formula, m.IUPAC = line.IUPAC, m.SMILES = line.SMILES,
            m.InChI = line.InChI, m.InChIKey = line["InChI Key"]
        MERGE (pa:Pathway {{ SMPDB_ID:line["SMPDB ID"] }} )
        SET pa.Name = line["Pathway Name"], pa.Subject = line["Pathway Subject"]
        MERGE (m)-[r:PART_OF_PATHWAY]->(pa)
        """)

def add_proteins(tx, filename):
    """
    Adds "Protein" nodes to the database, according to individual CSVs present in the SMPDB website

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: Some of the node's properties might be set to "null" (important in order to work with it)

    .. NOTE:: This database clearly differentiates Metabolites and Proteins, so no overlap is accounted for

    .. TODO:: Why is the SMPDB_ID property called like that and not SMDB_ID?

    .. WARNING:: Since no unique identifier was found, CREATE had to be used (instead of merge). This might create
             duplicates. which should be accounted for.
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
        CREATE (p:Protein )
        SET p.Name = line["Protein Name"], p.HMDB_ID = line["HMDBP ID"], p.DrugBank_ID = line["DrugBank ID"],
            p.Genbank_Protein_ID = line["GenBank ID"], p.Gene_Name = line["Gene Name"], p.Locus = line["Locus"],
            p.UniProt_ID = line["Uniprot ID"]
        MERGE (pa:Pathway {{ SMPDB_ID:line["SMPDB ID"] }} )
        SET pa.Name = line["Pathway Name"], pa.Subject = line["Pathway Subject"]
        MERGE (p)-[r:PART_OF_PATHWAY]->(pa)
        """)

def add_pathways(tx, filename):
    """
    Adds "Pathways" nodes to the database, according to individual CSVs present in the SMPDB website
    Since this is done after the creation of said pathways in the last step, this will most likely just annotate them.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. TODO:: This file is really big. It could be divided into smaller ones.
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
        MERGE (pa:Pathway {{ SMPDB_ID:line["SMPDB ID"] }} )
        SET pa.Name = line["Name"], pa.Category = line["Subject"], pa.PW_ID = line["PW ID"], pa.Description = line["Description"]
        """)

def add_sequence(tx, seq_id, seq_name, seq_type, seq, seq_format="FASTA"):
    """
    Adds "Pathways" nodes to the database, according to the sequences presented in FASTA files from the SMPDB website

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        seq_id      (str): The UniProt Database Identifier for the sequence that is been imported
        seq_name    (str): The Name (i.e. FASTA header) of the Sequence that is been imported
        seq_type    (bool): The type of the sequence; can be either of ["DNA", "PROT"]
        seq         (str): The seuqnce that is been imported; a text chain of nucleotides or aminoacids, identified by their acronyms
        seq_format  (str): The format the sequence is provided under; default is "FASTA", but its optional

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.
    """
    return tx.run(f"""
        MERGE (s:Sequence {{ UniProt_ID:"{seq_id}", Name:"{seq_name}", Type:"{seq_type}", Sequence:"{seq}", Format:"{seq_format}" }} )
        MERGE (p:Protein {{ UniProt_ID:"{seq_id}" }} )
        MERGE (p)-[r:SEQUENCED_AS]->(s)
        """)

def build_from_file(databasepath, filepath, Neo4JImportPath, driver, filetype):
    """
    A function able to build a portion of the SMPDB in graph format, provided that one CSV is supplied to it.
    This CSVs are downloaded from the website, and can be presented either as the full file, or as a splitted
    version of it, with just one item per file (which is recommended due to memory limitations)

    Since file title represents a different pathway, the function automatically picks up and import the relative pathway node.

    Args:
        databasepath (str): The path to the database where all SMPDB CSVs are stored
        filepath (str): The path to the current file being imported
        Neo4JImportPath (str): The path from which Neo4J is importing data
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        filetype (bool): The type of file being imported; one of ether ["Metabolite", "Protein"]- If the file is a FASTA sequence store, this will be auto-detected.

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.

     .. NOTE:: Since this adds a ton of low-resolution nodes, maybe have this db run first?
    """
    shutil.copyfile(os.path.abspath(filepath), f"{Neo4JImportPath}/{os.path.basename(filepath)}")
    pathway_id = filepath.split("/")[-1].split("_")[0]
    all_the_proteins = pd.read_csv(f"{os.path.abspath(databasepath)}/SMPDB/smpdb_proteins/{pathway_id}_proteins.csv")

    if filetype == "Metabolite":
        with driver.session() as session:
            session.write_transaction(add_metabolites, f"{os.path.basename(filepath)}")
            shutil.copyfile(f"{os.path.abspath(databasepath)}/SMPDB/smpdb_proteins/{pathway_id}_proteins.csv", f"{Neo4JImportPath}/corresponding.csv")
            session.write_transaction(add_proteins, "corresponding.csv")
    elif filetype == "Protein":
        with driver.session() as session:
            session.write_transaction(add_proteins, f"{os.path.basename(filepath)}")
            shutil.copyfile(f"{os.path.abspath(databasepath)}/SMPDB/smpdb_metabolites/{pathway_id}_metabolites.csv", f"{Neo4JImportPath}/corresponding.csv")
            session.write_transaction(add_metabolites, "corresponding.csv")

    # We then import pathway info
    all_pathways = pd.read_csv(f"{os.path.abspath(databasepath)}/SMPDB/smpdb_pathways.csv", delimiter=',', header=0)
    just_this_pathway = all_pathways.loc[all_pathways["SMPDB ID"] == f"{pathway_id}"]
    just_this_pathway.to_csv(f"{Neo4JImportPath}/just_this_pathway.csv", index=False)


    # And, finally, we add the protein sequences (They only reference proteins, not metabolites!):
    genomic_sequences = SeqIO.parse(open(f"{os.path.abspath(databasepath)}/SMPDB/smpdb_gene.fasta"),'fasta')

    for fasta in genomic_sequences:
        if fasta.description.split("(")[-1].split(")")[0] in all_the_proteins["Uniprot ID"].values:
            UniProt_ID = fasta.description.split("(")[-1].split(")")[0]
            Name =  " ".join(fasta.description.split(" ")[1:]).replace("("+UniProt_ID+")", "")
            Sequence = str(fasta.description) + "\n" + str(fasta.seq); Format = "FASTA"; Type = "DNA"
            with driver.session() as session:
                session.write_transaction(add_sequence, UniProt_ID, Name, Type, Sequence, Format)

    proteic_sequences = SeqIO.parse(open(f"{os.path.abspath(databasepath)}/SMPDB/smpdb_protein.fasta"),'fasta')
    for fasta in proteic_sequences:
        if fasta.description.split("(")[-1].split(")")[0] in all_the_proteins["Uniprot ID"].values:
            UniProt_ID = fasta.description.split("(")[-1].split(")")[0]
            Name =  " ".join(fasta.description.split(" ")[1:]).replace("("+UniProt_ID+")", "")
            Sequence = str(fasta.description) + "\n" + str(fasta.seq); Format = "FASTA"; Type = "PROT"
            with driver.session() as session:
                session.write_transaction(add_sequence, UniProt_ID, Name, Type, Sequence, Format)
    os.remove(f"{Neo4JImportPath}/just_this_pathway.csv"); os.remove(f"{Neo4JImportPath}/corresponding.csv"); os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")
