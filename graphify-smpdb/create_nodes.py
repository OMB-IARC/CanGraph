#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This is just a collection of functions used by the "main" script

def add_metabolites(tx, filename):
    """
    Adds "Metabolite" nodes to the database, according to individual CSVs present in the SMPDB website
    NOTE: Some of the node's properties might be set to "null" (important in order to work with it)
    NOTE: This database clearly differentiates Metabolites and Proteins, so no overlap is accounted for
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
        MERGE (m:Metabolite {{ Metabolite_ID:line["Metabolite ID"] }})
        SET m.Name = line["Metabolite Name"], m.HMDB_ID = line["HMDB ID"], m.KEGG_ID = line["KEGG ID"], m.ChEBI_ID = line["ChEBI ID"],
            m.DrugBank_ID = line["DrugBank ID"], m.CAS = line.CAS, m.Formula = line.Formula, m.IUPAC = line.IUPAC, m.SMILES = line.SMILES,
            m.InChI = line.InChI, m.InChI_Key = line["InChI Key"]
        MERGE (pa:Pathway {{ SMPDB_ID:line["SMPDB ID"] }} )
        SET pa.Name = line["Pathway Name"], pa.Subject = line["Pathway Subject"]
        MERGE (m)-[r:PART_OF_PÀTHWAY]->(pa)
        """)

def add_proteins(tx, filename):
    """
    Adds "Protein" nodes to the database, according to individual CSVs present in the SMPDB website
    NOTE: Some of the node's properties might be set to "null" (important in order to work with it)
    NOTE: This database clearly differentiates Metabolites and Proteins, so no overlap is accounted for
    TODO: Why is the SMPDB_ID property called like that and not SMDB_ID?
    WARNING: Since no unique identifier was found, CREATE had to be used (instead of merge). This might create
             duplicates. which should be accounted for.
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
        CREATE (p:Protein )
        SET p.Name = line["Protein Name"], p.HMDBP_ID = line["HMDBP ID"], p.DrugBank_ID = line["DrugBank ID"],
            p.GenBank_ID = line["GenBank ID"], p.Gene_Name = line["Gene Name"], p.Locus = line["Locus"],
            p.UniProt_ID = line["Uniprot ID"]
        MERGE (pa:Pathway {{ SMPDB_ID:line["SMPDB ID"] }} )
        SET pa.Name = line["Pathway Name"], pa.Subject = line["Pathway Subject"]
        MERGE (p)-[r:PART_OF_PÀTHWAY]->(pa)
        """)

def add_pathways(tx, filename):
    """
    Adds "Pathways" nodes to the database, according to individual CSVs present in the SMPDB website
    Since this is done after the creation of said pathways in the last step, this will most likely just annotate them.
    TODO: This file is really big. It could be divided into smaller ones.
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line
        MERGE (pa:Pathway {{ SMPDB_ID:line["SMPDB ID"] }} )
        SET pa.Name = line["Name"], pa.Subject = line["Subject"], pa.PW_ID = line["PW ID"], pa.Description = line["Description"]
        """)

def add_sequence(tx, seq_id, seq_name, seq_type, seq, seq_format="FASTA"):
    """
    Adds "Pathways" nodes to the database, according to the sequences presented in FASTA files from the SMPDB website
    """
    return tx.run(f"""
        MERGE (s:Sequence {{ UniProt_ID:"{seq_id}", Name:"{seq_name}", Type:"{seq_type}", Sequence:"{seq}", Format:"{seq_format}" }} )
        MERGE (p:Protein {{ UniProt_ID:"{seq_id}" }} )
        MERGE (p)-[r:SEQUENCED_AS]->(s)
        """)









