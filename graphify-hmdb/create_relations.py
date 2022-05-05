#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This is just a collection of functions used by the "main" script

def add_protein_associations(tx, filename):
    """
    Creates "Protein" nodes based on XML files obtained from the HMDB website.
    NOTE: Unlike the "self.add_protein" function, this creates Proteins based on info on the
    "Metabolite" files, not on the "Protein" files themselves. This could mean node duplication, but,
    hopefully, the MERGE by Acession will mean that this duplicates will be catched.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "protein_associations"] AS protein_associations

        MATCH (m:Metabolite {{Accession:accession}})

        WITH protein_associations, m
        UNWIND protein_associations AS protein_association
        WITH protein_association, m
        UNWIND protein_association["_children"] AS my_protein

        WITH
            [X in my_protein._children WHERE X._type = "protein_accession"][0]._text AS protein_accession,
            [X in my_protein._children WHERE X._type = "name"][0]._text AS name,
            [X in my_protein._children WHERE X._type = "uniprot_id"][0]._text AS uniprot_id,
            [X in my_protein._children WHERE X._type = "gene_name"][0]._text AS gene_name,
            [X in my_protein._children WHERE X._type = "protein_type"][0]._text AS protein_type,
            m

        MERGE (p:Protein {{ Acession:protein_accession }})
        ON CREATE SET p.Gene_Name = gene_name, p.Protein_Type = protein_type, p.Uniprot_ID = uniprot_id

        MERGE (m)-[r:ASSOCIATED_WITH]-(p)
        """)

def add_metabolite_associations(tx, filename):
    """
    WARNING HACEER PROTEIN Y METYBOLITE ASSOCIATIONS NO DIRECCIONALES PARA PDOER BORRAR REPES
    name might overwrite byt necessary on reate
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "metabolite_associations"] AS metabolite_associations

        MATCH (p:Protein {{ Acession:accession }})

        WITH metabolite_associations, p
        UNWIND metabolite_associations AS metabolite_association
        WITH metabolite_association, p
        UNWIND metabolite_association["_children"] AS my_metabolite

        WITH
            [X in my_metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in my_metabolite._children WHERE X._type = "name"][0]._text AS name,
            p

        MERGE (m:Metabolite {{ Acession:protein_accession }})
        ON CREATE SET p.Name = name

        MERGE (m)-[r:ASSOCIATED_WITH]-(p)
        """)

def add_metabolite_references(tx, filename):
    """



    TERMINAR ESTO

    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "metabolite_reference"] AS metabolite_reference

        MATCH (p:Protein {{ Acession:accession }})

        WITH metabolite_reference, p
        UNWIND metabolite_reference AS metabolite_ref
        WITH metabolite_ref, p
        UNWIND metabolite_ref["_children"] AS my_refs

        WITH
            [X in my_metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in my_metabolite._children WHERE X._type = "name"][0]._text AS name,
            p

        MERGE (m:Metabolite {{ Acession:protein_accession }})
        ON CREATE SET p.Name = name

        MERGE (m)-[r:ASSOCIATED_WITH]-(p)
        """)




