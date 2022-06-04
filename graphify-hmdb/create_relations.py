#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This is just a collection of functions used by the "main" script

def add_protein_associations(tx, filename):
    """
    Creates "Protein" nodes based on XML files obtained from the HMDB website.
    NOTE: Unlike the "create_nodes.add_protein" function, this creates Proteins based on info on the
    "Metabolite" files, not on the "Protein" files themselves. This could mean node duplication, but,
    hopefully, the MERGE by Accession will mean that this duplicates will be catched.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "metabolite"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "protein_associations"] AS protein_associations

        MERGE (m:Metabolite {{Accession:accession}})

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

        MERGE (p:Protein {{ Accession:protein_accession }})
        ON CREATE SET p.Gene_Name = gene_name, p.Protein_Type = protein_type, p.Uniprot_ID = uniprot_id

        MERGE (m)-[r:ASSOCIATED_WITH]-(p)
        """)

def add_metabolite_associations(tx, filename):
    """
    Adds associations contained in the "protein" file, between proteins and metabolites.
    NOTE: Like he "create_nodes.add_metabolite_associations" function, this creates non-directional
    relationships (m)-[r:ASSOCIATED_WITH]-(p) ; this helps duplicates be detected.
    NOTE: The "ON CREATE SET" clause for the "Name" param ensures no overwriting
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///{filename}")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "metabolite_associations"] AS metabolite_associations

        MERGE (p:Protein {{ Accession:accession }})

        WITH metabolite_associations, p
        UNWIND metabolite_associations AS metabolite_association
        WITH metabolite_association, p
        UNWIND metabolite_association["_children"] AS my_metabolite

        WITH
            [X in my_metabolite._children WHERE X._type = "accession"][0]._text AS metabolite_accession,
            [X in my_metabolite._children WHERE X._type = "name"][0]._text AS name,
            p

        MERGE (m:Metabolite {{ Accession:metabolite_accession }})
        ON CREATE SET m.Name = name

        MERGE (m)-[r:ASSOCIATED_WITH]-(p)
        """)

def add_metabolite_references(tx, filename):
    """
    Creates references for relations betweens Protein nodes and Metabolite nodes
    WARNING: Unfortunately, Neo4J makes it really, really, really difficult to work with XML,
    and so, this time, a r.Pubmed_ID list with the references could not be created. Nonetheless,
    I considered adding this useful.
    """
    return tx.run(f"""
        CALL apoc.load.xml("file:///test2.xml")
        YIELD value
        WITH [x in value._children WHERE x._type = "protein"] AS metabolites
        UNWIND metabolites AS metabolite
        WITH
            [X in metabolite._children WHERE X._type = "accession"][0]._text AS accession,
            [X in metabolite._children WHERE X._type = "metabolite_references"] AS metabolite_references

        MERGE (p:Protein {{ Accession:accession }})

        WITH metabolite_references, p
        UNWIND metabolite_references AS metabolite_reference
        WITH metabolite_reference, p
        UNWIND metabolite_reference["_children"] AS my_reference
        WITH my_reference, p
        UNWIND my_reference["_children"] AS my_ref

        WITH
            [X in my_ref._children WHERE X._type = "accession"][0]._text AS metabolite_accession,
            [X in my_ref._children WHERE X._type = "name"][0]._text AS name,
            [X in my_ref._children WHERE X._type = "reference_text"][0]._text AS reference_text,
            [X in my_ref._children WHERE X._type = "pubmed_id"][0]._text AS pubmed_id,
            p

        FOREACH(ignoreMe IN CASE WHEN metabolite_accession IS NOT null THEN [1] ELSE [] END |
            MERGE (m:Metabolite {{ Accession:metabolite_accession }})
            ON CREATE SET m.name = name
            MERGE (m)-[r:ASSOCIATED_WITH]-(p)
            )

        FOREACH(ignoreMe IN CASE WHEN reference_text IS NOT null THEN [1] ELSE [] END |
            MERGE (pu:Publication {{Authors:split(reference_text, ":")[0]}})
            SET pu.Abstract = split(replace(reference_text, split(reference_text, ":")[0]+": ", ""), ".")[0]
            SET pu.Publication = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[0]
            SET pu.Notes = split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[2]
            SET pu.Date = split(split(replace(reference_text, split(reference_text, ".")[0]+". ",""), ".")[1],";")[0]
            SET pu.Volume = split(split(reference_text, ";")[1], "(")[0]
            SET pu.Number = split(split(reference_text, "(")[1], ")")[0]
            SET pu.Pages = split(split(reference_text, ":")[-1], ".")[0]
            SET pu.Pubmed_ID = pubmed_id
            MERGE (p)-[r1:CITED_IN]->(pu)
            )
        """)




