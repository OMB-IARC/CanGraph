#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This is just a collection of functions used by the "main" script

def clean_database(tx):
    """ Gets all the nodes in a Neo4J database and removes them """
    return tx.run("""
        MATCH (n) DETACH DELETE n;
        """)

def initial_cancer_discovery(tx, number=None):
    """
    A Neo4J Cypher Statment that queries wikidata for Human Cancers. Since using the "afflicts:human"
    tag didnt have much use here, I used a simple workaround: Query wikidata for all humans, and, among them,
    find all of this for which their cause of death was a subclass of "Cancer" (Q12078). Unfortunaltely,
    some of them were diagnosed "Cancer" (Q12078), which is too general, so I removed it.
    """
    return tx.run("""
        WITH 'SELECT DISTINCT ?cause_of_death ?cause_of_death_name
        WHERE {
                ?human wdt:P31 wd:Q5;
                wdt:P509 ?cause_of_death.
                ?cause_of_death wdt:P279* wd:Q12078.
        SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                 ?cause_of_death rdfs:label ?cause_of_death_name. }
                }' AS sparql

        CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                { Accept: "application/sparql-results+json"}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row
        WITH row['cause_of_death_name']['value'] as label,
            row['cause_of_death']['value'] as url,
            split(row['cause_of_death']['value'],'/')[-1] as id

        CREATE (c:Disease)
        SET c.Name = label,
            c.URL = url,
            c.WikiData_ID = id

        WITH c
        MATCH (n:Disease)
            WHERE n.WikiData_ID = "Q12078"
        DETACH DELETE n
        """)

def find_subclass_of_cancer(tx, number=None):
    """
    A Neo4J Cypher Statment that queries wikidata for subclasses of "Cancer" nodes already
    present on the Database. Since this are expected to only affect humans, this subclasses
    should also, only affect humans
    """
    return tx.run("""
        MATCH (c:Disease)
        WITH 'SELECT DISTINCT ?cancer ?cancer_name
        WHERE {
                ?cancer wdt:P279 wd:' + c.WikiData_ID + ' .
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                 ?cancer rdfs:label ?cancer_name. }
                }' AS sparql, c

        CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                { Accept: "application/sparql-results+json"}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        CREATE (s:Disease)
        SET s.Name = row['cancer_name']['value'],
            s.URL = row['cancer']['value'],
            s.WikiData_ID = split(row['cancer']['value'],'/')[-1]

        CREATE (s)-[:SUBCLASS_OF]->(c)
        """)

def find_instance_of_cancer(tx, number=None):
    """
    A Neo4J Cypher Statment that queries wikidata for instances of "Cancer" nodes already
    present on the Database. Since this are expected to only affect humans, this subclasses
    should also, only affect humans
    """
    return tx.run("""
        MATCH (c:Disease)
        WITH 'SELECT DISTINCT ?cancer ?cancer_name
        WHERE {
                ?cancer wdt:P31 wd:' + c.WikiData_ID + ' .
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                 ?cancer rdfs:label ?cancer_name. }
                }' AS sparql, c

        CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                { Accept: "application/sparql-results+json"}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        CREATE (i:Disease)
        SET i.Name = row['cancer_name']['value'],
            i.URL = row['cancer']['value'],
            i.WikiData_ID = split(row['cancer']['value'],'/')[-1]

        CREATE (i)-[:INSTANCE_OF]->(c)
        """)

def add_cancer_info(tx, number):
    """
    Adds info to "Cancer" nodes for which its WikiData_ID ends in a given numbe. This way, only some of the nodes
    are targeted, and the Java Machine does not run out of memory
    """
    return tx.run(f"""
        MATCH (c:Disease)
            WHERE toFloat(split(c.WikiData_ID,"")[-1]) = toFloat({number})
        WITH 'SELECT *
            WHERE{{
                filter (?item = wd:' + c.WikiData_ID + ')

            OPTIONAL{{
                ?item wdt:P665 ?kegg_id }}
            OPTIONAL{{
                ?item wdt:P699 ?Disease_Ont }}
            OPTIONAL{{
                ?item wdt:P927 ?Anatomical_Location
                SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                         ?Anatomical_Location rdfs:label ?Anatomical_Location_name. }} }}
            OPTIONAL{{
                ?item wdt:P3841 ?HPO_ID }}
            OPTIONAL{{
                ?item wdt:P494 ?ICD_1O }}
            OPTIONAL{{
                ?item wdt:P7807 ?ICD_11 }}
            OPTIONAL{{
                ?item wdt:P492 ?OMIM_ID }}
            }}' AS sparql, c

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET c.KEGG_ID = row['kegg_id']['value'],
            c.Disease_Ontology_ID = row['Disease_Ont']['value'],
            c.Anatomical_Location = row['Anatomical_Location_name']['value'],
            c.HPO_ID = row['HPO_ID']['value'],
            c.ICD_1O = row['ICD_1O']['value'],
            c.OMIM_ID = row['OMIM_ID']['value'],
            c.ICD_11 = row['ICD_11']['value']
        """)

def add_drugs(tx, number):
    """
    Creates drug nodes related with each of the "Cancer" nodes already on the database
    """
    return tx.run(f"""
        MATCH (c:Disease)
            WHERE toFloat(split(c.WikiData_ID,"")[-1]) = toFloat({number})
        WITH 'SELECT *
            WHERE{{
                filter (?item = wd:' + c.WikiData_ID + ')

            OPTIONAL{{
                ?item wdt:P2176 ?Drugs
                SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                         ?Drugs rdfs:label ?DrugName. }} }}
            OPTIONAL{{
                ?item wdt:P2888 ?Exact_Matches }}
            }}' AS sparql, c

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        FOREACH(ignoreme in case when row['Drugs'] is not null then [1] else [] end |
                MERGE (d:Drug{{URL:row['Drugs']['value'],
                               WikiData_ID:split(row['Drugs']['value'],'/')[-1],
                               Name:row['DrugName']['value'] }} )
                MERGE (c)-[:TARGETED_BY]->(d))

        FOREACH(ignoreme in case when row['Exact_Matches'] is not null then [1] else [] end |
                MERGE (e:ExternalEquivalent{{URL:row['Exact_Matches']['value'],
                                            WikiData_ID:split(row['Exact_Matches']['value'],'/')[-1]}})
                MERGE (c)-[:EQUALS]-(e))
        """)

def add_causes(tx, number):
    """
    Creates drug nodes related with each of the "Cancer" nodes already on the database
    """
    return tx.run(f"""
        MATCH (c:Disease)
            WHERE toFloat(split(c.WikiData_ID,"")[-1]) = toFloat({number})
        WITH 'SELECT ?Causes ?CauseName
            WHERE{{
                filter (?item = wd:' + c.WikiData_ID + ')

            OPTIONAL{{
                ?item wdt:P828 ?Causes
                SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                         ?Causes rdfs:label ?CauseName. }} }}
            }}' AS sparql, c

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        FOREACH(ignoreme in case when row['Causes'] is not null then [1] else [] end |
                MERGE (ca:Cause{{URL:row['Causes']['value'],
                                WikiData_ID:split(row['Causes']['value'],'/')[-1],
                                Name:row['CauseName']['value'] }} )
                MERGE (c)-[:CAUSED_BY]->(ca))
        """)

def add_genes(tx, number):
    """
    Creates gene nodes related with each of the "Cancer" nodes already on the database
    """
    return tx.run(f"""
        MATCH (c:Disease)
            WHERE toFloat(split(c.WikiData_ID,"")[-1]) = toFloat({number})
        WITH 'SELECT *
            WHERE{{
                filter (?item = wd:' + c.WikiData_ID + ')

            OPTIONAL{{
                ?item wdt:P2293 ?Genetic_Associations.
                ?Genetic_Associations wdt:P703 wd:Q15978631
                SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                         ?Genetic_Associations rdfs:label ?GeneName. }} }}
            }}' AS sparql, c

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        FOREACH(ignoreme in case when row['Genetic_Associations'] is not null and row['GeneName'] is not null then [1] else [] end |
                MERGE (g:Gene{{URL:row['Genetic_Associations']['value'],
                              WikiData_ID:split(row['Genetic_Associations']['value'],'/')[-1],
                              Name:row['GeneName']['value'] }})
                MERGE (g)-[:ASSOCIATED_CANCER_GENE]->(c))
        """)

def add_drug_external_ids(tx, number=None):
    """
    Adds some external IDs to any "Drug" nodes already present on the database.
    Since the PDB information had too much values which caused triple duplicates that overcharged the system,
    they were intentionally left out.
    """
    return tx.run("""
        MATCH (d:Drug)
        WITH 'SELECT *
            WHERE{
                filter (?item = wd:' + d.WikiData_ID + ')

            OPTIONAL{
                ?item wdt:P592 ?CHEMBL_ID }
            OPTIONAL{
                ?item wdt:P683 ?ChEBI_ID }
            OPTIONAL{
                ?item wdt:P662 ?PubChem_ID }
            OPTIONAL{
                ?item wdt:P665 ?KEGG_ID }
            OPTIONAL{
                ?item wdt:P715 ?DrugBank }
            OPTIONAL{
                ?item wdt:P231 ?CAS_Number }
            OPTIONAL{
                ?item wdt:P661 ?ChemSpider_ID }
            OPTIONAL{
                ?item wdt:P267 ?ATC_Code }
            OPTIONAL{
                ?item wdt:P486 ?MeSH }
            OPTIONAL{
                ?item wdt:P234 ?InChI }
            OPTIONAL{
                ?item wdt:P235 ?InChIKey }

            OPTIONAL{
                ?item wdt:P2067 ?Mass }
            OPTIONAL{
                ?item wdt:P274 ?Chemical_Formula }
            OPTIONAL{
                ?item wdt:P233 ?Canonical_Smiles }

            OPTIONAL{
                ?item wdt:P3489 ?Pregnancy_Category
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                         ?Pregnancy_Category rdfs:label ?Pregnancy_Category_Name. } }
            OPTIONAL{
                ?item wdt:P7830 ?LiverTox }

            }' AS sparql, d

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                { Accept: "application/sparql-results+json"}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET d.CHEMBL_ID = row['CHEMBL_ID']['value'],
            d.ChEBI_ID = row['ChEBI_ID']['value'],
            d.PubChem_ID = row['PubChem_ID']['value'],
            d.KEGG_ID = row['KEGG_ID']['value'],
            d.DrugBank_ID = row['DrugBank']['value'],
            d.CAS_Number = row['CAS_Number']['value'],
            d.ChemSpider_ID = row['ChemSpider_ID']['value'],
            d.InChI = row['InChI']['value'],
            d.InChIKey = row['InChIKey']['value'],
            d.Average_Mass = row['Mass']['value'],
            d.Formula = row['Chemical_Formula']['value'],
            d.SMILES = row['Canonical_Smiles']['value'],
            d.Pregnancy_Category = row['Pregnancy_Category_Name']['value'],
            d.LiverTox = row['LiverTox']['value']

        MERGE (pri:ATC { Code:ATC_Code })
            MERGE (d)-[r:RELATED_ATC]->(pri)

        MERGE (c:MeSH { MeSH_ID:row['MeSH']['value'] })
        MERGE (d)-[r:RELATED_MESH]->(c)
        """)

def add_more_drug_info(tx, number=None):
    """
    Creates some nodes that are related with each of the "Drug" nodes already existing
    on the database: routes of administration, targeted metabolites and approved drugs
    that tehy are been used in
    TODO: ADD ROLE to metabolite interactions
    NOTE: This transaction has been separated in order to keep response times low
    """
    return tx.run("""
        MATCH (d:Drug)
        WITH 'SELECT *
            WHERE{
                filter (?item = wd:' + d.WikiData_ID + ')

                OPTIONAL{
                    ?item wdt:P636 ?Route_of_Admin
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                            ?Route_of_Admin rdfs:label ?Route_of_Admin_name. } }
                OPTIONAL{
                    ?item wdt:P3780 ?Active_ingredient_in
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                            ?Active_ingredient_in rdfs:label ?Active_ingredient_in_name. } }
                OPTIONAL{
                    ?item wdt:P129 ?Interacts_with
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                            ?Interacts_with rdfs:label ?Interacts_with_name. } }
            }' AS sparql, d

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                { Accept: "application/sparql-results+json"}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        FOREACH(ignoreme in case when row['Interacts_with'] is not null then [1] else [] end |
                    MERGE (m:Metabolite{URL:row['Interacts_with']['value'],
                                WikiData_ID:split(row['Interacts_with']['value'],'/')[-1],
                                Name:row['Interacts_with_name']['value'] })
                    MERGE (d)-[:INTERACTS_WITH]-(m))

        FOREACH(ignoreme in case when row['Route_of_Admin'] is not null then [1] else [] end |
                    MERGE (r:AdministrationRoute{URL:row['Route_of_Admin']['value'],
                                WikiData_ID:split(row['Route_of_Admin']['value'],'/')[-1],
                                Name:row['Route_of_Admin_name']['value'] })
                    MERGE (d)-[:ADMINISTERED_VIA]->(r))

        FOREACH(ignoreme in case when row['Active_ingredient_in'] is not null then [1] else [] end |
                    MERGE (me:Product{URL:row['Active_ingredient_in']['value'],
                                WikiData_ID:split(row['Active_ingredient_in']['value'],'/')[-1],
                                Name:row['Active_ingredient_in_name']['value'] })
                    MERGE (d)-[:PART_OF_PRODUCT]->(me))
        """)

def add_gene_info(tx, number=None):
    """
    A Cypher Query that adds some external IDs and properties to "Gene" nodes already existing on
    the database. This query forces the genes to have a "found_in_taxon:homo_sapiens" label. This means
    that any non-human genes will not be annotated (TODO: delete those)
        * Genomic Start and ends keep just the 2nd position, as reported in wikidata
        * TODO: Might include P684 "Orthologues" for more info (it crashed java)
    """
    return tx.run("""
        MATCH (g:Gene)
        WITH 'SELECT *
            WHERE{
                filter (?item = wd:' + g.WikiData_ID + ')
                ?item wdt:P703 wd:Q15978631

            OPTIONAL{
                ?item wdt:P351 ?Entrez_ID }
            OPTIONAL{
                ?item wdt:P594 ?Ensembl_Gene_ID }
            OPTIONAL{
                ?item wdt:P354 ?HGNC_ID }
            OPTIONAL{
                ?item wdt:P492 ?OMIM_ID }
            OPTIONAL{
                ?item wdt:P1057 ?Chromosome_Location
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                        ?Chromosome_Location rdfs:label ?Chromosome_Location_name. } }
            OPTIONAL{
                ?item wdt:P2548 ?Strand_Orientation
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                        ?Strand_Orientation rdfs:label ?Strand_Orientation_name. } }
            OPTIONAL{
                ?item wdt:P644 ?Genomic_Start }
            OPTIONAL{
                ?item wdt:P645 ?Genomic_End }
            OPTIONAL{
                ?item wdt:P4196 ?Cytogenetic_Location  }
            OPTIONAL{
                ?item wdt:P593 ?HomoloGene_ID }

            OPTIONAL{
                ?item wdt:P688 ?Encodes
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                        ?Encodes rdfs:label ?Encodes_name. } }
            OPTIONAL{
                ?item wdt:P5572 ?Expressed_in
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                        ?Expressed_in rdfs:label ?Expressed_in_name. } }
            OPTIONAL{
                ?item wdt:P2888 ?Exact_Matches }

            }' AS sparql, g

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                { Accept: "application/sparql-results+json"}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET g.Entrez_ID = row['Entrez_ID']['value'],
            g.Ensembl_Gene_ID = row['Ensembl_Gene_ID']['value'],
            g.HGNC_ID = row['HGNC_ID']['value'],
            g.OMIM_ID = row['OMIM_ID']['value'],
            g.HomoloGene_ID = row['HomoloGene_ID']['value'],
            g.Cytogenetic_Location = row['Cytogenetic_Location']['value'],
            g.Chromosome_Location = row['Chromosome_Location_name']['value'],
            g.Strand_Orientation = row['Strand_Orientation_name']['value'],
            g.Genomic_Start = row['Genomic_Start']['value'],
            g.Genomic_End = row['Genomic_End']['value']

        FOREACH(ignoreme in case when row['Encodes'] is not null then [1] else [] end |
                    MERGE (m:Metabolite{URL:row['Encodes']['value'],
                                WikiData_ID:split(row['Encodes']['value'],'/')[-1],
                                Name:row['Encodes_name']['value'] })
                    MERGE (g)-[:ENCODES]->(m))

        FOREACH(ignoreme in case when row['Expressed_in'] is not null then [1] else [] end |
                    MERGE (t:Tissue{URL:row['Expressed_in']['value'],
                                WikiData_ID:split(row['Expressed_in']['value'],'/')[-1],
                                Name:row['Expressed_in_name']['value'] })
                    MERGE (g)-[:EXPRESSED_IN]->(t))

        FOREACH(ignoreme in case when row['Exact_Matches'] is not null then [1] else [] end |
                MERGE (e:ExternalEquivalent{URL:row['Exact_Matches']['value'],
                                            WikiData_ID:split(row['Exact_Matches']['value'],'/')[-1]})
                MERGE (g)-[:EQUALS]-(e))
        """)

def add_metabolite_info(tx, number=None):
    """
    A Cypher Query that adds some external IDs and properties to "Metabolite" nodes already existing on
    the database. Two kind of metabolites exist: those that are encoded by a given gene, and those that interact
    with a given drug. Both are adressed here, since they are similar, and, most likely, instances of proteins.
        * This function forces all metabolites to have a "found_in_taxon:human" target
        * The metabolites are not forced to be proteins, but if they are, this is kept in the "instance_of" record
        * TODO: Might include P527 "has part or parts" for more info (it crashed java)
    """
    return tx.run("""
        MATCH (m:Metabolite)
        WITH 'SELECT *
            WHERE{
                filter (?item = wd:' + m.WikiData_ID + ')
                ?item wdt:P703 wd:Q15978631

            OPTIONAL{
                ?item wdt:P31 ?Instance_of
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                            ?Instance_of rdfs:label ?Instance_of_name. } }
            OPTIONAL{
                ?item wdt:P352 ?UniProt_ID }
            OPTIONAL{
                ?item wdt:P486 ?MeSH }
            OPTIONAL{
                ?item wdt:P7260 ?Transporter_DB_ID }
            OPTIONAL{
                ?item wdt:P683 ?ChEBI_ID }
            OPTIONAL{
                ?item wdt:P231 ?CAS_Number }
            OPTIONAL{
                ?item wdt:P274 ?Chemical_Formula }
            OPTIONAL{
                ?item wdt:P234 ?InChI }
            OPTIONAL{
                ?item wdt:P235 ?InChIKey }
            OPTIONAL{
                ?item wdt:P592 ?CHEMBL_ID }
            OPTIONAL{
                ?item wdt:P683 ?ChEBI_ID }
            OPTIONAL{
                ?item wdt:P662 ?PubChem_ID }
            OPTIONAL{
                ?item wdt:P233 ?Canonical_Smiles }
            OPTIONAL{
                ?item wdt:P8117 ?FooDB_Compound_ID }

            OPTIONAL{
                    ?item wdt:P129 ?Interacts_with
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                            ?Interacts_with rdfs:label ?Interacts_with_name. } }
            OPTIONAL{
                ?item wdt:P2888 ?Exact_Matches }

            }' AS sparql, m

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                { Accept: "application/sparql-results+json"}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET m.Function = row['Instance_of_name']['value'],
            m.UniProt_ID = row['UniProt_ID']['value'],
            m.ChEBI_ID = row['ChEBI_ID']['value'],
            m.CAS_Number = row['CAS_Number']['value'],
            m.Formula = row['Chemical_Formula']['value'],
            m.InChI = row['InChI']['value'],
            m.InChIKey = row['InChIKey']['value'],
            m.CHEMBL_ID = row['CHEMBL_ID']['value'],
            m.ChEBI_ID = row['ChEBI_ID']['value'],
            m.PubChem_ID = row['PubChem_ID']['value'],
            m.SMILES = row['Canonical_Smiles']['value'],
            m.FooDB_Compound_ID = row['FooDB_Compound_ID']['value'],
            m.Transporter_DB_ID = row['Transporter_DB_ID']['value']

        MERGE (c:MeSH { MeSH_ID:row['MeSH']['value'] })
        MERGE (m)-[r:RELATED_MESH]->(c)

        FOREACH(ignoreme in case when row['Interacts_with'] is not null then [1] else [] end |
                    MERGE (me:Metabolite{URL:row['Interacts_with']['value'],
                                WikiData_ID:split(row['Interacts_with']['value'],'/')[-1],
                                Name:row['Interacts_with_name']['value'] })
                    MERGE (m)-[:INTERACTS_WITH]-(me))

        FOREACH(ignoreme in case when row['Exact_Matches'] is not null then [1] else [] end |
                MERGE (e:ExternalEquivalent{URL:row['Exact_Matches']['value'],
                                            WikiData_ID:split(row['Exact_Matches']['value'],'/')[-1]})
                MERGE (m)-[:EQUALS]-(e))
        """)

def add_toomuch_metabolite_info(tx, number=None):
    """
    A function that adds loads of info to existing "Metabolite" nodes. This was left out, first because
    it might be too much information, (specially when it is already availaible by clicking the "url" field),
    and because, due to it been so much, it crashes the JVM.
    """
    return tx.run("""
        MATCH (m:Metabolite)
        WITH 'SELECT *
            WHERE{
                filter (?item = wd:' + m.WikiData_ID + ')
                ?item wdt:P703 wd:Q15978631

            OPTIONAL{
                ?item wdt:P680 ?Molecular_Function }
            OPTIONAL{
                ?item wdt:P681 ?Cell_Component }
            OPTIONAL{
                ?item wdt:P682 ?Biological_Process }
            OPTIONAL{
                ?item wdt:P705 ?Ensembl_Protein_ID }
            OPTIONAL{
                ?item wdt:P637 ?RefSeq_Protein_ID }
            OPTIONAL{
                ?item wdt:P638 ?PDB_Structure_ID }

        }' AS sparql, m

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
                { Accept: "application/sparql-results+json"}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        WITH m, row, CASE WHEN EXISTS(m.A) THEN m.A ELSE [] END AS nA
            WHERE NONE (x IN nA WHERE x = row['Molecular_Function']['value'])
        SET m.A = nA + row['Molecular_Function']['value']

        WITH m, row, CASE WHEN EXISTS(m.A) THEN m.A ELSE [] END AS nA
            WHERE NONE (x IN nA WHERE x = row['Cell_Component']['value'])
        SET m.A = nA + row['Cell_Component']['value']

        WITH m, row, CASE WHEN EXISTS(m.A) THEN m.A ELSE [] END AS nA
            WHERE NONE (x IN nA WHERE x = row['Biological_Process']['value'])
        SET m.A = nA + row['Biological_Process']['value']

        WITH m, row, CASE WHEN EXISTS(m.A) THEN m.A ELSE [] END AS nA
            WHERE NONE (x IN nA WHERE x = row['Ensembl_Protein_ID']['value'])
        SET m.A = nA + row['Ensembl_Protein_ID']['value']

        WITH m, row, CASE WHEN EXISTS(m.A) THEN m.A ELSE [] END AS nA
            WHERE NONE (x IN nA WHERE x = row['PDB_Structure_ID']['value'])
        SET m.A = nA + row['PDB_Structure_ID']['value']
        """)

def remove_ExternalEquivalent(tx, number=None):
    """
    Removes all nodes of type: ExternalEquivalent from he DataBase; since this do not add
    new info, one might consider them not useful.
    """
    return tx.run("""
        MATCH (e:ExternalEquivalent)
        DETACH DELETE e
        """)
