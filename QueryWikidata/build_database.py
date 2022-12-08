#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that provides the necessary functions to transition selected parts of the Wikidata database to graph format,
either from scratch importing all the nodes (as showcased in :obj:`CanGraph.QueryWikidata.main`) or in a case-by-case basis,
to annotate existing metabolites (as showcased in :obj:`CanGraph.main`).

.. NOTE:: You may notice some functions here present the ``**kwargs`` arguments option.
    This is in order to make the functions compatible with the
    :obj:`CanGraph.miscelaneous.manage_transaction` function, which might send back a variable
    number of arguments (although technically it could work without the ``**kwargs`` option)
"""

def initial_cancer_discovery():
    """
    A Neo4J Cypher Statment that queries wikidata for Human Cancers. Since using the "afflicts:human"
    tag didnt have much use here, I used a simple workaround: Query wikidata for all humans, and, among them,
    find all of this for which their cause of death was a subclass of "Cancer" (Q12078). Unfortunaltely,
    some of them were diagnosed "Cancer" (Q12078), which is too general, so I removed it.

    Args:
        tx          (neo4j.Session): The session under which the driver is running

    Returns:
        neo4j.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.
    """
    return ("""
        WITH 'SELECT DISTINCT ?cause_of_death ?cause_of_death_name
        WHERE {
                ?human wdt:P31 wd:Q5;
                wdt:P509 ?cause_of_death.
                ?cause_of_death wdt:P279* wd:Q12078.
        SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                 ?cause_of_death rdfs:label ?cause_of_death_name. }
                }' AS sparql

        CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                { Accept: "application/sparql-results+json" }, null )
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

def find_subclass_of_disease():
    """
    A Neo4J Cypher Statment that queries wikidata for subclasses of "Disease" nodes already
    present on the Database. Since these are expected to only affect humans, this subclasses
    should also, only affect humans

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.

    .. NOTE:: We are forcing c.WikiData_ID to not be null or "". This is not necessary if we are just building the wikidata database,
        because there will always be a WikiData_ID, but it is useful in the rest of the cases
    """
    return """
        CALL {
        MATCH (c:Disease)
        WHERE c.WikiData_ID IS NOT NULL AND c.WikiData_ID <> ""
        WITH 'SELECT DISTINCT ?cancer ?cancer_name
        WHERE {
                ?cancer wdt:P279 wd:' + c.WikiData_ID + ' .
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                 ?cancer rdfs:label ?cancer_name. }
                }' AS sparql, c

        CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                { Accept: "application/sparql-results+json" }, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        CREATE (s:Disease)
        SET s.Name = row['cancer_name']['value'],
            s.URL = row['cancer']['value'],
            s.WikiData_ID = split(row['cancer']['value'],'/')[-1]

        CREATE (s)-[:SUBCLASS_OF]->(c)
        } IN TRANSACTIONS OF 10 rows
        """

def find_instance_of_disease():
    """
    A Neo4J Cypher Statment that queries wikidata for instances of "Disease" nodes already
    present on the Database. Since these are expected to only affect humans, this subclasses
    should also, only affect humans

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.

    .. NOTE:: We are forcing c.WikiData_ID to not be null or "". This is not necessary if we are just building the wikidata database,
        because there will always be a WikiData_ID, but it is useful in the rest of the cases
    """
    return """
        CALL {
        MATCH (c:Disease)
        WHERE c.WikiData_ID IS NOT NULL AND c.WikiData_ID <> ""
        WITH 'SELECT DISTINCT ?cancer ?cancer_name
        WHERE {
                ?cancer wdt:P31 wd:' + c.WikiData_ID + ' .
                SERVICE wikibase:label { bd:serviceParam wikibase:language "en".
                                 ?cancer rdfs:label ?cancer_name. }
                }' AS sparql, c

        CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                { Accept: "application/sparql-results+json" }, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        CREATE (i:Disease)
        SET i.Name = row['cancer_name']['value'],
            i.URL = row['cancer']['value'],
            i.WikiData_ID = split(row['cancer']['value'],'/')[-1]

        CREATE (i)-[:INSTANCE_OF]->(c)
        } IN TRANSACTIONS OF 10 rows
        """

def add_disease_info(number, **kwargs):
    """
    Adds info to "Disease" nodes for which its WikiData_ID ends in a given number. This way, only some of the nodes
    are targeted, and the Java Virtual Machine does not run out of memory

    Args:
        number (int): From 0 to 9, the number under which the WikiData_IDs to process should ends.
            This allows us tho divide the work, although its not very elegant.
        **kwargs: Any number of arbitrary keyword arguments

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.

    .. NOTE:: Here, there is no need to force c.WikiData_ID to not be null
        or "" because it will already be = ``number`` (and, thus, exist)
    """
    return f"""
        CALL {{
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
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
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

        }} IN TRANSACTIONS OF 10 rows
        """

def add_drugs(number, **kwargs):
    """
    Creates drug nodes related with each of the "Cancer" nodes already on the database

    Args:
        number (int): From 0 to 9, the number under which the WikiData_IDs to process should ends.
            This allows us tho divide the work, although its not very elegant.
        **kwargs: Any number of arbitrary keyword arguments

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.

    .. NOTE:: Here, there is no need to force c.WikiData_ID to not be null
        or "" because it will already be = ``number`` (and, thus, exist)
    """
    return f"""
        CALL {{

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
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        FOREACH(ignoreMe IN CASE WHEN row['Drugs'] IS NOT NULL THEN [1] ELSE [] END |
                MERGE (d:Drug{{URL:row['Drugs']['value'],
                               WikiData_ID:split(row['Drugs']['value'],'/')[-1],
                               Name:row['DrugName']['value'] }} )
                MERGE (c)-[:TARGETED_BY]->(d))

        FOREACH(ignoreMe IN CASE WHEN row['Exact_Matches'] IS NOT NULL THEN [1] ELSE [] END |
                MERGE (e:ExternalEquivalent{{URL:row['Exact_Matches']['value'],
                                            WikiData_ID:split(row['Exact_Matches']['value'],'/')[-1]}})
                MERGE (c)-[:EQUALS]-(e))

        }} IN TRANSACTIONS OF 10 rows
        """

def add_causes(number, **kwargs):
    """
    Creates drug nodes related with each of the "Cancer" nodes already on the database

    Args:
        number (int): From 0 to 9, the number under which the WikiData_IDs to process should ends.
            This allows us tho divide the work, although its not very elegant.
        **kwargs: Any number of arbitrary keyword arguments

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.

    .. NOTE:: Here, there is no need to force c.WikiData_ID to not be null
        or "" because it will already be = ``number`` (and, thus, exist)
    """
    return f"""
        CALL {{

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
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        FOREACH(ignoreMe IN CASE WHEN row['Causes'] IS NOT NULL THEN [1] ELSE [] END |
                MERGE (ca:Cause{{URL:row['Causes']['value'],
                                WikiData_ID:split(row['Causes']['value'],'/')[-1],
                                Name:row['CauseName']['value'] }} )
                MERGE (c)-[:CAUSED_BY]->(ca))

        }} IN TRANSACTIONS OF 10 rows
        """

def add_genes(number, **kwargs):
    """
    Creates gene nodes related with each of the "Cancer" nodes already on the database

    Args:
        number (int): From 0 to 9, the number under which the WikiData_IDs to process should ends.
            This allows us tho divide the work, although its not very elegant.
        **kwargs: Any number of arbitrary keyword arguments

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.

    .. NOTE:: Here, there is no need to force c.WikiData_ID to not be null
        or "" because it will already be = ``number`` (and, thus, exist)
    """
    return f"""
        CALL {{

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
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        FOREACH(ignoreMe IN CASE WHEN row['Genetic_Associations'] is not null and row['GeneName'] IS NOT NULL THEN [1] ELSE [] END |
                MERGE (g:Gene{{URL:row['Genetic_Associations']['value'],
                              WikiData_ID:split(row['Genetic_Associations']['value'],'/')[-1],
                              Name:row['GeneName']['value'] }})
                MERGE (g)-[:ASSOCIATED_DISEASE_GENE]->(c))

        }} IN TRANSACTIONS OF 10 rows
        """

def add_drug_external_ids(query = "Wikidata_ID", **kwargs):
    """
    Adds some external IDs to any "Drug" nodes already present on the database.
    Since the PDB information had too much values which caused triple duplicates that overcharged the system,
    they were intentionally left out.

    Args:
        query (str): One of ["DrugBank_ID","WikiData_ID"], a way to identify the nodes for which external IDs will be added.
        **kwargs: Any number of arbitrary keyword arguments

    Returns:
        neo4j.Result:
            A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: We are forcing c.WikiData_ID to not be null or "". This is not necessary if we are just building the wikidata database,
        because there will always be a WikiData_ID, but it is useful in the rest of the cases
    """
    if query == "DrugBank_ID":
      query_text = [""" WHERE d.DrugBank_ID IS NOT NULL AND d.DrugBank_ID <> "" """,
                    """?item (wdt:P715) \"' + apoc.text.replace(d.DrugBank_ID, "DB", "") + '\" """]
    else:
      query_text = ["""WHERE d.WikiData_ID IS NOT NULL AND d.WikiData_ID <> "" """,
                    """filter (?item = wd:' + d.WikiData_ID + ')"""]


    return f"""
        CALL {{

        MATCH (d:Drug)
        {query_text[0]}
        WITH 'SELECT DISTINCT *
        WHERE {{
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGUAGE]". }}

        {query_text[1]}

        OPTIONAL{{
            ?item wdt:P715 ?DrugBank }}

        OPTIONAL{{
            ?item wdt:P683 ?ChEBI_ID }}
        OPTIONAL{{
            ?item wdt:P592 ?CHEMBL_ID }}
        OPTIONAL{{
            ?item wdt:P661 ?ChemSpider_ID }}
        OPTIONAL{{
            ?item wdt:P662 ?PubChem_ID }}
        OPTIONAL{{
            ?item wdt:P665 ?KEGG_ID }}
        OPTIONAL{{
            ?item wdt:P231 ?CAS_Number }}
        OPTIONAL{{
            ?item wdt:P267 ?ATC_Code }}

        OPTIONAL{{
            ?item wdt:P2175 ?Medical_condition_treated
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                    ?Medical_condition_treated rdfs:label ?Medical_condition_treated_name. }} }}

            }}' AS sparql, d

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                {{ Accept: "application/sparql-results+json" }}, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET d.WikiData_ID = split(row['item']['value'],'/')[-1],
            d.CHEMBL_ID = row['CHEMBL_ID']['value'],
            d.ChEBI_ID = row['ChEBI_ID']['value'],
            d.PubChem_ID = row['PubChem_ID']['value'],
            d.KEGG_ID = row['KEGG_ID']['value'],
            d.DrugBank_ID = row['DrugBank']['value'],
            d.CAS_Number = row['CAS_Number']['value'],
            d.ChemSpider_ID = row['ChemSpider_ID']['value']

        FOREACH(ignoreMe IN CASE WHEN row['ATC_Code'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (pri:ATC {{ Code:row['ATC_Code']['value'] }})
            MERGE (d)-[r:RELATED_ATC]->(pri)
        )

        FOREACH(ignoreMe IN CASE WHEN row['MeSH_Descriptor_Name'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (c:MeSH {{ MeSH_ID:row['MeSH_Descriptor_ID']['value'] ,
                             Type:"Descriptor", Name: row['MeSH_Descriptor_Name']['value']}})
            MERGE (d)-[r2:RELATED_MESH]->(c)
        )

        FOREACH(ignoreMe IN CASE WHEN row['MeSH_Concept_Name'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (cc:MeSH {{ MeSH_ID:row['MeSH_Concept_ID']['value'] ,
                              Type:"Concept", Name: row['MeSH_Concept_Name']['value']}})
            MERGE (d)-[r3:RELATED_MESH]->(cc)
        )

        }} IN TRANSACTIONS OF 10 rows
        """

def add_more_drug_info(query = "WikiData_ID", **kwargs):
    """
    Creates some nodes that are related with each of the "Drug" nodes already existing
    on the database: routes of administration, targeted metabolites and approved drugs
    that tehy are been used in

    Args:
        query (str): One of ["DrugBank_ID","WikiData_ID"], a way to identify the nodes for
            which external IDs will be added; default is "WikiData_ID"
        **kwargs: Any number of arbitrary keyword arguments

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.


    .. TODO:: ADD ROLE to metabolite interactions

    .. NOTE:: This transaction has been separated in order to keep response times low

    .. NOTE:: We are forcing c.WikiData_ID to not be null or "". This is not necessary if we are just building the wikidata database,
        because there will always be a WikiData_ID, but it is useful in the rest of the cases
    """
    if query == "DrugBank_ID":
      query_text = [""" WHERE d.DrugBank_ID IS NOT NULL AND d.DrugBank_ID <> "" """,
                    """?item (wdt:P715) \"' + apoc.text.replace(d.DrugBank_ID, "DB", "") + '\" """]
    else:
      query_text = ["""WHERE d.WikiData_ID IS NOT NULL AND d.WikiData_ID <> "" """,
                    """filter (?item = wd:' + d.WikiData_ID + ')"""]


    return f"""
        CALL {{

        MATCH (d:Drug)
        {query_text[0]}
        WITH 'SELECT DISTINCT *
        WHERE {{
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGUAGE]". }}

        {query_text[1]}

        OPTIONAL{{
            ?item wdt:P234 ?InChI }}
        OPTIONAL{{
            ?item wdt:P235 ?InChIKey }}
        OPTIONAL{{
            ?item wdt:P2067 ?Mass }}
        OPTIONAL{{
            ?item wdt:P274 ?Chemical_Formula }}
        OPTIONAL{{
            ?item wdt:P233 ?Canonical_Smiles }}
        OPTIONAL{{
            ?item wdt:P7830 ?LiverTox }}
            }}' AS sparql, d

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                {{ Accept: "application/sparql-results+json" }}, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET d.WikiData_ID = split(row['item']['value'],'/')[-1],
            d.InChI = row['InChI']['value'],
            d.InChIKey = row['InChIKey']['value'],
            d.Average_Mass = row['Mass']['value'],
            d.Formula = row['Chemical_Formula']['value'],
            d.SMILES = row['Canonical_Smiles']['value'],
            d.Pregnancy_Category = row['Pregnancy_Category_Name']['value'],
            d.LiverTox = row['LiverTox']['value']

        }} IN TRANSACTIONS OF 10 rows
        """

def add_yet_more_drug_info(query = "WikiData_ID", **kwargs):
    """
    Creates some nodes that are related with each of the "Drug" nodes already existing
    on the database: routes of administration, targeted metabolites and approved drugs
    that tehy are been used in

    Args:
        query (str): One of ["DrugBank_ID","WikiData_ID"], a way to identify the nodes for
            which external IDs will be added; default is "WikiData_ID"
        **kwargs: Any number of arbitrary keyword arguments

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.
    """
    if query == "DrugBank_ID":
      query_text = [""" WHERE d.DrugBank_ID IS NOT NULL AND d.DrugBank_ID <> "" """,
                    """?item (wdt:P715) ' + apoc.text.replace(d.DrugBank_ID, "DB", "") + '"""]
    else:
      query_text = ["""WHERE d.WikiData_ID IS NOT NULL AND d.WikiData_ID <> "" """,
                    """filter (?item = wd:' + d.WikiData_ID + ')"""]

    return f"""
        CALL {{

        MATCH (d:Drug)
        {query_text[0]}
        WITH 'SELECT DISTINCT *
        WHERE {{
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGUAGE]". }}

        {query_text[1]}

        OPTIONAL{{
            ?item wdt:P636 ?Route_of_Admin
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                    ?Route_of_Admin rdfs:label ?Route_of_Admin_name. }} }}
        OPTIONAL{{
            ?item wdt:P3780 ?Active_ingredient_in
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                    ?Active_ingredient_in rdfs:label ?Active_ingredient_in_name. }} }}
        OPTIONAL{{
            ?item wdt:P129 ?Interacts_with
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                    ?Interacts_with rdfs:label ?Interacts_with_name. }} }}
            }}' AS sparql, d

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET d.WikiData_ID = split(row['item']['value'],'/')[-1]

        FOREACH(ignoreMe IN CASE WHEN row['Interacts_with'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (m:Metabolite{{URL:row['Interacts_with']['value'],
                        WikiData_ID:split(row['Interacts_with']['value'],'/')[-1],
                        Name:row['Interacts_with_name']['value'] }})
            MERGE (d)-[:INTERACTS_WITH]-(m))

        FOREACH(ignoreMe IN CASE WHEN row['Route_of_Admin'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (r:AdministrationRoute{{URL:row['Route_of_Admin']['value'],
                        WikiData_ID:split(row['Route_of_Admin']['value'],'/')[-1],
                        Name:row['Route_of_Admin_name']['value'] }})
            MERGE (d)-[:ADMINISTERED_VIA]->(r))

        FOREACH(ignoreMe IN CASE WHEN row['Active_ingredient_in'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (me:Product{{URL:row['Active_ingredient_in']['value'],
                        WikiData_ID:split(row['Active_ingredient_in']['value'],'/')[-1],
                        Name:row['Active_ingredient_in_name']['value'] }})
            MERGE (d)-[:PART_OF_PRODUCT]->(me))

        }} IN TRANSACTIONS OF 10 rows
        """

def add_gene_info():
    """
    A Cypher Query that adds some external IDs and properties to "Gene" nodes already existing on
    the database. This query forces the genes to have a "found_in_taxon:homo_sapiens" label. This means
    that any non-human genes will not be annotated (.. TODO:: delete those)

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.

    .. NOTE:: Genomic Start and ends keep just the 2nd position, as reported in wikidata

    .. NOTE:: We are forcing c.WikiData_ID to not be null or "". This is not necessary if we are just building the wikidata database,
        because there will always be a WikiData_ID, but it is useful in the rest of the cases

    .. TODO:: Might include P684 "Orthologues" for more info (it crashed java)
    """
    return """
        CALL {

        MATCH (g:Gene)
        WHERE g.WikiData_ID IS NOT NULL AND g.WikiData_ID <> ""
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
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                { Accept: "application/sparql-results+json" }, null )
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

        FOREACH(ignoreMe IN CASE WHEN row['Encodes'] IS NOT NULL THEN [1] ELSE [] END |
                    MERGE (m:Metabolite{URL:row['Encodes']['value'],
                                WikiData_ID:split(row['Encodes']['value'],'/')[-1],
                                Name:row['Encodes_name']['value'] })
                    MERGE (g)-[:ENCODES]->(m))

        FOREACH(ignoreMe IN CASE WHEN row['Expressed_in'] IS NOT NULL THEN [1] ELSE [] END |
                    MERGE (t:Tissue{URL:row['Expressed_in']['value'],
                                WikiData_ID:split(row['Expressed_in']['value'],'/')[-1],
                                Name:row['Expressed_in_name']['value'] })
                    MERGE (g)-[:EXPRESSED_IN]->(t))

        FOREACH(ignoreMe IN CASE WHEN row['Exact_Matches'] IS NOT NULL THEN [1] ELSE [] END |
                MERGE (e:ExternalEquivalent{URL:row['Exact_Matches']['value'],
                                            WikiData_ID:split(row['Exact_Matches']['value'],'/')[-1]})
                MERGE (g)-[:EQUALS]-(e))

        } IN TRANSACTIONS OF 10 rows
        """

def add_metabolite_info(query = "ChEBI_ID", **kwargs):
    """
    A Cypher Query that adds some external IDs and properties to "Metabolite" nodes
    already existing on the database. Two kind of metabolites exist: those that are
    encoded by a given gene, and those that interact with a given drug. Both are
    adressed here, since they are similar, and, most likely, instances of proteins.

    * This function forces all metabolites to have a "found_in_taxon:human" target
    * The metabolites are not forced to be proteins, but if they are, this is kept in the "instance_of" record

    Args:
        query (str): One of ["DrugBank_ID","WikiData_ID"], a way to identify the nodes for
            which external IDs will be added; default is "WikiData_ID"
        **kwargs: Any number of arbitrary keyword arguments

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.

    .. TODO:: Might include P527 "has part or parts" for more info (it crashed java)

    .. NOTE:: We are forcing c.WikiData_ID to not be null or "". This is not necessary if we are just building the wikidata database,
        because there will always be a WikiData_ID, but it is useful in the rest of the cases
    """
    if query == "ChEBI_ID":
      query_text = [""" WHERE m.ChEBI_ID IS NOT NULL AND m.ChEBI_ID <> "" """,
                    """?item (wdt:P683) \"' +  m.ChEBI_ID + '\" ."""]
    else:
      query_text = ["""WHERE m.WikiData_ID IS NOT NULL AND m.WikiData_ID <> "" """,
                    """filter (?item = wd:' + m.WikiData_ID + ')"""]

    return f"""
        CALL {{

        MATCH (m:Metabolite)
        {query_text[0]}
        WITH 'SELECT DISTINCT *
        WHERE {{
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGUAGE]". }}

        {query_text[1]}

        OPTIONAL{{
            ?item wdt:P31 ?Instance_of
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                                     ?Instance_of rdfs:label ?Instance_of_name. }} }}

        OPTIONAL{{
            ?item (wdt:P683) ?ChEBI_ID }}
        OPTIONAL{{
            ?item wdt:P231 ?CAS_Number }}

        OPTIONAL{{
            ?item wdt:P234 ?InChI }}
        OPTIONAL{{
            ?item wdt:P235 ?InChIKey }}
        OPTIONAL{{
            ?item wdt:P352 ?UniProt_ID }}
        OPTIONAL{{
            ?item wdt:P7260 ?Transporter_DB_ID }}
        OPTIONAL{{
            ?item wdt:P274 ?Chemical_Formula }}
        OPTIONAL{{
            ?item wdt:P592 ?CHEMBL_ID }}
        OPTIONAL{{
            ?item wdt:P662 ?PubChem_ID }}
        OPTIONAL{{
            ?item wdt:P233 ?Canonical_Smiles }}
        OPTIONAL{{
            ?item wdt:P8117 ?FooDB_Compound_ID }}

        OPTIONAL{{
            ?item p:P486 ?statement1 .
            ?statement1 ps:P486 ?MeSH_Descriptor_ID .
            ?statement1 pq:P1810 ?MeSH_Descriptor_Name .
            }}
        OPTIONAL{{
            ?item p:P6694 ?statement2 .
            ?statement2 ps:P6694 ?MeSH_Concept_ID .
            ?statement2 pq:P1810 ?MeSH_Concept_Name . }}

        OPTIONAL{{
            ?item wdt:P129 ?Interacts_with
            SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en".
                    ?Interacts_with rdfs:label ?Interacts_with_name. }} }}
        OPTIONAL{{
            ?item wdt:P2888 ?Exact_Matches }}

            }}' AS sparql, m

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                {{ Accept: "application/sparql-results+json"}}, null)
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET m.WikiData_ID = split(row['item']['value'],'/')[-1],
            m.Function = row['Instance_of_name']['value'],
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

        FOREACH(ignoreMe IN CASE WHEN row['MeSH_Descriptor_Name'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (c:MeSH {{ MeSH_ID:row['MeSH_Descriptor_ID']['value'] ,
                             Type:"Descriptor", Name: row['MeSH_Descriptor_Name']['value']}})
            MERGE (m)-[r:RELATED_MESH]->(c)
        )

        FOREACH(ignoreMe IN CASE WHEN row['MeSH_Concept_Name'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (cc:MeSH {{ MeSH_ID:row['MeSH_Concept_ID']['value'] ,
                              Type:"Concept", Name: row['MeSH_Concept_Name']['value']}})
            MERGE (m)-[r2:RELATED_MESH]->(cc)
        )

        FOREACH(ignoreMe IN CASE WHEN row['Interacts_with'] IS NOT NULL THEN [1] ELSE [] END |
                    MERGE (me:Metabolite{{URL:row['Interacts_with']['value'],
                                WikiData_ID:split(row['Interacts_with']['value'],'/')[-1],
                                Name:row['Interacts_with_name']['value'] }})
                    MERGE (m)-[:INTERACTS_WITH]-(me))

        FOREACH(ignoreMe IN CASE WHEN row['Exact_Matches'] IS NOT NULL THEN [1] ELSE [] END |
                MERGE (e:ExternalEquivalent{{URL:row['Exact_Matches']['value'],
                                            WikiData_ID:split(row['Exact_Matches']['value'],'/')[-1]}})
                MERGE (m)-[:EQUALS]-(e))

        }} IN TRANSACTIONS OF 10 rows
        """

def add_toomuch_metabolite_info():
    """
    A function that adds loads of info to existing "Metabolite" nodes. This was left out, first because
    it might be too much information, (specially when it is already availaible by clicking the "url" field),
    and because, due to it been so much, it crashes the JVM.

    Returns:
        str: A CYPHER query that modifies the DB according to the CYPHER statement contained in the function.
    """
    return """
        CALL {
        MATCH (m:Metabolite)
        WHERE m.WikiData_ID IS NOT NULL AND m.WikiData_ID <> ""
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
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                { Accept: "application/sparql-results+json" }, null )
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

        } IN TRANSACTIONS OF 10 rows
        """

def add_wikidata_and_mesh_by_name():
    """
    A function that adds some MeSH nodes and WikiData_IDs
    to existing nodes, based on their Wikipedia Article Title.

    Args:
        tx (neo4j.Session): The session under which the driver is running

    Returns:
        neo4j.Result:
            A Neo4J connexion to the database that modifies it according
            to the CYPHER statement contained in the function.
    """
    return """
        CALL {
        MATCH (d)
        WHERE (d:AdministrationRoute OR d:Cause OR d:Disease OR d:Drug OR d:ExternalEquivalent
                                     OR d:Gene OR d:Metabolite OR d:Product OR d:Tissue)
        AND (d.Name IS NOT null)
        WITH 'SELECT DISTINCT *
        WHERE{
            ?sitelink schema:about ?WikiData_ID;
                    schema:isPartOf <https://en.wikipedia.org/>;
                    schema:name \"' +  replace(d.Name, "  ", " ") + '\"@en.
        OPTIONAL{
            ?WikiData_ID p:P486 ?statement1 .
            ?statement1 ps:P486 ?MeSH_Descriptor_ID .
            ?statement1 pq:P1810 ?MeSH_Descriptor_Name .
            }
        OPTIONAL{
            ?WikiData_ID p:P6694 ?statement2 .
            ?statement2 ps:P6694 ?MeSH_Concept_ID .
            ?statement2 pq:P1810 ?MeSH_Concept_Name . }
        }' AS sparql, d

            CALL apoc.load.jsonParams(
                replace("https://query.wikidata.org/sparql?query=" + apoc.text.urlencode(sparql), "\n", ""),
                { Accept: "application/sparql-results+json" }, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET d.WikiData_ID = split(row['WikiData_ID']['value'],'/')[-1]

        FOREACH(ignoreMe IN CASE WHEN row['MeSH_Descriptor_Name'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (c:MeSH { MeSH_ID:row['MeSH_Descriptor_ID']['value'] ,
                            Type:"Descriptor", Name: row['MeSH_Descriptor_Name']['value']})
            MERGE (d)-[r:RELATED_MESH]->(c)
        )

        FOREACH(ignoreMe IN CASE WHEN row['MeSH_Concept_Name'] IS NOT NULL THEN [1] ELSE [] END |
            MERGE (cc:MeSH { MeSH_ID:row['MeSH_Concept_ID']['value'] ,
                            Type:"Concept", Name: row['MeSH_Concept_Name']['value']})
            MERGE (d)-[r2:RELATED_MESH]->(cc)
        )

        } IN TRANSACTIONS OF 10 rows
        """
