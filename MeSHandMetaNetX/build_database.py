#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

# This is just a collection of functions used by the "main" script

def add_mesh_by_name(tx, query = None ):
    """
    A function that adds some MeSH nodes to any existing nodes, based on their Name property.
    NOTE: Only exact matches work here, which is not ideal.
    """
    return tx.run("""
        MATCH (n)
        WITH '
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
        PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
        PREFIX mesh2022: <http://id.nlm.nih.gov/mesh/2022/>
        PREFIX mesh2021: <http://id.nlm.nih.gov/mesh/2021/>
        PREFIX mesh2020: <http://id.nlm.nih.gov/mesh/2020/>

        SELECT DISTINCT  *
        FROM <http://id.nlm.nih.gov/mesh>
        WHERE
        { ?MeSH_Descriptor_ID
                     a              meshv:Descriptor ;
                     meshv:concept  ?MeSH_Concept_ID ;
                     rdfs:label     ?MeSH_Descriptor_Name .
          ?MeSH_Concept_ID
                      rdfs:label     ?MeSH_Concept_Name
          FILTER ( ( ?MeSH_Descriptor_Name = \"' +  n.Name + '\"@en ) || ( ?MeSH_Concept_Name = \"' +  n.Name + '\"@en ) )
        }' AS sparql, n

            CALL apoc.load.jsonParams(
                replace("https://id.nlm.nih.gov/mesh/sparql?query=" + sparql + "&format=JSON&limit=50&offset=0&inference=true", "\n", ""),
                { Accept: "application/sparql-results+json" }, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        MERGE (c:MeSH { MeSH_ID:split(row['MeSH_Descriptor_ID']['value'],'/')[-1] , Type:"Descriptor", Name: row['MeSH_Descriptor_Name']['value']})
        MERGE (n)-[r:RELATED_MESH]->(c)
        MERGE (cc:MeSH { MeSH_ID:split(row['MeSH_Concept_ID']['value'],'/')[-1] , Type:"Concept", Name: row['MeSH_Concept_Name']['value']})
        MERGE (n)-[r2:RELATED_MESH]->(cc)
        """)

def write_synonyms_in_metanetx(tx, query = None ):
    """
    A SPARQL function that finds synonyms for metabolites, proteins or drugs based on a given `query`, using MetaNetX.
    This query can either be Name, KEGG_ID, CHEBI, HMDB ID, InChI or InChIKey
    NOTE: The `RETURN databasename, databaseid` statement means we will` be able to use this in a read transaction too!
    """
    if query == "Name":
      query_text = """?mnx_url  rdfs:comment  \"' +  n.Name + '\" . """
    elif query == "KEGG_ID":
      query_text = """?mnx_url  mnx:chemXref  keggC:' +  n.KEGG_ID + '  . """
    elif query == "ChEBI":
      query_text = """?mnx_url  mnx:chemXref  chebi:' +  n.ChEBI_ID + ' . """
    elif query == "HMDB":
      query_text = """?mnx_url  mnx:chemXref  hmbd:' +  n.HMDB_ID + ' . """
    elif query == "InChI":
      query_text = """?mnx_url  mnx:inchi  \"' +  n.InChI + '\" . """
    elif query == "InChIKey":
      query_text = """?mnx_url  mnx:inchikey  \"InChIKey=' +  n.InChIKey + '\" . """

    else:
      query_text = ["" ""]

    return tx.run(f"""
        MATCH (n)
        WHERE (n:Metabolite OR n:Protein OR n:Drug)
        WITH '

        { add_prefixes() }

        SELECT DISTINCT *
        WHERE
            {{  ?mnx_url            rdf:type        mnx:CHEM .
                {query_text}
                ?mnx_url            rdfs:comment    ?Name;
                                    mnx:inchi       ?InChI ;
                                    mnx:inchikey    ?InChIKey ;
                                    mnx:smiles      ?SMILES ;
                                    mnx:formula     ?Formula ;
                                    mnx:mass        ?Mass ;
                                    mnx:chemXref    ?cross_refference .
                ?cross_refference
                                    rdfs:label      ?external_identifier
            }}' AS sparql, n

        CALL apoc.load.jsonParams(
                replace("https://rdf.metanetx.org/sparql/?query=" + sparql, "\n", ""),
                {{ Accept: "application/sparql-results+json" }}, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        FOREACH(ignoreme in case when n.Name IS NULL OR n.Name = "" then [1] else [] end |
                SET n.Name = row["Name"]["value"]
        )

        { get_identifiers() }
        """)

def read_synonyms_in_metanetx(tx, querytype = None, query = None ):
    """


    TODO: HACER :param :return

    """
    if querytype == "Name":
      query_text = f"""?mnx_url  rdfs:comment  \"{ query }\" . """
    elif querytype == "ChEBI":
      query_text = f"""?mnx_url  mnx:chemXref  chebi:{ query } . """
    elif querytype == "HMDB":
      query_text = f"""?mnx_url  mnx:chemXref  hmdb:{ query } . """
    elif querytype == "InChI":
      query_text = f"""?mnx_url  mnx:inchi  \"{ query }\" . """
    elif querytype == "InChIKey":
      query_text = f"""?mnx_url  mnx:inchikey  \"InChIKey={ query }\" . """

    graph_response = tx.run(f"""
        WITH '

        { add_prefixes() }

        SELECT DISTINCT *
        WHERE
            {{  ?mnx_url            rdf:type        mnx:CHEM .
                {query_text}
                ?mnx_url            rdfs:comment    ?Name;
                                    mnx:inchi       ?InChI ;
                                    mnx:inchikey    ?InChIKey ;
                                    mnx:chemXref    ?cross_refference .
                ?cross_refference
                                    rdfs:label      ?external_identifier
            }}' AS sparql

        CALL apoc.load.jsonParams(
                replace("https://rdf.metanetx.org/sparql/?query=" + sparql, "\n", ""),
                {{ Accept: "application/sparql-results+json" }}, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        WITH
            toLower(split(row["external_identifier"]["value"], ":")[0]) as databasename,
            split(row["external_identifier"]["value"], ":")[-1] as databaseid,
            row["InChI"]["value"] as InChI, row["InChIKey"]["value"] as InChIKey,
            row["Name"]["value"] as Name

        RETURN databasename, databaseid, InChI, InChIKey, Name

        """)

    return [record.data() for record in graph_response]


def find_protein_interactions_in_metanetx(tx):
    """
    A SPARQL function that finds the Metabolites a given Protein (based on its UniProt_ID) interacts with, using MetaNetX.
    NOTE: We are not using peptXref: since all proteins in MetaNetX come from UniProt, there is no use here
    """
    return tx.run(f"""
        MATCH (p)
        WHERE (p:Metabolite OR p:Protein OR p:Drug)
        WITH '
        { add_prefixes() }

        SELECT DISTINCT  *
        WHERE
            {{ ?pept    mnx:peptXref  up:' +  p.UniProt_ID + ' ;
                        mnx:geneName  ?geneName ;
                        rdfs:comment  ?Name.
               ?cata    mnx:pept      ?pept .
               ?gpr     mnx:cata      ?cata ;
                        mnx:reac      ?reac .
               ?reac    ?side         ?part .
               ?part    mnx:chem      ?chem ;
                        mnx:comp      ?component .
               ?chem    rdfs:comment  ?Name ;
                        mnx:inchi     ?InChI ;
                        mnx:inchikey  ?InChIKey ;
                        mnx:smiles    ?SMILES ;
                        mnx:formula   ?Formula ;
                        mnx:mass      ?Mass ;
                        mnx:chemXref  ?cross_refference .
               ?cross_refference
                         rdfs:label    ?external_identifier .
               ?component
                         rdfs:comment  ?CelularLocation
            }}' AS sparql, p

        CALL apoc.load.jsonParams(
                replace("https://rdf.metanetx.org/sparql/?query=" + sparql, "\n", ""),
                {{ Accept: "application/sparql-results+json" }}, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        MERGE (c:CelularLocation)
        SET c.Name = row["CelularLocation"]["value"]
        MERGE (p)-[r:LOCATED_INSIDE_CELL]->(c)

        FOREACH(ignoreme in case when p.Name IS NULL OR p.Name = "" then [1] else [] end |
                SET p.Name = row["Name"]["value"]
        )
        MERGE (g:Gene {{ Name:row["geneName"]["value"]}})
        MERGE (g)-[:ENCODES]->(p))

        MERGE (n:Metabolite {{ Name: row["Name"]["value"] }})
        MERGE (p)-[r:INTERACTS_WITH]->(n)

        { get_identifiers() }
        """)

def get_kegg_pathways_for_metabolites(tx):
    """
    A function that finds the Pathways a given Metabolite (based on its Kegg_ID) is a part of, using KEGG.
    This uses genome.jp's dbget web service, since I honestly could not find a way to use KEGG's SPARQL
    service (https://www.genome.jp/linkdb/linkdb_rdf.html) for that.
    Another possibility could be using https://rest.kegg.jp/get/R00703/
    """
    return tx.run("""
        MATCH (n)
        WHERE (n:Metabolite OR n:Protein OR n:Drug)

        LOAD CSV FROM "https://www.genome.jp/dbget-bin/get_linkdb?targettype=all&keywords=cpd:" +  n.KEGG_ID + "&targetformat=text" AS line FIELDTERMINATOR '\t'
        WITH
                split(line[1],":")[0] AS databasename,
                split(line[1],":")[1] AS databaseid,
                n
        FOREACH(ignoreme in case when databasename = "path" then [1] else [] end |
                MERGE (p:Pathway {KEGG_ID:databaseid})
                MERGE (n)-[r:PART_OF_PATHWAY]-(p)
        )
        """)

def add_prefixes():
    """Add some prefixes necessary for all MetaNetX queries to work. This are kept together since adding extra prefixes does not increase compu time"""
    return("""
        PREFIX mnx: <https://rdf.metanetx.org/schema/>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

        PREFIX up: <http://purl.uniprot.org/uniprot/>
        PREFIX chebi: <http://purl.obolibrary.org/obo/CHEBI_>
        PREFIX hmdb: <https://identifiers.org/hmdb:>
        PREFIX keggC: <https://identifiers.org/kegg.compound:>
        PREFIX keggR: <https://identifiers.org/kegg.reaction:>
           """)


def get_identifiers():
    """
    Part of a CYPHER query that processes the outcome from a SPARQL query that searches for information on MetaNetX
    It takes an original metabolite (n) and a row variable, which should have columns named external_identifier,
    cross_refference, InChIKey, InChI, SMILES, Formula and Mass with the adequated format; it is basically a code-reuser,
    not intended to be used separately.
    NOTE: All HMDB matches might create a Metabolite without CHEBI_ID or CAS_Number, which would violate our schema. This will be later on accounted for.
    NOTE: Some keys, such as VMH_ID, are not merged into their own node, but rather added to an existing one.
          This is because this do not prevously exist in our Schema, and might be changed in the future.
    NOTE: We dont care about overwriting InChI and InChIKey because they are necessarily unique; the same is true for Mass and Formula,
          as they are not all that important. However, for HMDB ID and others, we will take care not to overwrite, which could mess up the DB
    """
    return ("""
        WITH
            toLower(split(row["external_identifier"]["value"], ":")[0]) as databasename,
            split(row["external_identifier"]["value"], ":")[-1] as databaseid,
            row["cross_refference"]["value"] as url, row["InChI"]["value"] as InChI,
            row["InChIKey"]["value"] as InChIKey, row["SMILES"]["value"] as SMILES,
            row["Formula"]["value"] as Formula, row["Mass"]["value"] as Mass,
            n

        FOREACH(ignoreme in case when InChIKey IS NOT NULL AND InChIKey <> "" then [1] else [] end |
            SET n.InChIKey = InChIKey
        )
        FOREACH(ignoreme in case when InChI IS NOT NULL AND InChI <> "" then [1] else [] end |
            SET n.InChI = InChI
        )

        SET n.Formula = Formula, n.Average_Mass = Mass, n.SMILES = SMILES

        FOREACH(ignoreme in case when databasename = "chebi" AND n.ChEBI_ID <> databaseid then [1] else [] end |
            MERGE (m {ChEBI_ID: databaseid})
            MERGE (n)-[r:SYNONYM_OF]-(m)
            FOREACH(ignoreme in case when size(labels(m)) < 1 then [1] else [] end |
                SET m:Metabolite
            )
        )

        FOREACH(ignoreme in case when databasename = "kegg.compound" OR databasename = "keggc" AND n.KEGG_ID <> databaseid then [1] else [] end |
            MERGE (m {KEGG_ID: split(databaseid, "M_")[-1]})
            MERGE (n)-[r:SYNONYM_OF]-(m)
            FOREACH(ignoreme in case when size(labels(m)) < 1 then [1] else [] end |
                SET m:Metabolite
            )
        )

        FOREACH(ignoreme in case when databasename = "hmdb" AND n.HMDB_ID <> databaseid then [1] else [] end |
            MERGE (m {HMDB_ID: split(databaseid, "M_")[-1]})
            MERGE (n)-[r:SYNONYM_OF]-(m)
            FOREACH(ignoreme in case when size(labels(m)) < 1 then [1] else [] end |
                SET m:Metabolite
            )
        )

        FOREACH(ignoreme in case when databasename = "vmhmetabolite" then [1] else [] end |
            SET n.VMH_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "biggm" OR databasename = "bigg.metabolite" then [1] else [] end |
            SET n.BiGG_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "sabiork.compound" OR databasename = "sabiorkm" then [1] else [] end |
            SET n.SabioRK_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "seed.compound" then [1] else [] end |
            SET n.ModelSeed_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "seedm" then [1] else [] end |
            SET n.ModelSeed_ID = split(databaseid, "M_")[-1]
        )

        FOREACH(ignoreme in case when databasename = "metacyc.compound" OR databasename = "metacycm" then [1] else [] end |
            SET n.MetaCyc_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "slm" then [1] else [] end |
            SET n.SwissLipids_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "envipath" OR databasename = "envipathm" then [1] else [] end |
            SET n.EnviPath_ID = databaseid
        )

        WITH n, databasename, databaseid
        MATCH (m) WHERE (m:Metabolite OR m:Protein)

        FOREACH(ignoreme in case when databasename = "hmdb" then [1] else [] end |
            FOREACH(ignoreme in case when databaseid IN n.Secondary_HMDB_IDs then [1] else [] end |
                MERGE (n)-[r:SYNONYM_OF]-(m)
            )
        )
        """)
