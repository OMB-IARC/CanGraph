#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that provides the necessary functions to transition the MetaNetX database (and related MeSH terms and KEGG IDs) to graph format,
either from scratch importing all the nodes (as showcased in :obj:`CanGraph.MeSHandMetaNetX.main`) or in a case-by-case basis,
to annotate existing metabolites (as showcased in :obj:`CanGraph.main`).
"""

# ********* SPARQL queries to annotate existing nodes using MeSH ********* #

def add_mesh_by_name():
    """
    A function that adds some MeSH nodes to any existing nodes, based on their Name property.

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`

    .. NOTE:: Only exact matches work here, which is not ideal.
    """
    return """
        CALL {

        MATCH (n)
        WHERE n.Name IS NOT null AND n.Name <> ""
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

        } IN TRANSACTIONS OF 100 rows
        """

# ********* SPARQL queries to annotate existing nodes using MetaNetX ********* #

def add_prefixes():
    """
    Add some prefixes necessary for all MetaNetX queries to work.
    This are kept together since adding extra prefixes does not increase computation time

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`
    """
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


def get_identifiers(from_sparql=False):
    """
    Part of a CYPHER query that processes the outcome from a SPARQL query that searches for information on MetaNetX
    It takes an original metabolite (n) and a row variable, which should have columns named external_identifier,
    cross_refference, InChIKey, InChI, SMILES, Formula and Mass with the adequated format; it is basically a code-reuser,
    not intended to be used separately.

    Args:
        from_sparql (bool): A True/False param defining whether the identifiers are being parsed from a SPARQL query;
            default is False (i.e. imported from file)

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`

    .. NOTE:: All HMDB matches might create a Metabolite without CHEBI_ID or CAS_Number, which would violate our schema. This will be later on accounted for.

    .. NOTE:: Some keys, such as VMH_ID, are not merged into their own node, but rather added to an existing one.
          This is because this do not prevously exist in our Schema, and might be changed in the future.

    .. NOTE:: We dont care about overwriting InChI and InChIKey because they are necessarily unique; the same is true for Mass and Formula,
          as they are not all that important. However, for HMDB ID and others, we will take care not to overwrite, which could mess up the DB
    """
    if from_sparql == True:
        sparql_parser = ["""

        WITH
            toLower(split(row["external_identifier"]["value"], ":")[0]) as databasename,
            split(row["external_identifier"]["value"], ":")[-1] as databaseid,
            row["cross_refference"]["value"] as url, row["InChI"]["value"] as InChI,
            row["InChIKey"]["value"] as InChIKey, row["SMILES"]["value"] as SMILES,
            row["Formula"]["value"] as Formula, row["Mass"]["value"] as Mass,
            split(row['mnx_url']['value'],'/')[-1] as MetaNetX_ID,
            split(row['Isomers']['value'],'/')[-1] as Isomer,
            n

        SET n.MetaNetX_ID = MetaNetX_ID

        FOREACH(ignoreme in case when InChIKey IS NOT NULL AND InChIKey <> "" then [1] else [] end |
            SET n.InChIKey = InChIKey
        )
        FOREACH(ignoreme in case when InChI IS NOT NULL AND InChI <> "" then [1] else [] end |
            SET n.InChI = InChI
        )

        SET n.Formula = Formula, n.Average_Mass = Mass, n.SMILES = SMILES
            """,
            """
        WITH n, databasename, databaseid, MetaNetX_ID, Isomer
        MATCH (m) WHERE (m:Metabolite OR m:Protein)

        FOREACH(ignoreme in case when databasename = "hmdb" then [1] else [] end |
            FOREACH(ignoreme in case when databaseid IN n.Secondary_HMDB_IDs then [1] else [] end |
                MERGE (n)-[r:SYNONYM_OF]-(m)
                SET m.MetaNetX_ID = MetaNetX_ID
            )
        )

        MERGE (z:Metabolite { MetaNetX_ID:Isomer })
        MERGE (z)-[r:ISOMER_OF]-(n)
            """]
    else: sparql_parser = ["",""]

    return f"""

        {sparql_parser[0]}

        FOREACH(ignoreme in case when databasename = "chebi" AND n.ChEBI_ID <> databaseid then [1] else [] end |
            MERGE (m {{ChEBI_ID: databaseid}})
            MERGE (n)-[r:SYNONYM_OF]-(m)
            FOREACH(ignoreme in case when size(labels(m)) < 1 then [1] else [] end |
                SET m:Metabolite
            )
        )

        FOREACH(ignoreme in case when databasename = "kegg.compound" OR databasename = "keggc" AND n.KEGG_ID <> databaseid then [1] else [] end |
            MERGE (m {{KEGG_ID: split(databaseid, "M_")[-1]}})
            MERGE (n)-[r:SYNONYM_OF]-(m)
            FOREACH(ignoreme in case when size(labels(m)) < 1 then [1] else [] end |
                SET m:Metabolite
            )
        )

        FOREACH(ignoreme in case when databasename = "hmdb" AND n.HMDB_ID <> databaseid then [1] else [] end |
            MERGE (m {{HMDB_ID: split(databaseid, "M_")[-1]}})
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

        {sparql_parser[1]}
        """

def write_synonyms_in_metanetx(query = None ):
    """
    A SPARQL function that finds synonyms for metabolites, proteins or drugs in an existing Neo4J database, using MetaNetX.
    At the same time, it is able to annotate them a bit, adding Name, InChI, InChIKey, SMILES, Formula, Mass, some External IDs,
    and finding whether the metabolite in question has any known isomers, anootating if so.

    Args:
        query (str): The type of query that is being searched for. One of ["Name","KEGG_ID","ChEBI_ID","HMDB_ID","InChI","InChIKey"];
            default is None (i.e. the function will not process anything).

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`

    .. NOTE:: This is intended to be run as a write_transaction, modifying the existing database.
    """
    if query == "Name":
      query_text = """?mnx_url  rdfs:comment  \"' +  n.Name + '\" . """
    elif query == "KEGG_ID":
      query_text = """?mnx_url  mnx:chemXref  keggC:' +  n.KEGG_ID + '  . """
    elif query == "ChEBI_ID":
      query_text = """?mnx_url  mnx:chemXref  chebi:' +  n.ChEBI_ID + ' . """
    elif query == "HMDB_ID":
      query_text = """?mnx_url  mnx:chemXref  hmbd:' +  n.HMDB_ID + ' . """
    elif query == "InChI":
      query_text = """?mnx_url  mnx:inchi  \"' +  n.InChI + '\" . """
    elif query == "InChIKey":
      query_text = """?mnx_url  mnx:inchikey  \"InChIKey=' +  n.InChIKey + '\" . """

    else:
      query_text = ["" ""]

    return f"""
        CALL {{

        MATCH (n)
        WHERE (n:Metabolite OR n:Protein OR n:Drug) AND n.{query} IS NOT null AND n.{query} <> ""
        WITH '

        { add_prefixes() }

        SELECT DISTINCT *
        WHERE
            {{      ?mnx_url    rdf:type        mnx:CHEM .

                    {query_text}

        OPTIONAL{{  ?mnx_url    rdfs:comment    ?Name       }}
        OPTIONAL{{  ?mnx_url    mnx:inchi       ?InChI      }}
        OPTIONAL{{  ?mnx_url    mnx:inchikey    ?InChIKey   }}
        OPTIONAL{{  ?mnx_url    mnx:smiles      ?SMILES     }}
        OPTIONAL{{  ?mnx_url    mnx:formula     ?Formula    }}
        OPTIONAL{{  ?mnx_url    mnx:mass        ?Mass       }}
        OPTIONAL{{  ?mnx_url    mnx:hasIsomericChild    ?Isomers }}
                    ?mnx_url    mnx:chemXref    ?cross_refference .

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

        { get_identifiers(from_sparql = True) }

        }} IN TRANSACTIONS OF 100 rows
        """

def read_synonyms_in_metanetx(tx, querytype = None, query = None ):
    """
    A SPARQL function that finds synonyms for metabolites, proteins or drugs based on a given `query`, using MetaNetX.
    At the same time, it is able to annotate them a bit, adding Name, InChI, InChIKey, SMILES, Formula, Mass, some External IDs,
    and finding whether the metabolite in question has any known isomers, anootating if so.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        querytype   (str): The type of query that is being searched for. One of ["Name","KEGG_ID","ChEBI_ID","HMDB_ID","InChI","InChIKey"];
            default is :obj:`None` (i.e. the function will not process anything).
        query       (str): The query we are searching for; must be of type ```querytype```


    .. NOTE:: This is intended to be run as a read_transaction, only returning synonyms present in the DB. No modifications will be applied.
    """
    if querytype == "Name":
      query_text = f"""?mnx_url  rdfs:comment  \"{ query }\" . """
    elif querytype == "ChEBI_ID":
      query_text = f"""?mnx_url  mnx:chemXref  chebi:{ query } . """
    elif querytype == "HMDB_ID":
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
            {{      ?mnx_url            rdf:type        mnx:CHEM .

                    {query_text}

        OPTIONAL{{  ?mnx_url    rdfs:comment    ?Name       }}
        OPTIONAL{{  ?mnx_url    mnx:inchi       ?InChI      }}
        OPTIONAL{{  ?mnx_url    mnx:inchikey    ?InChIKey   }}
        OPTIONAL{{  ?mnx_url    mnx:smiles      ?SMILES     }}
        OPTIONAL{{  ?mnx_url    mnx:formula     ?Formula    }}
        OPTIONAL{{  ?mnx_url    mnx:mass        ?Mass       }}
        OPTIONAL{{  ?mnx_url    mnx:hasIsomericChild    ?Isomers }}
                    ?mnx_url    mnx:chemXref    ?cross_refference .

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


def find_protein_interactions_in_metanetx():
    """
    A SPARQL function that finds the Metabolites a given Protein (based on its UniProt_ID) interacts with, using MetaNetX.

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`

    .. NOTE:: We are not using peptXref: since all proteins in MetaNetX come from UniProt, there is no use here
    """
    return f"""
        CALL {{

        MATCH (p)
        WHERE (p:Metabolite OR p:Protein OR p:Drug) AND p.UniProt_ID IS NOT null AND p.UniProt_ID <> ""
        WITH '

        { add_prefixes() }

        SELECT DISTINCT  *
        WHERE
            {{      ?pept       mnx:peptXref  up:' +  p.UniProt_ID + ' .

                    ?cata       mnx:pept      ?pept .
                    ?gpr        mnx:cata      ?cata ;
                                mnx:reac      ?reac .
                    ?reac       ?side         ?part .
                    ?part       mnx:chem      ?mnx_url .

        OPTIONAL{{  ?mnx_url    rdfs:comment    ?Name       }}
        OPTIONAL{{  ?mnx_url    mnx:inchi       ?InChI      }}
        OPTIONAL{{  ?mnx_url    mnx:inchikey    ?InChIKey   }}
        OPTIONAL{{  ?mnx_url    mnx:smiles      ?SMILES     }}
        OPTIONAL{{  ?mnx_url    mnx:formula     ?Formula    }}
        OPTIONAL{{  ?mnx_url    mnx:mass        ?Mass       }}
        OPTIONAL{{  ?mnx_url    mnx:hasIsomericChild    ?Isomers }}
                    ?mnx_url    mnx:chemXref  ?cross_refference .
               ?cross_refference
                            rdfs:label    ?external_identifier .

            }}' AS sparql, p

        CALL apoc.load.jsonParams(
            replace("https://rdf.metanetx.org/sparql/?default-graph-uri=https://rdf.metanetx.org/&query=" + sparql, "\n", ""),
            {{ Accept: "application/sparql-results+json" }}, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET p.MetaNetX_ID = split(row['pept']['value'],'/')[-1]

        MERGE (n:Metabolite {{ MetaNetX_ID: split(row['mnx_url']['value'],'/')[-1] }})
        SET n.Name = row["Name"]["value"]
        MERGE (p)-[r:INTERACTS_WITH]->(n)

        { get_identifiers(from_sparql = True) }

        }} IN TRANSACTIONS OF 100 rows
        """

def find_protein_data_in_metanetx():
    """
    A SPARQL function that annotates Protein nodes in an exiting Neo4J database by using the information provided by MetaNetX

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`

    .. NOTE:: This function is partly a duplicate of `self.find_protein_interactions_in_metanetx()`, which was split to prevent timeouts
    """

    return f"""
        CALL {{

        MATCH (p)
        WHERE (p:Metabolite OR p:Protein OR p:Drug) AND p.UniProt_ID IS NOT null AND p.UniProt_ID <> ""

        WITH '

        { add_prefixes() }

        SELECT DISTINCT  *
        WHERE
            {{      ?pept       mnx:peptXref  up:' +  p.UniProt_ID + ' .
        OPTIONAL{{  ?pept       mnx:geneName       ?geneName      }}
        OPTIONAL{{  ?pept       rdfs:comment       ?peptName      }}

                    ?cata       mnx:pept      ?pept .
                    ?gpr        mnx:cata      ?cata ;
                                mnx:reac      ?reac .
                    ?reac       ?side         ?part .

        OPTIONAL{{
                    ?part       mnx:comp      ?component .
                    ?component  rdfs:comment  ?CelularLocation }}
            }}' AS sparql, p

        CALL apoc.load.jsonParams(
            replace("https://rdf.metanetx.org/sparql/?query=" + sparql, "\n", ""),
            {{ Accept: "application/sparql-results+json" }}, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        SET p.MetaNetX_ID = split(row['pept']['value'],'/')[-1]

        FOREACH(ignoreme in case when row["CelularLocation"]["value"] IS NOT NULL then [1] else [] end |
            MERGE (c:CelularLocation {{Name:row["CelularLocation"]["value"]}})
            MERGE (p)-[r:LOCATED_INSIDE_CELL]->(c)
        )

        FOREACH(ignoreme in case when p.Name IS NULL OR p.Name = "" then [1] else [] end |
            SET p.Name = row["peptName"]["value"]
        )

        FOREACH(ignoreme in case when row["geneName"]["value"] IS NOT NULL then [1] else [] end |
            MERGE (g:Gene {{ Name:row["geneName"]["value"]}})
            MERGE (g)-[:ENCODES]->(p)
        )

        }} IN TRANSACTIONS OF 100 rows
        """

# ********* Build the entire MetaNetX DB as a graph under our format ********* #


def add_chem_xref(tx, filename):
    """
    A CYPHER query that loads the `chem_xref.tsv` file availaible at the MetaNetX site, using a graph format.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: For performance, it is recommended to split the file in 1 subfile for each row in the DataBase
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line FIELDTERMINATOR '\t'

        MERGE (n:Metabolite {{ MetaNetX_ID:line.ID }})
        SET n.Description = line.description

        WITH toLower(split(line["#source"], ":")[0]) as databasename,
        split(line["#source"], ":")[-1] as databaseid,
        n

        { get_identifiers(from_sparql = False) }
        """)

def add_chem_prop(tx, filename):
    """
    A CYPHER query that loads the `chem_prop.tsv` file availaible at the MetaNetX site, using a graph format.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: For performance, it is recommended to split the file in 1 subfile for each row in the DataBase
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line FIELDTERMINATOR '\t'

        MERGE (n:Metabolite {{ MetaNetX_ID:line["#ID"] }})
        SET n.Name = line.name, n.Formula = line.Formula, n.Charge = line.charge, n.Average_Mass = line.mass,
            n.InChI = line.InChI, n.InChIKey = line.InChIKey, n.SMILES = line.SMILES
        """)

def add_chem_isom(tx, filename):
    """
    A CYPHER query that loads the `chem_isom.tsv` file availaible at the MetaNetX site, using a graph format.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: For performance, it is recommended to split the file in 1 subfile for each row in the DataBase
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line FIELDTERMINATOR '\t'

        MERGE (n:Metabolite {{ MetaNetX_ID:line["#parent"] }})
        MERGE (m:Metabolite {{ MetaNetX_ID:line["child"] }})

        MERGE (m)-[r:ISOMER_OF]-(m)

        SET n.Alternative_names = split(line.description," -> ")[0] + "," + n.Alternative_names
        SET m.Alternative_names = split(line.description," -> ")[1] + "," + m.Alternative_names
        """)

def add_comp_xref(tx, filename):
    """
    A CYPHER query that loads the `comp_xref.tsv` file availaible at the MetaNetX site, using a graph format.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: For performance, it is recommended to split the file in 1 subfile for each row in the DataBase

    .. NOTE:: Some identifiers present the CL/cl prefix. Since I could not find what this prefix refers to,
          and since it only pertains to one single MetaNetX ID, we did not take them into account

    .. NOTE:: The "description" field in the DataBase is ignored, since it seems to be quite similar, but less useful,
          than the "name" field from comp_prop, which is more coherent with our pre-existing schema
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line FIELDTERMINATOR '\t'

        MERGE (n:CelularLocation {{ MetaNetX_ID:line.ID }})

        WITH toLower(split(line["#source"], ":")[0]) as databasename,
        split(line["#source"], ":")[-1] as databaseid,
        n

        FOREACH(ignoreme in case when databasename = "biggc" OR databasename = "bigg.compartment" then [1] else [] end |
            SET n.BiGG_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "go" then [1] else [] end |
            SET n.GO_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "seedc" then [1] else [] end |
            SET n.ModelSeed_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "cl" then [1] else [] end |
            SET n.ModelSeed_ID = databaseid
        )

        FOREACH(ignoreme in case when databasename = "cco" then [1] else [] end |
            SET n.Cell_Component_Ontology_ID = databaseid
        )
        """)

def add_comp_prop(tx, filename):
    """
    A CYPHER query that loads the `comp_prop.tsv` file availaible at the MetaNetX site, using a graph format.

    Args:
        tx          (neo4j.work.simple.Session): The session under which the driver is running
        filename    (str): The name of the CSV file that is being imported

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.

    .. NOTE:: For performance, it is recommended to split the file in 1 subfile for each row in the DataBase
    """
    return tx.run(f"""
        LOAD CSV WITH HEADERS FROM ('file:///{filename}') AS line FIELDTERMINATOR '\t'

        MERGE (n:CelularLocation {{ MetaNetX_ID:line["#ID"] }})
        SET n.Name = line.name
        """)

def add_pept():
    """
    A CYPHER query that all the protein availaible at the MetaNetX site, using a graph format and SPARQL.

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`

    .. NOTE:: SPARQL was only used here because, unlike with the other files, there is no download available;
          also, given there are few proteins, Neo4J is able to process it without running out of memory
          (unlike what happened with the other fields)

    .. NOTE:: This is an **autocommit transaction**. This means that, in order to not keep data in memory
          (and make running it with a huge amount of data) more efficient, you will need to add ```:auto ```
          when calling it from the Neo4J browser, or call it as ```session.run( clean_database() )``` from the driver.
    """
    return f"""
        WITH '

        { add_prefixes() }

        SELECT DISTINCT  *
        WHERE{{
                    ?protein    a               mnx:PEPT ;
                                mnx:peptXref    ?cross_refference .
             }}' AS sparql

        CALL apoc.load.jsonParams(
            replace("https://rdf.metanetx.org/sparql/?query=" + sparql, "\n", ""),
            {{ Accept: "application/sparql-results+json" }}, null )
        YIELD value

        UNWIND value['results']['bindings'] as row

        CALL {{
            WITH row
            MERGE (p:Metabolite {{MetaNetX_ID:split(row['protein']['value'],'/')[-1]}})
            SET p:Protein, p.UniProt_ID = split(row['cross_refference']['value'],'/')[-1]
        }} IN TRANSACTIONS OF 1000 rows
        """

# ********* Annotate existing nodes using KEGG Pathways and Component IDs ********* #

def get_kegg_pathways_for_metabolites():
    """
    A function that finds the Pathways a given Metabolite (based on its Kegg_ID) is a part of, using KEGG.
    This uses genome.jp's dbget web service, since I honestly could not find a way to use KEGG's SPARQL
    service (https://www.genome.jp/linkdb/linkdb_rdf.html) for that.

    .. seealso:: Another possibility could be using `Kegg's Rest API <https://rest.kegg.jp/get/R00703/>`_

    Returns:
        str: A text chain that represents the CYPHER query with the desired output. This can be run using: :obj:`neo4j.Session.run`
    """
    return """
        CALL {

        MATCH (n)
        WHERE (n:Metabolite OR n:Protein OR n:Drug) AND n.KEGG_ID IS NOT null AND n.KEGG_ID <> ""

        LOAD CSV FROM "https://www.genome.jp/dbget-bin/get_linkdb?targettype=all&keywords=cpd:" +  n.KEGG_ID + "&targetformat=text" AS line FIELDTERMINATOR '\t'
        WITH
                split(line[1],":")[0] AS databasename,
                split(line[1],":")[1] AS databaseid,
                n
        FOREACH(ignoreme in case when databasename = "path" then [1] else [] end |
                MERGE (p:Pathway {KEGG_ID:databaseid})
                MERGE (n)-[r:PART_OF_PATHWAY]-(p)
        )

        } IN TRANSACTIONS OF 100 rows
        """

# ********* Build from file ********* #

def build_from_file(filename, driver):
    """
    A function able to build a portion of the MetaNetX database in graph format, provided that one MetaNetX CSV is supplied to it.
    This CSVs are downloaded from the website, and can be presented either as the full file, or as a splitted
    version of it, with just one item per file (which is recommended due to memory limitations). If you want all the database to be
    imported, you should run this function with all the CSVs that form it, as portrayed in the :obj:`~CanGraph.MeSHandMetaNetX.main` module

    Args:
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        filename    (str): The name of the CSV file that is being imported

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    if "chem_xref" in filename:
        with driver.session() as session:
            session.write_transaction(add_chem_xref, filename)
    elif "chem_prop" in filename:
        with driver.session() as session:
            session.write_transaction(add_chem_prop, filename)
    elif "chem_isom" in filename:
        with driver.session() as session:
            session.write_transaction(add_chem_isom, filename)
    elif "comp_xref" in filename:
        with driver.session() as session:
            session.write_transaction(add_comp_xref, filename)
    elif "comp_prop" in filename:
        with driver.session() as session:
            session.write_transaction(add_comp_prop, filename)
