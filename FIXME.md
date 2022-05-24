Avisar que se instalen APOC y NeoSemantics
CALL db.schema.visualization()
SOURCE: https://towardsdatascience.com/lord-of-the-wiki-ring-importing-wikidata-into-neo4j-and-analyzing-family-trees-da27f64d675e

PRINT("LEETE LAS CONVENCIONES DE NOMBRADO DE NEO4J PARA EL SCRIPT FINAL")

// Prepare a SPARQL query
WITH 'SELECT ?item ?itemLabel
WHERE{ ?item wdt:P31 wd:Q12078 . SERVICE wikibase:label {bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en" }}' AS sparql
// make a request to Wikidata
CALL apoc.load.jsonParams(
        replace("https://query.wikidata.org/sparql?query=" + sparql, "\n", ""),
        { Accept: "application/sparql-results+json"}, null)
YIELD value
// Unwind results to row
UNWIND value['results']['bindings'] as row
// Prepare data
WITH row['itemLabel']['value'] as race,
     row['item']['value'] as url,
     split(row['item']['value'],'/')[-1] as id
// Store to Neo4j
CREATE (r:Race)
SET r.race = race,
    r.url = url,
    r.id = id
FIX REQUIREMENTS: POR EJEMPLO, AÑADIR BIO. BASICALLY LEETE TODO Y LO VAS METIENDO
    URL https://github.com/neo4j-contrib/neo4j-apoc-procedures/issues/2377

    Maybe all sequences should be in a "SEQUENCED_AS"?
purge_database(tx)
HACER MAS USO DEL SESSION.RUN PARA AUTOCOMMIT COSILLAS PEQUEÑAS, COMO PURGE DATABASES
ejemplo:

def purge_database(tx):
    """
    This function purges and prepares the database for export:
        * Removes "Publication" nodes without the Pubmed_ID property
        * For relations with a "Pubmed_ID" property, it deletes the last character. This is an unnecesary ",".
        * In case any "Concentration" node was created without measurements, it is removed.
    WARNING: This should be run ONLY ONCE PER SCRIPT EXECUTION, or else the r.Pubmed_ID property will be affected
    """
    return tx.run(f"""
        MATCH (p:Publication)
            WHERE p.Pubmed_ID IS null
        DETACH DELETE p
        WITH p
        MATCH ()-[r]-() SET r.Pubmed_ID = substring(r.Pubmed_ID, 0, size(r.Pubmed_ID) -1 )
        WITH p
        MATCH (c:Concentration)
        WHERE size(keys(properties(c))) < 2
        DETACH DELETE c
        """)

        NO deberia haber nodos blancos ya ue invocamos CREATE con cabeza. Still pasa un session.run
