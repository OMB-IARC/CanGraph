Avisar que se instalen APOC y NeoSemantics
CALL db.schema.visualization()
SOURCE: https://towardsdatascience.com/lord-of-the-wiki-ring-importing-wikidata-into-neo4j-and-analyzing-family-trees-da27f64d675e

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

    URL https://github.com/neo4j-contrib/neo4j-apoc-procedures/issues/2377
