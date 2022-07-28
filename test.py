from neo4j import GraphDatabase      # The Neo4J python driver
import sys
import miscelaneous as misc

def trial(tx, title):
    return tx.run("""

        MATCH (p:Person)-[:ACTED_IN]->(m:Movie)
        WHERE m.title = $title // (1)
        RETURN p.name AS a

        """, title=title)

instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
driver = GraphDatabase.driver(instance, auth=(user, passwd))

Neo4JImportPath = misc.get_import_path(driver)

with driver.session() as session:
    result = session.write_transaction(trial, "The Matrix")
    for record in result:
        print("{} knows {}".format(record["a"] ,record["a"]))
