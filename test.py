from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import rdkit                         # Cheminformatics and ML package
import pandas as pd                  # Analysis of tabular data
import os, sys, shutil               # Vital modules to interact with the filesystem
import rdkit                         # Cheminformatics and ML package
from rdkit.Chem import MACCSkeys     # MACCS fingerprint calculation
from Bio import SeqIO                # Bioinformatics package
import re                            # Regular expression search
import miscelaneous as misc

instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
driver = GraphDatabase.driver(instance, auth=(user, passwd))

Neo4JImportPath = misc.get_import_path(driver)
print("Connected to Neo4J")

def test(tx):
    print( tx.run("""MATCH (n) RETURN n LIMIT 1""") )

with driver.session() as session:
    print(driver)
    session.write_transaction(test)


