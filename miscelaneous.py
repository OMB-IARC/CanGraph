##!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Import modules necessary for this subscript

from urllib.request import urlopen   # To open a URL in python
from io import BytesIO               # To process the opened file
from zipfile import ZipFile          # To unzip things from Exposome Explorer

def get_neo4j_path(tx):
    result= tx.run("Call dbms.listConfig() YIELD name, value WHERE name='dbms.directories.neo4j_home' RETURN value")
    return [record["value"] for record in result]

def get_import_path(current_driver):
    with current_driver.session() as session:
        Neo4JImportPath = session.read_transaction(get_neo4j_path)[0]+'/import'
    return Neo4JImportPath

def download_and_unzip(url, path):
    http_response = urlopen(url) # Open the URL
    zipfile = ZipFile(BytesIO(http_response.read())) # Read its contents
    zipfile.extractall(path=path) # Extract CSV

