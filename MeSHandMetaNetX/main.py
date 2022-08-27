#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that leverages the functions present in the :obj:`~CanGraph.MeSHandMetaNetX.build\_database`
module to recreate `the MetaNetX database <https://www.metanetx.org/>`_ using a graph format and Neo4J,
and then provides an GraphML export file. It also annotates related MeSH_IDs and KEGG Pathway IDs

Please note that, to work, the functions here pre-suppose you have internet access, which will be used to download
MetaNetX's TSVs under a folder provided as ```sys.argv[4]```. (please ensure you have read-write access there)
and query some web SPARQL and REST web services.

For more details on how to run this script, please consult the package's README
"""

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import os, sys, shutil               # Vital modules to interact with the filesystem

# Import subscripts for the program
import build_database
# A hack that allows us to de-duplicate the miscleaneous script in this less-used script
sys.path.append("../")
# NOTE: Please beware that, if using this module by itself, you might need to copy "miscelaneous.py" into your path
# This is not the most elegant, but simplifies code maintenance, and this script shouldnt be used much so...
import miscelaneous as misc

def main():
    """
    The function that executes the code
    """

    mnx_urls = ["https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",
                "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv",
                "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_isom.tsv",
                "https://www.metanetx.org/cgi-bin/mnxget/mnxref/comp_xref.tsv",
                "https://www.metanetx.org/cgi-bin/mnxget/mnxref/comp_prop.tsv"
                ]

    instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
    driver = GraphDatabase.driver(instance, auth=(user, passwd))

    Neo4JImportPath = misc.get_import_path(driver)

    print("Connected to Neo4J")

    #with driver.session() as session:
        #session.run( misc.clean_database() )

    print("Cleaned DataBase")

    datafolder = os.path.abspath(sys.argv[4])

    # We download all the necessary DB files
    #for index, url in enumerate(mnx_urls):
        #print(f"Downloading files: {index + 1}/{len(mnx_urls)}")
        #misc.download(url, f"{datafolder}")
        #print(f"Splitting files: {index + 1}/{len(mnx_urls)}")
        #misc.split_csv(f"{url.split('/')[-1]}", f"{datafolder}", sep='\t', sep_out='\t', startFrom=351, withStepsOf=1000)

    all_files = []
    for root,dirs,files in os.walk(sys.argv[4]):
        for filename in files:
            all_files.append( os.path.abspath(os.path.join(root, filename)) )

    print("Database ready for import. Commencing process...")
    with alive_bar(len(all_files) + 26) as bar:
        # First, we import all the files
        for filename in os.listdir(datafolder):
            #shutil.copyfile(f"{os.path.abspath(sys.argv[4])}/{filename}", f"{Neo4JImportPath}/{filename}")
            #build_database.build_from_file(filename, driver)
            #os.remove(f"{Neo4JImportPath}/{filename}")
            bar()

        # Then, we add proteins and their properties:
        with driver.session() as session:
            #session.run(build_database.add_pept())
            bar(); bar(); bar(); bar(); bar() # Add some bulk bars because this will take forever
            #session.run(build_database.find_protein_data_in_metanetx())
            bar(); bar(); bar(); bar(); bar()

        # And their interactions with other MetaNetX metabolites
        with driver.session() as session:
            session.run(build_database.find_protein_interactions_in_metanetx())
            bar(); bar(); bar(); bar(); bar()

        # Finally, we try to find KEGG Patways for all metabolites in the DB
        with driver.session() as session:
            #session.run(build_database.get_kegg_pathways_for_metabolites())
            bar(); bar(); bar(); bar(); bar()

        # And MeSH IDs, too (by name):
        with driver.session() as session:
            #session.run(build_database.add_mesh_by_name())
            bar(); bar(); bar(); bar(); bar()

        # At the end, purge the database
        #misc.purge_database(driver)

        # And export it:
        with driver.session() as session:
            session.write_transaction(misc.export_graphml, "graph.graphml")
            bar()

    print(f"You can find the exported graph at {Neo4JImportPath}/graph.graphml")
    shutil.copyfile(f"{Neo4JImportPath}/graph.graphml", f"./graph.graphml")
    print(f"A copy of the file has been saved in this project's work directory")

if __name__ == '__main__':

    main()
