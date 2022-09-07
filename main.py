#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that leverages the functions present in the :obj:`~CanGraph.miscelaneous`
module and all other subpackages to annotate metabolites using a graph format and Neo4J,
and then provides an GraphML export file.

CanGraph.main Usage
---------------------

To use this module:

.. argparse::
   :module: CanGraph.main
   :func: args_parser
   :prog: python3 main.py
   :nodefault:

You may find more info in the package's README.

.. NOTE:: For this program to work, the Git environment **has to be set up first**.
    You can ensure this by using: :obj:`CanGraph.setup.setup_git`

CanGraph.main Functions
-------------------------

This module is comprised of:
"""

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
from alive_progress import alive_bar # A cute progress bar that shows the script is still running
import rdkit                         # Cheminformatics and ML package
import pandas as pd                  # Analysis of tabular data
import os, sys, shutil               # Vital modules to interact with the filesystem
import rdkit                         # Cheminformatics and ML package
from rdkit.Chem import MACCSkeys     # MACCS fingerprint calculation
from Bio import SeqIO                # Bioinformatics package
import re                            # Regular expression search
import argparse                      # Arguments pàrser for Python

# Import internal modules for the program
import miscelaneous as misc
from GraphifyDrugBank import build_database as DrugBankDataBase
from GraphifyHMDB import build_database as HumanMetabolomeDataBase
from GraphifySMPDB import build_database as SmallMoleculePathWayDataBase
from ExposomeExplorer import build_database as ExposomeExplorerDataBase
from QueryWikidata import build_database as WikiDataBase
from MeSHandMetaNetX import build_database as MeSHandMetaNetXDataBases

def args_parser():
    """
    Parses the command line arguments into a more usable form, providing help and more

    Returns:
        argparse.ArgumentParser:
            A dictionary of the different possible options for the program as keys, specifying their set value.
            If no command-line arguments are provided, the help message is shown and the program exits.

    .. NOTE:: Note that, in Google Docstrings, if you want a multi-line ``Returns`` comment,
        you have to start it in a different line :(
    .. NOTE:: The return **must** be of type :obj:`argparse.ArgumentParser` for the ``argparse``
        directive to work and auto-gen docs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("neo4jadress", help="the URL of the database, in neo4j:// or bolt:// format")
    parser.add_argument("username", help="the username of the neo4j database in use")
    parser.add_argument("password", help="the passowrd for the neo4j database in use. NOTE: "
                                         "Since passed through bash, you may need to escape special characters")
    parser.add_argument("databasefolder", help="The folder indicated to ```setup.py``` as the one where your databases "
                                               "will be stored")
    parser.add_argument("inputfile", help="The location of the CSV file in which the program will search for metabolites")

    # If no args are provided, show the help message
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser

def scan_folder(folder_path):
    """
    Scans a folder and finds all the files present in it

    Args:
        folder_path (str): The folder that is to be scanned
    Returns:
        list: A list of all the files in the folder, listed by their absolute path
    """
    all_files = []
    for root,dirs,files in os.walk(folder_path):
        for filename in files:
            all_files.append( os.path.abspath(os.path.join(root, filename)) )
    return all_files

def find_reasons_to_import(filepath, chebi_ids, names, hmdb_ids, inchis):
    """
    Finds reasons to import a metabolite given a candidate filepath **with one metabolite per file**
    and a series of lists containing all synonyms of the values considered reasons for import

    Args:
        filepath (str): The path to the file in which we will search for reasons to import
        chebi_ids (list): A list of all the ChEBI_IDs which are considered a reason to import
        names (list): A list of all the Names which are considered a reason to import
        hmdb_ids (list): A list of all the HMDB_IDs which are considered a reason to import
        inchis (list): A list of all the InChIs which are considered a reason to import

    Returns:
        list: A list of the methods that turned out to be valid for import, such as Name, ChEBI_ID...
    """
    import_based_on = []; text = ""
    relpath = os.path.relpath(filepath, ".")

    with open(f'{filepath}', "r") as f:
        text = f.read()

    # We try to find exact InChI matches:
    if any(inchi in text for inchi in inchis):
        import_based_on.append("Exact InChI")
    # If none if found, we use the "Similarity Evaluator" metric
    elif "InChI=" in text:
        for inchi in inchis:
            # SOURCE: https://chemistry.stackexchange.com/questions/82144/what-is-the-correct-regular-expression-for-inchi
            results = re.search("InChI\=1S?\/[A-Za-z0-9\.]+(\+[0-9]+)?(\/[cnpqbtmsih][A-Za-z0-9\-\+\(\)\,\/\?\;\.]+)*(\"|\<)", text)
            if results:
                result = results.group(0).replace("<","").replace("\"","")

                found_error = False
                try:
                    Query = rdkit.Chem.MolFromInchi(f"{result}")
                    MACCSQuery = MACCSkeys.GenMACCSKeys(Query)
                except Exception as error:
                    found_error = True
                    pass

                if found_error == False:
                    Subject = rdkit.Chem.MolFromInchi(inchi)
                    MACCSSubject = MACCSkeys.GenMACCSKeys(Subject)


                    DICE_MACCS = rdkit.DataStructs.DiceSimilarity(MACCSQuery, MACCSSubject)
                    if DICE_MACCS > 0.95:
                        import_based_on.append(f"DICE-MACCS {100*round(DICE_MACCS, 4)} % similarity")

    # For CHEBI, if we are using E-E, and since they dont have a prefix (i.e. they are only a number) we have to process the files.
    if "ExposomeExplorer/components" in relpath:
        component = pd.read_csv(os.path.abspath(filepath))
        chebi_query = component["chebi_id"]
        # NOTE: Here, we remove the optional CHEBI: prefix
        if chebi_query in [chebi_id.replace("CHEBI:", "").replace("chebi:", "") for chebi_id in chebi_ids]:
            import_based_on.append("ChEBI")
    # And, even if its not E-E, we still need to add the tag before for things to match
    for chebi_query in chebi_ids:
        if f"<chebi_id>{chebi_query.replace('CHEBI:', '')}" in text:
            import_based_on.append("ChEBI")

    # For the rest of the databases, we simply search for exact matches in our list and the texts:
    if any(hmdb in text for hmdb in hmdb_ids):
        import_based_on.append("HMDB ID")
    if any(name in text for name in names):
        import_based_on.append("Name")

    # We return a list of a dict from keys to remove duplicates from the "reasons to import" list
    return list(dict.fromkeys(import_based_on))

def build_from_file(filepath):
    """
    Imports a given metabolite from a sigle-metabolite containing file by checking its type
    and calling the appropriate import functions.

    Args:
        filepath (str): The path to the file in which will be imported

    Returns:
        This function does not provide a particular return, but rather imports the requested file

    .. NOTE:: The ``filepath`` may be absolute or relative, but it is transformed to a relative ``relpath``
        in order to remove possible influence of higher-name folders in the import type selection. This is
        also why the condition is stated as a big "if/elif/else" instead of a series of "ifs"
    """
    relpath = os.path.relpath(filepath, ".")

    if "DrugBank" in relpath:
        shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
        DrugBankDataBase.build_from_file(f"{os.path.basename(filepath)}", driver)
        os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")

    elif "HMDB" in relpath:
        if "protein" in relpath:
            shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
            HumanMetabolomeDataBase.build_from_protein_file(f"{os.path.basename(filepath)}", driver)
            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")
        elif "metabolite" in relpath:
            shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
            HumanMetabolomeDataBase.build_from_metabolite_file(f"{os.path.basename(filepath)}", driver)
            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")

    elif "SMPDB" in relpath:
        # NOTE: Since this adds a ton of low-resolution nodes, maybe have this db run first?
        # We will ignore the smpdb_pathways file because it doesnt have "real" identifiers
        if "proteins" in relpath:
            SmallMoleculePathWayDataBase.build_from_file(sys.argv[4], filepath, Neo4JImportPath, driver, "Protein")
        if "metabolites" in relpath:
            SmallMoleculePathWayDataBase.build_from_file(sys.argv[4], filepath, Neo4JImportPath, driver, "Metabolite")

    elif "ExposomeExplorer/components" in relpath:
            # NOTE: Since only "components" can result in a match based on our current criteria,
            #   we will build the DB starting with the components only.
            # Here, instead of using shutil.copyfile, we will use pandas to purge the _count columns when copying
            original_file = pd.read_csv(filepath)[[x for x in open(f"{filepath}").readline().rstrip().split(",")
                                                if not x.endswith('_count')]]
            original_file.to_csv(f"{Neo4JImportPath}/{os.path.basename(filepath)}", index=False)
            with driver.session() as session:
                    session.write_transaction(ExposomeExplorerDataBase.add_components, os.path.basename(filepath))
            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")

            ExposomeExplorerDataBase.build_from_file( os.path.dirname(filepath), Neo4JImportPath, driver, False)


def link_to_original_data(driver, original_ids, import_based_on):
    """
    Links a recently-imported metabolite to the original data (that which caused it to be imported) by creating an
    ``ÒriginalMetabolite`` node that is ``(n)-[r:ORIGINALLY_IDENTIFIED_AS]->(a)`` related to the imported data

    Args:
         driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
         original_ids (dict): A dictionary of the Original IDs in our query file, which was the basis for import
         import_based_on (list): A list of the methods that turned out to be valid for import, such as Name, ChEBI_ID...

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.
    """
    with driver.session() as session:
        session.run( f"""
                    MATCH (a) WHERE a.InChI = "{row["InChI"]}"
                    MATCH (c) WHERE c.Name = "{row["Name"]}"
                    MATCH (d) WHERE d.SMILES = "{row["SMILES"]}"
                    MATCH (e) WHERE e.InChI = "{row["InChI"]}"
                    MATCH (f) WHERE f.HMDB_ID = "{row["Identifier"]}"
                    MATCH (g) WHERE g.Monisotopic_Molecular_Weight = "{row["MonoisotopicMass"]}"
                    CREATE (n:OriginalMetabolite)
                    SET n.InChI = "{row["InChI"]}", n.Name = "{row["Name"]}", n.ChEBI = "{row["ChEBI"]}",
                        n.SMILES = "{row["SMILES"]}", n.HMDB_ID = "{row["Identifier"]}",
                        n.Monisotopic_Molecular_Weight = "{row["MonoisotopicMass"]}"
                    MERGE (n)-[r1:ORIGINALLY_IDENTIFIED_AS]->(a)
                    MERGE (n)-[r2:ORIGINALLY_IDENTIFIED_AS]->(b)
                    MERGE (n)-[r3:ORIGINALLY_IDENTIFIED_AS]->(c)
                    MERGE (n)-[r4:ORIGINALLY_IDENTIFIED_AS]->(d)
                    MERGE (n)-[r5:ORIGINALLY_IDENTIFIED_AS]->(e)
                    MERGE (n)-[r6:ORIGINALLY_IDENTIFIED_AS]->(f)
                    MERGE (n)-[r7:ORIGINALLY_IDENTIFIED_AS]->(g)
                    SET r1.Identified_By = {import_based_on}
                    SET r2.Identified_By = {import_based_on}
                    SET r3.Identified_By = {import_based_on}
                    SET r4.Identified_By = {import_based_on}
                    SET r5.Identified_By = {import_based_on}
                    SET r6.Identified_By = {import_based_on}
                    SET r7.Identified_By = {import_based_on}
                    """ )

def annotate_using_wikidata(driver):
    """
    Once we finish the search, we annotate the nodes added to the database using WikiData

    Args:
         driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.

    .. TODO:: When fixing queries, fix the main subscript also
    """
    with driver.session() as session:
        misc.repeat_transaction(WikiDataBase.add_wikidata_to_mesh, 10, driver)
        # The ``query`` param is, remember, so as to remove the wikidata_id search which is by default
        misc.repeat_transaction(WikiDataBase.add_metabolite_info, 10, driver, query = "ChEBI_ID")
        misc.repeat_transaction(WikiDataBase.add_drug_external_ids, 10, driver, query = "DrugBank_ID")
        misc.repeat_transaction(WikiDataBase.add_more_drug_info, 10, driver, query = "DrugBank_ID")

        misc.repeat_transaction(WikiDataBase.find_subclass_of_cancer, 10, driver)
        misc.repeat_transaction(WikiDataBase.find_subclass_of_cancer, 10, driver)
        misc.repeat_transaction(WikiDataBase.find_subclass_of_cancer, 10, driver)
        misc.repeat_transaction(WikiDataBase.find_instance_of_cancer, 10, driver)

        # For each of the 10 numbers a wikidata_id may have as ending
        for number in range(10):
            misc.repeat_transaction(WikiDataBase.add_cancer_info, 10, driver, number)
            misc.repeat_transaction(WikiDataBase.add_drugs, 10, driver, number)
            misc.repeat_transaction(WikiDataBase.add_causes, 10, driver, number)
            misc.repeat_transaction(WikiDataBase.add_genes, 10, driver, number)

        misc.repeat_transaction(WikiDataBase.add_drug_external_ids, 10, driver)
        misc.repeat_transaction(WikiDataBase.add_more_drug_info, 10, driver)
        misc.repeat_transaction(WikiDataBase.add_gene_info, 10, driver)
        misc.repeat_transaction(WikiDataBase.add_metabolite_info, 10, driver)

def add_mesh_and_metanetx(driver):
    """
    Add MeSH Term IDs, Synonym relations and Protein interactions to existing nodes using MeSH and MetaNetX
    Also, adds Kegg Pathway IDs

    Args:
         driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    # We will also add MeSH terms to all nodes:
    with driver.session() as session:
        misc.repeat_transaction(MeSHandMetaNetXDataBases.add_mesh_by_name(), 10, driver, bar)
    # We also add synonyms:
    with driver.session() as session:
        session.run(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("Name"))
        session.run(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("KEGG_ID"))
        session.run(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("ChEBI_ID"))
        session.run(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("HMDB_ID"))
        session.run(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("InChI"))
        session.run(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("InChIKey"))

    # And some protein interactions, together with their pathways, too:
    with driver.session() as session:
        session.run(MeSHandMetaNetXDataBases.find_protein_data_in_metanetx())
        session.run(MeSHandMetaNetXDataBases.find_protein_interactions_in_metanetx())
        session.run(MeSHandMetaNetXDataBases.get_kegg_pathways_for_metabolites())

def main():
    """
    The function that executes the code

    .. TODO:: CAMBIAR NOMBRE A LOS MESH PARA INDICAR EL TIPO. AÑADIR NAME A LOS WIKIDATA
    .. TODO:: FIX THE REPEAT TRANSACTION FUNCTION
    .. TODO:: Match partial InChIs based on DICE-MACCS
    .. TODO:: QUE FUNCIONE -> ACTUALMENTE ESTA SECCION RALENTIZA MAZO
    .. TODO:: CHECK APOC IS INSTALLED
    .. TODO:: FIX MAIN
    .. TODO:: MERGE BY INCHI, METANETX ID
    .. TODO:: Fix find_protein_interactions_in_metanetx
    .. TODO:: Mover esa funcion de setup a misc
    .. TODO:: EDIT conf.py
    .. TODO:: Detail Schema Changes

    .. TODO:: Document the following Schema Changes:
        * For Subject, we have a composite PK: Exposome_Explorer_ID, Age, Gender e Information
        * Now, more diseases will have a WikiData_ID and a related MeSH. This will help with networking. And, this diseases dont even need to be a part of a cancer!
        * The Gene nodes no longer exist in the full db? -> They do
    """

    # Parse the command line arguments
    # Done first in order to show errors if bad commands are issued
    parser = args_parser(); args = vars(parser.parse_args())

    # First, we prepare a scan of all the files available on our "DataBases" folder
    # We will cycle through them later on to try and find matches
    all_files = scan_folder(args["databasefolder"])

    raw_database = pd.read_csv(sys.argv[5], delimiter=',', header=0)

    # Set up the progress bar
    with alive_bar(len(all_files)*len(raw_database)) as bar:

        # And connect to the Neo4J database
        instance = args["neo4jadresss"]; user = args["username"]; passwd = args["password"]
        driver = GraphDatabase.driver(instance, auth=(user, passwd))

        Neo4JImportPath = misc.get_import_path(driver)
        print("Connected to Neo4J")

        # For each item in our "query" file, we will try to find matches:
        for index, row in raw_database.iterrows():

            # We start by cleaning the database (important if this is not the first run)
            with driver.session() as session:
                session.run( misc.clean_database() )
                if index == 0: print("Cleaned DataBase") # Only the first time so as not to repeat

            print(f"Importing metabolites: {index+1}/{len(raw_database)}")

            # And then, we can calculate synonyms for each metabolite we will search
            names = misc.find_synonyms(driver, "Name", row.Name)
            inchis = misc.find_synonyms(driver, "InChI", row.InChI)
            hmdb_ids = misc.find_synonyms(driver, "HMDB_ID", row.Identifier)
            chebi_ids = misc.find_synonyms(driver, "ChEBI_ID", row.ChEBI.replace("CHEBI:", ""))

             # And search for them in the all_files list we created earlier on based on a series of criteria:
            for filepath in all_files:
                import_based_on = find_reasons_to_import(filepath, chebi_ids, names, hmdb_ids, inchis)

                # Once we know the reasons to import (this is done so that it only cycles one
                # time through the code), we import the files themselves
                if len(import_based_on) != 0:
                    build_from_file(filepath)

                    #Each time we import some nodes, link them to our original data
                    #NOTE: This WILL duplicate nodes and relations, but we will fix this later when we purge the DB
                    original_ids = row.to_dict('records')[0] # Convert to dict first for the functions
                    link_to_original_data(driver, original_ids, import_based_on)

                # And advance, of course
                bar()

            # Finally, we purge the database by removing nodes considered as duplicated
            # We only purge once at the end, to limit processing time
            misc.purge_database(driver)

            # And save it in GraphML format
            with driver.session() as session:
                session.write_transaction(misc.export_graphml, f"metabolite_{index+1}.graphml")

            print(f"Metabolite {index+1}/{len(raw_database)} processed. You can find a copy of the associated knowledge graph at {os.path.abspath(sys.argv[4])}/metabolite_{index+1}.graphml")
            shutil.copyfile(f"{Neo4JImportPath}/metabolite_{index+1}.graphml", f"{os.path.abspath(sys.argv[4])}/metabolite_{index+1}.graphml")

if __name__ == '__main__':

    main()
