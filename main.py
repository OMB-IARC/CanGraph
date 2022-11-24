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
from contextlib import redirect_stdout # Redirect stdout, to not show things on the stdout
from copy import deepcopy            # Do deep copies of python objects
import bioservices                   # Query web bio-databases from python
import ijson                         # Read JSON files from Python in an iterative way
import logging                       # Make ``verbose`` messages easier to show

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
    .. NOTE:: By using :obj:`argparse.const` instead of :obj:`argparse.default`, the check_file function will check ""
        (the current dir, always exists) if the arg is not provided, not breaking the function; if it is, it checks it.
    """
    parser = argparse.ArgumentParser(
        description = "A python utility to study and analyse cancer-associated metabolites "
                      "using knowledge graphs ")
    parser.add_argument("-c", "--check_args", action="store_true",
                        help="Checks if the rest of the arguments are OK, then exits")
    parser.add_argument("-n", "--noindex", action="store_true",
                        help="Runs the program checking each file one-by-one, instead of using a JSON index")
    parser.add_argument("-s", "--similarity", action="store_true", default=True,
                        help="Deactivates the import of information based on Structural Similarity."
                             "This might dramatically increase processing time; default is True.")
    parser.add_argument("-w", "--webdbs", action="store_true", default=True,
                        help="Activates import of information based on web databases."
                             "This might dramatically increase processing time; default is True.")
    parser.add_argument("-i", "--interactive", action="store_true", default=True,
                    help=("tells the script if it wants interaction from the user "
                            "and more information shown to them; similar to --verbose"))

    parser.add_argument("--query", type=misc.check_file, required=True,
                        help="The location of the CSV file in which the program will search for metabolites")
    parser.add_argument("--dbfolder", type=misc.check_file, default="DataBases",
                        help="The folder indicated to ```setup.py``` as the one where your databases will be stored; "
                             "default is ``./DataBases``")
    parser.add_argument("--results", default="Results",
                        help="The folder where the resulting GraphML exports will be stored; "
                             "default is ``./Results``")
    parser.add_argument("--adress", help="the URL of the database, in neo4j:// or bolt:// format",
                        type=misc.check_neo4j_protocol, default="bolt://localhost:7687")
    parser.add_argument("--username", help="the username of the neo4j database in use",
                        default="neo4j")
    parser.add_argument("--password", help="the password for the neo4j database in use. NOTE: "
                                         "Since passed through bash, you may need to escape some chars",
                        default="neo4j")

    # If no args are provided, show the help message
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if "--check_args" in sys.argv:
        parser.parse_args()
        exit(0)

    return parser

def improve_search_terms_with_metanetx(query, query_type,
                                       driver, chebi_ids, names, hmdb_ids, inchis, mesh_ids):
    """
    Improves the search terms already provided to the CanGraph programme
    by using the MetaNetX web service to find synonyms in IDs

    Args:
        query (str): The term we are currently querying for
        query_type (str): The kind of query to search; one of ["ChEBI_ID", "HMDB_ID", "Name", "InChI", "MeSH_ID"]
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        chebi_ids (str): A string of ";" separated values of all the ChEBI_ID representing the current metabolite
        names (list): A string of ";" separated values of all the Name representing the current metabolite
        hmdb_ids (list): A string of ";" separated values of all the HMDB_ID representing the current metabolite
        inchis (list): A string of ";" separated values of all the InChI representing the current metabolite
        mesh_ids (list): A string of ";" separated values of all the MeSH_ID representing the current metabolite

    Returns:
        list: A list containing [ chebi_ids, names, hmdb_ids, inchis, mesh_ids ], with all their synonyms
    """
    # First in MetaNetX
    with driver.session() as session:
        # For MeSH, we cannot search for synonyms in MetaNetX (it doesn't index them), so we instead
        # look for related metabolites in MeSH itself, so that the import makes more sense
        if query_type == "MeSH_ID":
            graph_response = session.execute_read(
                MeSHandMetaNetXDataBases.find_metabolites_related_to_mesh, query)
        else:
            graph_response = session.execute_read(
                MeSHandMetaNetXDataBases.read_synonyms_in_metanetx, query_type, query)

    for element in graph_response:
        element = {key: element[key] for key in element if element[key] != None }
        # With .get, the lookup does not fail, even if query_type == "MeSH_ID"
        if element.get("databasename", "").lower() == "hmdb":
            if element["databaseid"] not in hmdb_ids: hmdb_ids.append(element.get("databaseid", ""))

        if element.get("databasename", "").lower() == "chebi":
            if element["databaseid"] not in chebi_ids: chebi_ids.append(element.get("databaseid", ""))

        if element.get("InChI", "").lower() not in inchis: inchis.append(element.get("InChI", ""))

        if element.get("Name", "").lower() not in names: names.append(element.get("Name", ""))

    return [ chebi_ids, names, hmdb_ids, inchis, mesh_ids ]

def improve_search_terms_with_cts(query, query_type,
                                  chebi_ids, names, hmdb_ids, inchis, mesh_ids):
    """
    Improves the search terms already provided to the CanGraph programme
    by using The Chemical Translation Service to find synonyms in IDs

    Args:
        query (str): The term we are currently querying for
        query_type (str): The kind of query to search; one of ["ChEBI_ID", "HMDB_ID", "Name", "InChI", "MeSH_ID"]
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        chebi_ids (str): A string of ";" separated values of all the ChEBI_ID representing the current metabolite
        names (list): A string of ";" separated values of all the Name representing the current metabolite
        hmdb_ids (list): A string of ";" separated values of all the HMDB_ID representing the current metabolite
        inchis (list): A string of ";" separated values of all the InChI representing the current metabolite
        mesh_ids (list): A string of ";" separated values of all the MeSH_ID representing the current metabolite

    Returns:
        list: A list containing [ chebi_ids, names, hmdb_ids, inchis, mesh_ids ], with all their synonyms
    """
    # We get the correct names for the identifiers we want to turn into:
    query_dict_keys = ["ChEBI_ID", "HMDB_ID", "Name", "InChI", "MeSH_ID"]
    oldIdentifier = [ x for x in query_dict_keys if x not in [query_type, "MeSH_ID"] ]

    replace_keys = {"ChEBI_ID":"ChEBI", "HMDB_ID":"Human Metabolome Database", "Name":"Chemical Name", "InChI":"InChIKey"}
    toIdentifierList = list(map(replace_keys.get, oldIdentifier, oldIdentifier))
    fromIdentifier = replace_keys.get(query_type, query_type)

    # Convert the InChI to InChIKeys (required by CTS)
    if query_type == "InChI":
        rdkit_mol = rdkit.Chem.MolFromInchi(query)
        searchTerm = rdkit.Chem.MolToInchiKey(rdkit_mol)
    else:
        searchTerm = query

    # And run the CTS search:
    for toIdentifier in toIdentifierList:
        result = MeSHandMetaNetXDataBases.find_synonyms_in_cts(fromIdentifier, toIdentifier, searchTerm)

        for element in result: # And append its results to the correct list
            if toIdentifier == "Human Metabolome Database" and element not in hmdb_ids:
                hmdb_ids.append(element)
            if toIdentifier == "ChEBI" and element.replace("CHEBI:", "") not in chebi_ids:
                chebi_ids.append(element.replace("CHEBI:", ""))
            if toIdentifier == "Chemical Name" and element not in names:
                names.append(element)
            if toIdentifier == "InChIKey": # Convert key to InChI
                unichem = bioservices.UniChem()
                inchis_from_unichem = unichem.get_inchi_from_inchikey(element)
                inchis.append(inchis_from_unichem[0]["standardinchi"])

    return [ chebi_ids, names, hmdb_ids, inchis, mesh_ids ]

def improve_search_terms(driver, chebi_ids, names, hmdb_ids, inchis, mesh_ids):
    """
    Improves the search terms already provided to the CanGraph programme by processing the
    text stings and finding synonyms in various platforms

    Args:
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        chebi_ids (str): A string of ";" separated values of all the ChEBI_ID representing the current metabolite
        names (list): A string of ";" separated values of all the Name representing the current metabolite
        hmdb_ids (list): A string of ";" separated values of all the HMDB_ID representing the current metabolite
        inchis (list): A string of ";" separated values of all the InChI representing the current metabolite
        mesh_ids (list): A string of ";" separated values of all the MeSH_ID representing the current metabolite

    Returns:
        list: A list containing [ chebi_ids, names, hmdb_ids, inchis, mesh_ids ], with all their synonyms
    """
    # First, we convert the strings we have received from the function call into lists, by splitting for ";":
    replace_chebis = re.compile(re.escape('ChEBI:'), re.IGNORECASE);
    replace_meshes = re.compile(re.escape('MeSH:'), re.IGNORECASE)

    chebi_ids = replace_chebis.sub('', chebi_ids).split(";")
    mesh_ids = replace_meshes.sub('', mesh_ids).split(";")
    names = names.split(";");    hmdb_ids = hmdb_ids.split(";");     inchis = inchis.split(";")

    # Then, we may search for synonyms for each field using MetaNetX
    # To do this, we organize a dict so that we can iterate more easily
    query_dict = {"ChEBI_ID": chebi_ids, "HMDB_ID": hmdb_ids, "InChI": inchis, "Name": names, "MeSH_ID": mesh_ids}
    query_dict = deepcopy(query_dict) # Make a deepcopy of the dict so that it doesn't update while on the loop

    # And then, proceed to search for synonyms and append the appropriate results if they are not already there
    with alive_bar( sum( [ len(v) for v in query_dict.values() if any(v) ] ) + 4, title="Finding Synonyms...") as bar:
        for query_type, query_list in query_dict.items():
            for query in query_list:
                if not query.isspace() and query:
                    chebi_ids, names, hmdb_ids, inchis, mesh_ids = (
                    improve_search_terms_with_metanetx(query, query_type, driver,
                                                    chebi_ids, names, hmdb_ids, inchis, mesh_ids))

                bar()

        # Once we have all the synonyms that we could find on MeSHandMetaNetX, we remove all duplicates
        chebi_ids = list(filter(None, set(chebi_ids)));   names = list(filter(None, set(names)));
        inchis = list(filter(None, set(inchis)));         mesh_ids = list(filter(None, set(mesh_ids)))
        hmdb_ids = list(filter(None, set(hmdb_ids)))

        # And re-generate the query_dict
        query_dict = {"ChEBI_ID": chebi_ids, "HMDB_ID": hmdb_ids, "InChI": inchis, "Name": names, "MeSH_ID": mesh_ids}
        query_dict = deepcopy(query_dict) # Make a deepcopy of the dict so that it doesn't update while on the loop

        # And search for even more synonyms on CTS, the Chemical Translation Service!
        # This is done separately because, while MeSHandMetaNetX do not overlap, they do overlap with CTS
        for query_type, query_list in {key: query_dict[key]
                                    for key in query_dict if key != 'MeSH_ID'
                                                            and any(query_dict[key])  }.items():
            for query in query_list:
                chebi_ids, names, hmdb_ids, inchis, mesh_ids = (
                    improve_search_terms_with_cts(query, query_type,
                                            chebi_ids, names, hmdb_ids, inchis, mesh_ids))

            bar()

        # Finally, we remove outdated HMDB IDs
        regex = re.compile(r'^HMDB\d\d\d\d\d$')
        hmdb_ids = [i for i in hmdb_ids if regex.match(str(i))]

    # And return the simplified versions of the lists
    return [ list(filter(None, set(chebi_ids))),    list(filter(None, set(names))),
             list(filter(None, set(hmdb_ids))),     list(filter(None, set(inchis))),
             list(filter(None, set(mesh_ids))) ]

def find_reasons_to_import_inchi(query, subject):
    """
    Takes two chains of text and finds if the ``query`` is present in the ``subject``,
    or if there are molecules common between them with at least 95% similarity

    Args:
        query (str or list): A string or list of strings describing valid InChI(s)
        subject (str): A valid InChI

    Returns:
        dict: A dict with each query as a key and the reason to import it as value, if there is one.

    .. seealso:: This approach was taken from `Chemistry StackExchange #82144
        <https://chemistry.stackexchange.com/questions/82144/what-is-the-correct-regular-expression-for-inchi>`_

    .. NOTE:: Since this is a one-to-one comparison, subject and query can be used interchangeably; however,
        bear in mind that only the query can be provided as a list
    """
    reasons_to_import_inchi = {}
    queries = list(query)

    try:
        MolSubject = rdkit.Chem.MolFromInchi(subject)
        MACCSSubject = MACCSkeys.GenMACCSKeys(MolSubject)
    except:
        return

    for each_query in queries:
        found_error = False
        try:
            Query = rdkit.Chem.MolFromInchi(f"{each_query}")
            MACCSQuery = MACCSkeys.GenMACCSKeys(Query)
        except Exception as error:
            found_error = True
            pass

        if found_error == False:
            DICE_MACCS = rdkit.DataStructs.DiceSimilarity(MACCSQuery, MACCSSubject)
            if DICE_MACCS > 0.95:
                dice_maccs_import = f"DICE-MACCS {100*round(DICE_MACCS, 4)} % similarity"
                reasons_to_import_inchi.setdefault(each_query, dice_maccs_import)

    return reasons_to_import_inchi

def find_reasons_to_import_all_files(filepath, similarity, chebi_ids, names, hmdb_ids, inchis, mesh_ids):
    """
    Finds reasons to import a metabolite given a candidate filepath **with one metabolite per file**
    and a series of lists containing all synonyms of the values considered reasons for import

    Args:
        filepath (str): The path to the file in which we will search for reasons to import
        similarity (bool): Whether to use similarity as a measure to import or not
        chebi_ids (list): A list of all the ChEBI_ID which are considered a reason to import
        names (list): A list of all the Name which are considered a reason to import
        hmdb_ids (list): A list of all the HMDB_ID which are considered a reason to import
        inchis (list): A list of all the InChI which are considered a reason to import
        mesh_ids (list): A list of all the MeSH_ID which are considered a reason to import

    Returns:
        list: A list of the methods that turned out to be valid for import, such as Name, ChEBI_ID...
    """
    import_based_on = []; text = ""
    importing_ids = {}

    with open(f'{filepath}', "r") as f:
        text = f.read()

    # We try to find exact InChI matches:
    if any((match:= inchi) in text for inchi in inchis):
        import_based_on.append("Exact InChI")
        importing_ids.setdefault("InChI", []).append(match)
    # If none if found, we use the "Similarity Evaluator" metric
    elif "InChI=" in text and any(inchis) and similarity:
        results = re.search("InChI\=1S?\/[A-Za-z0-9\.]+(\+[0-9]+)?(\/[cnpqbtmsih][A-Za-z0-9\-\+\(\)\,\/\?\;\.]+)*(\"|\<)", text)
        if results:
            result = results.group(0).replace("<","").replace("\"","")
            reason_to_import = find_reasons_to_import_inchi(inchis, result)
            if reason_to_import:
                import_based_on.append(list(set(reason_to_import.values())))
                importing_ids.setdefault("InChI", []).append(match)

    # For CHEBI, if we are using E-E or SMPDB, and since they dont have a prefix
    # (i.e. they are only a number) we have to process the files.
    if "ExposomeExplorer/components" in filepath or "SMPDB/smpdb_metabolites" in filepath:
        component = pd.read_csv(os.path.abspath(filepath), dtype = str)
        if "ExposomeExplorer/components" in filepath:
            if str(component["chebi_id"][0]) != "nan": chebi_query = list(component["chebi_id"][0])
            else: chebi_query = []
        elif "SMPDB/smpdb_metabolites" in filepath: chebi_query = list(component["ChEBI ID"])

        # NOTE: Here, we remove the optional CHEBI: prefix
        chebi_subject = [chebi_id.replace("CHEBI:", "").replace("chebi:", "") for chebi_id in chebi_ids]

        # And we check the two lists against each other
        if any((match:= item) in chebi_query for item in chebi_subject) and chebi_query:
            import_based_on.append("ChEBI_ID")
            importing_ids.setdefault("ChEBI_ID", []).append(match)
    else:
        # And, even if its not E-E, we still need to add the tag before for things to match
        for chebi_query in chebi_ids:
            replace_tag_exp = re.compile(re.escape('ChEBI:'), re.IGNORECASE)
            if f"<chebi_id>{replace_tag_exp.sub('', chebi_query)}" in text:
                import_based_on.append("ChEBI_ID")
                importing_ids.setdefault("ChEBI_ID", []).append(match)

    # For MeSH, we just check for the tags
    for mesh_id in mesh_ids:
        replace_tag_exp = re.compile(re.escape('MeSH:'), re.IGNORECASE)
        if f"<mesh-id>{replace_tag_exp.sub('', mesh_id)}" in text:
            import_based_on.append("MeSH_ID")
            importing_ids.setdefault("MeSH_ID", []).append(match)

    # For the rest of the databases, we simply search for exact matches in our list and the texts:
    if any((match:= hmdb) in text for hmdb in hmdb_ids):
        import_based_on.append("HMDB_ID")
        importing_ids.setdefault("HMDB_ID", []).append(match)
    if any((match:= name) in text for name in names):
        import_based_on.append("Name")
        importing_ids.setdefault("Name", []).append(match)

    # We return a list of a dict from keys to remove duplicates from the "reasons to import" list
    return list(dict.fromkeys(import_based_on)), importing_ids

def build_from_file(filepath, Neo4JImportPath, driver):
    """
    Imports a given metabolite from a sigle-metabolite containing file by checking its type
    and calling the appropriate import functions.

    Args:
        filepath (str): The path to the file in which will be imported
        Neo4JImportPath (str): The path which Neo4J will use to import data
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        This function does not provide a particular return, but rather imports the requested file

    .. NOTE:: The ``filepath`` may be absolute or relative, but it is transformed to a relative ``relpath``
        in order to remove possible influence of higher-name folders in the import type selection. This is
        also why the condition is stated as a big "if/elif/else" instead of a series of "ifs"
    """
    relpath = os.path.relpath(filepath, ".")
    filepath = os.path.abspath(filepath)
    fixedpath = os.path.basename(filepath).replace(" ", "_")
    shutil.copyfile(filepath, f"{Neo4JImportPath}/{fixedpath}")

    if "DrugBank" in relpath:
        DrugBankDataBase.build_from_file(f"{fixedpath}", driver)
        os.remove(f"{Neo4JImportPath}/{os.path.basename(fixedpath)}")

    elif "HMDB" in relpath:
        if "protein" in relpath:
            HumanMetabolomeDataBase.build_from_protein_file(f"{fixedpath}", driver)
        elif "metabolite" in relpath:
            HumanMetabolomeDataBase.build_from_metabolite_file(f"{fixedpath}", driver)
        os.remove(f"{Neo4JImportPath}/{os.path.basename(fixedpath)}")

    elif "SMPDB" in relpath:
        # NOTE: Since this adds a ton of low-resolution nodes, maybe have this db run first?
        # We will ignore the smpdb_pathways file because it doesnt have "real" identifiers
        if "proteins" in relpath:
            SmallMoleculePathWayDataBase.build_from_file(filepath, Neo4JImportPath, driver, "Protein")
        if "metabolites" in relpath:
            SmallMoleculePathWayDataBase.build_from_file(filepath, Neo4JImportPath, driver, "Metabolite")

    elif "ExposomeExplorer/components" in relpath:
            # NOTE: Since only "components" can result in a match based on our current criteria,
            #   we will build the DB starting with the components only.
            # Here, instead of using shutil.copyfile, we will use pandas to purge the _count columns when copying
            original_file = pd.read_csv(filepath)[[x for x in open(f"{filepath}").readline().rstrip().split(",")
                                                if not x.endswith('_count')]]
            original_file.to_csv(f"{Neo4JImportPath}/{os.path.basename(filepath)}", index=False)
            with driver.session() as session:
                    session.execute_write(ExposomeExplorerDataBase.add_components, os.path.basename(filepath))
            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")

            ExposomeExplorerDataBase.build_from_file( os.path.dirname(filepath),
                                                      Neo4JImportPath, driver, keep_counts_and_displayeds = False)

def import_based_on_all_files(all_files, Neo4JImportPath, driver, similarity, chebi_ids, names, hmdb_ids, inchis, mesh_ids):
    """
    A function that searches inside a series of lists, provided as arguments, and imports
    the metabolites matching those present in them iterating over a list of files which
    may contain relevant information to be imported

    Args:
        all_files(list): A list of all the posible files where we want to look for info
        Neo4JImportPath (str): The path which Neo4J will use to import data
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        similarity (bool): Whether to use similarity as a measure to import or not
        chebi_ids (list): A list of all the ChEBI_ID which are considered a reason to import
        names (list): A list of all the Name which are considered a reason to import
        hmdb_ids (list): A list of all the HMDB_ID which are considered a reason to import
        inchis (list): A list of all the InChI which are considered a reason to import
        mesh_ids (list): A list of all the MeSH_ID which are considered a reason to import
    """
    # Set up the progress bar
    with alive_bar( len(all_files), title="Scanning all files for matches...") as bar:
        i = 0 # Initialize counter for verbose messages on slurn
        # And search for them in the all_files list we created earlier on based on a series of criteria:
        for filepath in all_files:
            import_based_on, importing_ids = find_reasons_to_import_all_files(filepath,
                                                similarity, chebi_ids, names, hmdb_ids, inchis, mesh_ids)

            # Once we know the reasons to import (this is done so that it only cycles one
            # time through the code), we import the files themselves
            if import_based_on:
                build_from_file(filepath, Neo4JImportPath, driver)

                for item_type, items in importing_ids.items():
                    for item in items:
                        # We add the OriginalMetabolite label here; this will necessarily create duplicates,
                        # but this will be handled later on when the DB gets purged
                        with driver.session() as session:
                            session.execute_write(link_to_original_data, item_type, item, import_based_on)

            if i % 15000 == 0 and i > 1: logging.info(f"Scanned file: {i} / {len(all_files)}")
            i += 1; bar() # And advance, of course

def import_based_on_index(databasefolder, Neo4JImportPath, driver, similarity, chebi_ids, names, hmdb_ids, inchis, mesh_ids):
    """
    A function that searches inside a series of lists, provided as arguments, and imports
    the metabolites matching those present in them using a JSON file to map the bits of the databases where
    the relevant information lies

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found
            There *must* be an index.json file located in ``databasefolder``/index.json
        Neo4JImportPath (str): The path which Neo4J will use to import data
        driver (neo4j.Driver): Neo4J's Bolt Driver currently in use
        similarity (bool): Whether to use similarity as a measure to import or not
        chebi_ids (list): A list of all the ChEBI_ID which are considered a reason to import
        names (list): A list of all the Name which are considered a reason to import
        hmdb_ids (list): A list of all the HMDB_ID which are considered a reason to import
        inchis (list): A list of all the InChI which are considered a reason to import
        mesh_ids (list): A list of all the MeSH_ID which are considered a reason to import
    """
    # Set the databasefolder to be an absolute path
    databasefolder = os.path.abspath(databasefolder)

    # We assume there is an index at the ``databasefolder``
    index_path = f"{databasefolder}/index.json"

    # Define the function we will use for deciding on importing a given file
    def find_reasons_to_import_from_index(item_list, item_type, previous_files, driver):
        # We use ijson to prevent collapse when reading the huge index file
        index_file = open(index_path, "r")
        for record in ijson.items(index_file, item_type):
            for item in item_list:
                # We initialize the list of reasons to import this specific ID
                import_based_on = []
                # And get the files to import based on exact matches of the ID on index
                files_to_import = record.get(item, "")
                for each_file in list(record.get(item, "")):
                    if each_file not in previous_files:
                        previous_files.append(each_file)
                        import_based_on.append(item_type)

                # Only if exact inchi search fails, it makes sense to do similarity search
                if item_type == "InChI" and "InChI" not in import_based_on and similarity:
                    all_inchis_on_index = record.keys()

                    reason_to_import = find_reasons_to_import_inchi(all_inchis_on_index, item)
                    if reason_to_import:
                        import_based_on.append(list(set(reason_to_import.values())))
                        for each_file in list(record.get(item, "")):
                            if each_file not in previous_files:
                                previous_files.append(each_file)

                if import_based_on:
                    # We add the OriginalMetabolite label here; this will necessarily create duplicates,
                    # but this will be handled later on when the DB gets purged
                    with driver.session() as session:
                        session.execute_write(link_to_original_data, item_type, item, import_based_on)

        index_file.close()
        return previous_files

    all_files = [] # Calculate the files that we will need to import
    logging.info(f"Finding Reasons to Import...")
    all_files = find_reasons_to_import_from_index(chebi_ids, "ChEBI_ID", all_files, driver)
    all_files = find_reasons_to_import_from_index(hmdb_ids, "HMDB_ID", all_files, driver)
    all_files = find_reasons_to_import_from_index(mesh_ids, "MeSH_ID", all_files, driver)
    all_files = find_reasons_to_import_from_index(names, "Name", all_files, driver)
    all_files = find_reasons_to_import_from_index(inchis, "InChI", all_files, driver)

    i = 0 # And import them
    with alive_bar(len(all_files), title="Importing Selected Files...") as bar:
        for each_file in all_files:
            filepath = f"{databasefolder}/{each_file}"
            # We can only import a file if it exists, of course; files may not actually
            # exist if using an index generated with a different number of DBs / DB version
            if os.path.isfile(filepath):
                build_from_file(filepath, Neo4JImportPath, driver)
                if i % 150 == 0 and i > 1: logging.info(f"Importing file: {i} / {len(all_files)}")
            i += 1; bar() # And advance, of course

def link_to_original_data(tx, item_type, item, import_based_on):
    """
    Links a recently-imported metabolite to the original data (that which caused it to be imported) by creating an
    ``ÒriginalMetabolite`` node that is ``(n)-[r:ORIGINALLY_IDENTIFIED_AS]->(a)`` related to the imported data

    Args:
        tx (neo4j.work.simple.Session): The session under which the driver is running
        item_type (str): The property to match in the Neo4J DataBase
        item (dict): The value of property ```item_type```
        import_based_on (list): A list of the methods that turned out to be valid for import, such as Name, ChEBI_ID...

    Returns:
        neo4j.work.result.Result: A Neo4J connexion to the database that modifies it according to the CYPHER statement contained in the function.
    """
    return tx.run(f"""
                    MERGE (n {{ {item_type}:"{item}" }})
                        SET n:OriginalMetabolite
                        SET n.Reasons_To_Import =  "{",".join(import_based_on)}"
                    """)

def annotate_using_wikidata(driver):
    """
    Once we finish the search, we annotate the nodes added to the database using WikiData

    Args:
         driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.

    .. TODO:: When fixing queries, fix the main subscript also
    """
    with driver.session() as session, alive_bar( 53, title="Querying WikiData...") as bar:
        misc.repeat_transaction(WikiDataBase.add_wikidata_and_mesh_by_name(), driver); bar()
        # The ``query`` param is, remember, so as to remove the wikidata_id search which is by default
        misc.repeat_transaction(WikiDataBase.add_metabolite_info(query = "ChEBI_ID"), driver); bar()
        misc.repeat_transaction(WikiDataBase.add_drug_external_ids(query = "DrugBank_ID"), driver); bar()
        misc.repeat_transaction(WikiDataBase.add_more_drug_info(query = "DrugBank_ID"), driver); bar()

        misc.repeat_transaction(WikiDataBase.find_subclass_of_disease(), driver); bar()
        misc.repeat_transaction(WikiDataBase.find_subclass_of_disease(), driver); bar()
        misc.repeat_transaction(WikiDataBase.find_subclass_of_disease(), driver); bar()
        misc.repeat_transaction(WikiDataBase.find_instance_of_disease(), driver); bar()

        # For each of the 10 numbers a wikidata_id may have as ending
        for number in range(10):
            misc.repeat_transaction(WikiDataBase.add_disease_info(number=number), driver); bar()
            misc.repeat_transaction(WikiDataBase.add_drugs(number=number), driver); bar()
            misc.repeat_transaction(WikiDataBase.add_causes(number=number), driver); bar()
            misc.repeat_transaction(WikiDataBase.add_genes(number=number), driver); bar()

        misc.repeat_transaction(WikiDataBase.add_drug_external_ids(), driver); bar()
        misc.repeat_transaction(WikiDataBase.add_more_drug_info(), driver); bar()
        misc.repeat_transaction(WikiDataBase.add_yet_more_drug_info(), driver); bar()
        misc.repeat_transaction(WikiDataBase.add_gene_info(), driver); bar()
        misc.repeat_transaction(WikiDataBase.add_metabolite_info(), driver); bar()

def add_mesh_and_metanetx(driver):
    """
    Add MeSH Term IDs, Synonym relations and Protein interactions to existing nodes using MeSH and MetaNetX
    Also, adds Kegg Pathway IDs

    Args:
         driver (neo4j.Driver): Neo4J's Bolt Driver currently in use

    Returns:
        This function modifies the Neo4J Database as desired, but does not produce any particular return.
    """
    with driver.session() as session, alive_bar( 10, title="Querying MeSH and MetaNetX...") as bar:
    # We will also add MeSH terms to all nodes:
        misc.repeat_transaction(MeSHandMetaNetXDataBases.add_mesh_by_name(), driver); bar()
    # We also add synonyms:
        misc.repeat_transaction(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("Name"), driver); bar()
        misc.repeat_transaction(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("KEGG_ID"), driver); bar()
        misc.repeat_transaction(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("ChEBI_ID"), driver); bar()
        misc.repeat_transaction(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("HMDB_ID"), driver); bar()
        misc.repeat_transaction(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("InChI"), driver); bar()
        misc.repeat_transaction(MeSHandMetaNetXDataBases.write_synonyms_in_metanetx("InChIKey"), driver); bar()

    # And some protein interactions, together with their pathways, too:
        misc.repeat_transaction(MeSHandMetaNetXDataBases.find_protein_data_in_metanetx(), driver); bar()
        misc.repeat_transaction(MeSHandMetaNetXDataBases.find_protein_interactions_in_metanetx(), driver); bar()
        misc.repeat_transaction(MeSHandMetaNetXDataBases.get_kegg_pathways_for_metabolites(), driver); bar()

def main():
    """
    The function that executes the code

    .. NOTE:: This function disables rdkit's log messages, since rdkit seems to dislike the way
        some of the InChI strings it is getting from the databases are formatted

    .. TODO:: CAMBIAR NOMBRE A LOS MESH PARA INDICAR EL TIPO. AÑADIR NAME A LOS WIKIDATA
    .. TODO:: FIX THE REPEAT TRANSACTION FUNCTION
    .. TODO:: Match partial InChI based on DICE-MACCS
    .. TODO:: QUE FUNCIONE -> ACTUALMENTE ESTA SECCION RALENTIZA MAZO
    .. TODO:: CHECK APOC IS INSTALLED
    .. TODO:: FIX MAIN
    .. TODO:: MERGE BY INCHI, METANETX ID
    .. TODO:: Fix find_protein_interactions_in_metanetx
    .. TODO:: Mover esa funcion de setup a misc
    .. TODO:: EDIT conf.py

    .. TODO:: Document the following Schema Changes:
        * For Subject, we have a composite PK: Exposome_Explorer_ID, Age, Gender e Information
        * Now, more diseases will have a WikiData_ID and a related MeSH. This will help with networking. And, this diseases dont even need to be a part of a cancer!
        * The Gene nodes no longer exist in the full db? -> They do
    """

    # Parse the command line arguments
    # Done first in order to show errors if bad commands are issued
    parser = args_parser(); args = parser.parse_args()

    # If the session is set to be interactive, display the logging messages
    logging.basicConfig(format='%(message)s')
    if args.interactive: logging.getLogger().setLevel(logging.INFO)

    # We may also disable rdkit's log messages
    rdkit.RDLogger.DisableLog('rdApp.*')

    # And we read the query file
    raw_database = pd.read_csv(args.query, delimiter=',', header=0)

    # And connect to the Neo4J database
    driver = misc.connect_to_neo4j(args.adress, args.username, args.password)

    Neo4JImportPath = misc.get_import_path(driver)
    logging.info("Connected to Neo4J")

    # For each item in our "query" file, we will try to find matches:
    for index, row in raw_database.iterrows():

        # We start by cleaning the database (important if this is not the first run)
        with driver.session() as session:
            misc.repeat_transaction( misc.clean_database(), driver)
            logging.info("Cleaned DataBase")

        print(f"Searching Synonyms for Metabolite {index+1}/{len(raw_database)} ...")

        # Then, we fill the empty values so that the program does not crash
        row.fillna("", inplace=True)

        chebi_ids, names, hmdb_ids, inchis, mesh_ids = improve_search_terms(driver, row.ChEBI, row.Name,
                                                                            row.Identifier, row.InChI, row.MeSH)

        print(f"Annotating Metabolite {index+1}/{len(raw_database)} using Built-In DataBases...")

        if args.noindex:
            # We prepare a scan of all the files available on our "DataBases" folder
            # We will cycle through them to try and find matches
            all_files = misc.scan_folder(args.dbfolder)
            all_files = [ x for x in all_files if "index.json" not in x]

            import_based_on_all_files(all_files, Neo4JImportPath, driver, args.similarity,
                                      chebi_ids, names, hmdb_ids, inchis, mesh_ids)

        else:
            import_based_on_index(args.dbfolder, Neo4JImportPath, driver, args.similarity,
                                  chebi_ids, names, hmdb_ids, inchis, mesh_ids)

        print(f"Annotating Metabolite {index+1}/{len(raw_database)} using Web DataBases...")

        # We first purge the database by deleting useless noded that might overcharge the web queries
        misc.purge_database(driver, method = "delete")

        # Finally, we apply some functions that, although they could be run each time,
        # are so resource intensive that its better to just use once:

        if args.webdbs:
            # # We annotate the existing nodes using WikiData
            annotate_using_wikidata(driver)
            # Add their MeSH and MetaNetX IDs and synonyms
            add_mesh_and_metanetx(driver)

        # And purge any duplicates once again
        misc.purge_database(driver)

        # And save it in GraphML format
        with driver.session() as session:
            session.execute_write(misc.export_graphml,
                                  f"metabolite_{index+1}.graphml")

        if not os.path.exists(args.results): os.makedirs(os.path.abspath(args.results))

        shutil.copyfile(f"{Neo4JImportPath}/metabolite_{index+1}.graphml",
                        f"{os.path.abspath(args.results)}/metabolite_{index+1}.graphml")

        logging.info(f"Metabolite {index+1}/{len(raw_database)} processed. You can find "
                     f"a copy of the associated knowledge graph at "
                     f"{os.path.abspath(args.results)}/metabolite_{index+1}.graphml")

    logging.info("Metabolite Processing has finished")

if __name__ == '__main__':

    main()
