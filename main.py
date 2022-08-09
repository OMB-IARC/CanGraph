#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

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

# Import internal modules for the program
import miscelaneous as misc
from GraphifyDrugBank import build_database as DrugBankDataBase
from GraphifyHMDB import build_database as HumanMetabolomeDataBase
from GraphifySMPDB import build_database as SmallMoleculePathWayDataBase
from ExposomeExplorer import build_database as ExposomeExplorerDataBase
from QueryWikidata import build_database as WikiDataBase
from MeSHandMetaNetX import build_database as MeSHandMetaNetXDataBases


# First, we prepare a scan of all the files available on our "DataBases" folder
# We will cycle through them later on to try and find matches
all_files = []
for root,dirs,files in os.walk(sys.argv[4]):
    for filename in files:
        all_files.append( os.path.abspath(os.path.join(root, filename)) )

raw_database = pd.read_csv(sys.argv[5], delimiter=',', header=0)

# Set up the progress bar
with alive_bar(len(all_files)*len(raw_database)) as bar:

    # And connect to the Neo4J database
    instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
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

        # Then, we declare some lists to store the synonyms
        chebi_ids = []; names = []; hmdb_ids = []; inchis = []
        # And we store the original IDs there
        chebi_ids.append(row["ChEBI"].replace("CHEBI:", "")); names.append(row["Name"])
        hmdb_ids.append(row["Identifier"]); inchis.append(row["InChI"])

        # Now, we can calculate synonyms for each metabolite we will search
        misc.find_synonyms(driver, hmdb_ids, chebi_ids, inchis, names, "Name", row.Name)
        misc.find_synonyms(driver, hmdb_ids, chebi_ids, inchis, names, "InChI", row.InChI)
        misc.find_synonyms(driver, hmdb_ids, chebi_ids, inchis, names, "HMDB", row.Identifier)
        misc.find_synonyms(driver, hmdb_ids, chebi_ids, inchis, names, "ChEBI", row.ChEBI.replace("CHEBI:", ""))

        # And search for them in the all_files list we created earlier on based on a series of criteria:
        for filepath in all_files:
            with open(f'{filepath}', "r") as f:
                text = f.read(); relpath = os.path.relpath(filepath, ".")
                import_files = False; import_based_on = []

                # We try to find exact InChI matches:
                if any(inchi in text for inchi in inchis):
                    import_files = True; import_based_on.append("Exact InChI")
                # If none if found, we use the "Similarity Evaluator" metric
                elif "InChI=" in text:
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
                            Subject = rdkit.Chem.MolFromInchi(row["InChI"])
                            MACCSSubject = MACCSkeys.GenMACCSKeys(Subject)


                            DICE_MACCS = rdkit.DataStructs.DiceSimilarity(MACCSQuery, MACCSSubject)
                            if DICE_MACCS > 0.95:
                                import_files = True; import_based_on.append(f"DICE-MACCS {100*round(DICE_MACCS, 4)} % similarity")

                # For CHEBI, if we are using E-E, and since they dont have a prefix (i.e. they are only a number) we have to process the files.
                if "ExposomeExplorer/components" in relpath:
                    component = pd.read_csv(os.path.abspath(filepath))
                    chebi_query = component["chebi_id"]
                    # NOTE: Here, we remove the optional CHEBI: prefix
                    if chebi_query in [chebi_id.replace("CHEBI:", "").replace("chebi:", "") for chebi_id in chebi_ids]:
                        import_files = True; import_based_on.append("ChEBI")
                # And, even if its not E-E, we still need to convert it to tag format
                elif row["ChEBI"].replace("CHEBI:", "<chebi_id>") in text:
                    import_files = True; import_based_on.append("ChEBI")

                # For the rest of the databases, we simply search for exact matches in our list and the texts:
                if any(hmdb in text for hmdb in hmdb_ids):
                    import_files = True; import_based_on.append("HMDB ID")
                if any(name in text for name in names):
                    import_files = True; import_based_on.append("Name")

                # Once we know the reasons to import (this is done so that it only cycles one
                # time through the code), we import the files themselves
                if import_files == True:
                    if "DrugBank" in relpath:
                        shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                        DrugBankDataBase.build_from_file(f"{os.path.basename(filepath)}", driver)
                        os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")

                    if "HMDB" in relpath:
                        if "protein" in relpath:
                            shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                            HumanMetabolomeDataBase.build_from_protein_file(f"{os.path.basename(filepath)}", driver)
                            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                        elif "metabolite" in relpath:
                            shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                            HumanMetabolomeDataBase.build_from_metabolite_file(f"{os.path.basename(filepath)}", driver)
                            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")

                    if "SMPDB" in relpath:
                        # NOTE: Since this adds a ton of low-resolution nodes, maybe have this db run first?
                        # We will ignore the smpdb_pathways file because it doesnt have "real" identifiers
                        if "proteins" in relpath:
                            SmallMoleculePathWayDataBase.build_from_file(sys.argv[4], filepath, Neo4JImportPath, driver, "Protein")
                        if "metabolites" in relpath:
                            SmallMoleculePathWayDataBase.build_from_file(sys.argv[4], filepath, Neo4JImportPath, driver, "Metabolite")

                    if "ExposomeExplorer/components" in relpath:
                            # NOTE: Since only "components" can result in a match based on our current criteria, we will build the DB based on the components only.
                            # Here, instead of using shutil.copyfile, we will use pandas to purge the _count columns when copying
                            original_file = pd.read_csv(filepath)[[x for x in open(f"{filepath}").readline().rstrip().split(",")
                                                                if not x.endswith('_count')]]
                            original_file.to_csv(f"{Neo4JImportPath}/{os.path.basename(filepath)}", index=False)
                            with driver.session() as session:
                                    session.write_transaction(ExposomeExplorerDataBase.add_components, os.path.basename(filepath))
                            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")

                            ExposomeExplorerDataBase.build_from_file( os.path.dirname(filepath), Neo4JImportPath, driver, False)

                    #Each time we import some nodes, link them to our original data
                    #NOTE: This WILL duplicate nodes and relations, but we will fix this later when we purge the DB
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

                # TODO: CAMBIAR NOMBRE A LOS MESH PARA INDICAR EL TIPO. AÑADIR NAME A LOS WIKIDATA
                # TODO: FIX THE REPEAT TRANSACTION FUNCTION
                # TODO: Match partial InChIs based on DICE-MACCS
                # TODO: QUE FUNCIONE -> ACTUALMENTE ESTA SECCION RALENTIZA MAZO
                # TODO: CHECK APOC IS INSTALLED
                # TODO: FIX MAIN
                # TODO: MERGE BY INCHI, METANETX ID
                # TODO: Fix find_protein_interactions_in_metanetx

                # Section Schema Changes
                # NOTE: For Subject, we have a composite PK: Exposome_Explorer_ID, Age, Gender e Information
                # NOTE: Now, more diseases will have a WikiData_ID and a related MeSH. This will help with networking. And, this diseases dont even need to be a part of a cancer!
                # NOTE: The Gene nodes no longer exist in the full db

                # And advance, of course
                bar()

        #Once we finish the search, we annotate the existing nodes using WikiData
        #TODO: WHEN FIXING QUERIES, FIX THE MAIN SUBSCRIPT ALSO
        with driver.session() as session:
            misc.repeat_transaction(WikiDataBase.add_wikidata_to_mesh, 10, session, bar)
            misc.repeat_transaction(WikiDataBase.add_metabolite_info, 10, session, bar, number = None, query = "ChEBI_ID")
            misc.repeat_transaction(WikiDataBase.add_drug_external_ids, 10, session, bar, number = None, query = "DrugBank_ID")
            misc.repeat_transaction(WikiDataBase.add_more_drug_info, 10, session, bar, query = "DrugBank_ID")

            misc.repeat_transaction(WikiDataBase.find_subclass_of_cancer, 10, session, bar)
            misc.repeat_transaction(WikiDataBase.find_subclass_of_cancer, 10, session, bar)
            misc.repeat_transaction(WikiDataBase.find_subclass_of_cancer, 10, session, bar)
            misc.repeat_transaction(WikiDataBase.find_instance_of_cancer, 10, session, bar)
            for number in range(3):
                misc.repeat_transaction(WikiDataBase.add_cancer_info, 10, session, bar, number)
                misc.repeat_transaction(WikiDataBase.add_drugs, 10, session, bar, number)
                misc.repeat_transaction(WikiDataBase.add_causes, 10, session, bar, number)
                misc.repeat_transaction(WikiDataBase.add_genes, 10, session, bar, number)
            misc.repeat_transaction(WikiDataBase.add_drug_external_ids, 10, session, bar)
            misc.repeat_transaction(WikiDataBase.add_more_drug_info, 10, session, bar)
            misc.repeat_transaction(WikiDataBase.add_gene_info, 10, session, bar)
            misc.repeat_transaction(WikiDataBase.add_metabolite_info, 10, session, bar)

        # We will also add MeSH terms to all nodes:
        with driver.session() as session:
            misc.repeat_transaction(MeSHandMetaNetXDataBases.add_mesh_by_name(), 10, session, bar)
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

        # Finally, we purge the database by removing nodes considered as duplicated
        # We only purge once at the end, to limit processing time
        misc.purge_database(driver)

        # And save it in GraphML format
        with driver.session() as session:
            session.write_transaction(misc.export_graphml, f"metabolite_{index+1}.graphml")

        print(f"Metabolite {index+1}/{len(raw_database)} processed. You can find a copy of the associated knowledge graph at {os.path.abspath(sys.argv[4])}/metabolite_{index+1}.graphml")
        shutil.copyfile(f"{Neo4JImportPath}/metabolite_{index+1}.graphml", f"{os.path.abspath(sys.argv[4])}/metabolite_{index+1}.graphml")
        exit(0)
