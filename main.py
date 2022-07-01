#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: GPL-3.0-or-later
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

all_files = []
for root,dirs,files in os.walk(sys.argv[4]):
    for filename in files:
        all_files.append( os.path.abspath(os.path.join(root, filename)) )

raw_database = pd.read_csv(sys.argv[5], delimiter=',', header=0)

with alive_bar(len(all_files)*len(raw_database)) as bar:

    instance = f"{sys.argv[1]}"; user = f"{sys.argv[2]}"; passwd = f"{sys.argv[3]}"
    driver = GraphDatabase.driver(instance, auth=(user, passwd))

    Neo4JImportPath = misc.get_import_path(driver)
    print("Connected to Neo4J")

    for index, row in raw_database.iterrows():
        print(f"Importing metabolites: {index+1}/{len(raw_database)}")
        with driver.session() as session:
            session.run( misc.clean_database() )
        print("Cleaned DataBase")
        for filepath in all_files:
            with open(f'{filepath}', "r") as f:
                text = f.read(); relpath = os.path.relpath(filepath, ".")
                import_files = False; import_based_on = []

                # Similarity evaluator
                if row["InChI"] in text:
                    import_files = True; import_based_on.append("Exact InChI")
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

                # TODO: FIX CHEBI IMPORTS! NOW, THEY DONT WORK
                if row["ChEBI"] in text:
                    import_files = True; import_based_on.append("ChEBI")
                if row["Identifier"] in text:
                    import_files = True; import_based_on.append("HMDB ID")
                if row["Name"] in text:
                    import_files = True; import_based_on.append("Name")

                # This is so that it only cycles one time through the code
                if import_files == True:
                    if "DrugBank" in relpath:
                        shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                        DrugBankDataBase.build_from_file(f"{os.path.basename(filepath)}", driver)
                        os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                        DrugBankDataBase.purge_database(driver)

                    if "HMDB" in relpath:
                        if "protein" in relpath:
                            shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                            HumanMetabolomeDataBase.build_from_protein_file(f"{os.path.basename(filepath)}", driver)
                            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                            HumanMetabolomeDataBase.purge_database(driver)
                        elif "metabolite" in relpath:
                            shutil.copyfile(filepath, f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                            HumanMetabolomeDataBase.build_from_metabolite_file(f"{os.path.basename(filepath)}", driver)
                            os.remove(f"{Neo4JImportPath}/{os.path.basename(filepath)}")
                            HumanMetabolomeDataBase.purge_database(driver)

                    if "SMPDB" in relpath:
                        # NOTE: Since this adds a ton of low-resolution nodes, maybe have this db run first?
                        # We will ignore the smpdb_pathways file because it doesnt have "real" identifiers
                        if "proteins" in relpath:
                            SmallMoleculePathWayDataBase.build_from_file(sys.argv[4], filepath, Neo4JImportPath, driver, "Protein")
                            SmallMoleculePathWayDataBase.purge_database(driver)
                        if "metabolites" in relpath:
                            SmallMoleculePathWayDataBase.build_from_file(sys.argv[4], filepath, Neo4JImportPath, driver, "Metabolite")
                            SmallMoleculePathWayDataBase.purge_database(driver)

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

                    #Finally, we annotate the existing nodes using WikiData
                    #TODO: Fix the WikiData Script so that this actually does something!
                    #with driver.session() as session:
                        #misc.repeat_transaction(WikiDataBase.add_metabolite_info, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.add_genes, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.add_gene_info, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.add_drugs, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.add_drug_external_ids, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.add_more_drug_info, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.add_causes, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.find_subclass_of_cancer, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.find_instance_of_cancer, 10, session, bar)
                        #misc.repeat_transaction(WikiDataBase.add_cancer_info, 10, session, bar)

                    #Each time we import some nodes, link them to our original data
                    #NOTE: This WILL duplicate nodes and relations, but we will repeat this later

                    # TODO: Match partial InChIs based on DICE-MACCS
                    with driver.session() as session:
                        session.run( f"""
                                    MATCH (a) WHERE a.InChI = "{row["InChI"]}"
                                    MATCH (b) WHERE b.InChIKey = "{row["InChIKey"]}"
                                    MATCH (c) WHERE c.Name = "{row["Name"]}"
                                    MATCH (d) WHERE d.SMILES = "{row["SMILES"]}"
                                    MATCH (e) WHERE e.InChI = "{row["InChI"]}"
                                    MATCH (f) WHERE f.HMDB_ID = "{row["Identifier"]}"
                                    MATCH (g) WHERE g.Monisotopic_Molecular_Weight = "{row["MonoisotopicMass"]}"
                                    CREATE (n:OriginalMetabolite)
                                    SET n.InChI = "{row["InChI"]}", n.InChIKey = "{row["InChIKey"]}", n.Name = "{row["Name"]}",
                                        n.SMILES = "{row["SMILES"]}", n.HMDB_ID = "{row["Identifier"]}", n.ChEBI = "{row["ChEBI"]}",
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

                # And advance, of course
                bar()

        #Once we finish the search for each metabolite on the original table, we purge the database by removing EQUALLY EXACT metabolites
        # TODO: FIX PUBLICATIONS BEING MERGED!!
        # TODO: CHANGE NAME TO: ASSOCIATED_CANCER_METABOLITE TO ASSOCIATED DISEASE
        # TODO: DE DONDE SALEN LOS DISEASES SIN ID???
        # TODO: REMOVE EXPOSOME EXPLORER CROSS REFERENCERS
        # TODO: Fix publication Primary Key
        # TODO: Fix subject Primary Key
        # TODO: MIT LICENSE
        with driver.session() as session:
            session.write_transaction(misc.remove_duplicate_nodes, "", "n.InChI as sth", "WHERE n.InChI IS NOT null AND n:Metabolite or n:Protein")
            session.write_transaction(misc.remove_duplicate_nodes, "", "n.InChIKey as sth", "WHERE n.InChIKey IS NOT null AND n:Metabolite or n:Protein")
            session.write_transaction(misc.remove_duplicate_nodes, "", """n.InChI as inchi, n.InChIKey as inchikey, n.Name as name, n.SMILES as smiles,
                                                                          n.Identifier as hmdb_id, n.ChEBI as chebi, n.Monisotopic_Molecular_Weight as mass""",
                                                                       "WHERE n:Metabolite or n:Protein")
            session.write_transaction(misc.remove_duplicate_relationships)

        # And save it in GraphML format
        with driver.session() as session:
            session.write_transaction(misc.export_graphml, f"metabolite_{index+1}.graphml")

        print(f"Metabolite {index+1}/{len(raw_database)} processed. You can find a copy of the associated knowledge graph at {os.path.abspath(sys.argv[4])}/metabolite_{index+1}.graphml")
        shutil.copyfile(f"{Neo4JImportPath}/metabolite_{index+1}.graphml", f"{os.path.abspath(sys.argv[4])}/metabolite_{index+1}.graphml")
