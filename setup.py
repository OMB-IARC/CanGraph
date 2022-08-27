#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that prepares the local environment, to be able to run the :obj:`~CanGraph.main` function.
This can be either run in an interactive way, requiring user input; or in a automatic way, in order to pre-configure
things, for example, if you are using the singularity package

For more details on how to run this script, please consult the package's README
"""

# Import external modules necessary for the script
from neo4j import GraphDatabase      # The Neo4J python driver
import os, sys, shutil               # Vital modules to interact with the filesystem
from time import sleep               # Cute go slowly
from zipfile import ZipFile          # Work with ZIP files

import pip                           # Interact with the Python Package Index directly
import miscelaneous as misc          # A collection of useful functions


def check_file(filepath):
    """
    Checks for the presence of a file, and exits the program

    Args:
        filepath (str): The path of the file its existence is being checked

    Raises:
        ValueError: If the file does not exist

    Returns:
        bool: True if successful, False otherwise.
    """
    if not os.path.exists(filename):
        raise ValueError("Missing file: {filename}. Please add the file and run the script back")
    return True

def initial_message(interactive=False):
    """
    Prompts the user with an initial message if the session is set to be interactive.

    Args:
        interactive (bool): Whether the session is set to be interactive or not
    """
    if interactive == True:
        print("Hi! This setup script will guide you through the necessary steps to prepare the running of the main script")
        sleep(1)
        print("This is because most of the databases we are using are confidential, so we cannot bundle them in the script")
        sleep(1)
        print("Fortunately, you will only need to run the setup once: if everything is ready, you can just directly run main.py (please read the README for instructions)")
        sleep(2)

def install_packages(requirements_file = None, package_name = None):
    """
    Automates installing packages using PIP

    Args:
        requirements_file (str): The path to a "requirements.txt" file, containing one requirement per line
        package_name (str): A package to be installed

    Raises:
        ValueError: If neither a ``requirements_file`` nor a ``package_name`` is provided
    """
    print("Checking for and installing packages using PIP...")
    if requirements_file != None:
        pip.main(["install", "-r", requirements_file])
    elif package != None:
        pip.main(["install", package_name])
    else:
        raise ValueError("You must give at least one requirement to install")

def setup_folders(interactive=False):
    """
    Creates the "Databases" folder if it does not exist. If it does, it exits so as not to overwrite.

    Args:
        interactive (bool): Whether the session is set to be interactive or not

    Raises:
        ValueError: If the Databases folder already exists (so as not to overwrite)

    Returns:
        bool: True if successful, False otherwise.
    """
    if os.path.exists("./DataBases"):
        if interactive == True:
            print("I see you already have a \"DataBases\" folder in this directory. May we use it for the project? (We may ovewrite its contents!")
            print("Use existing folder? [Y/n]", end=""); response = input()
        else:
            response = "n"

        if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
            raise ValueError("Please fix the folder situation. Exiting..."); sleep(1); exit(1)
    else:
        os.mkdir("DataBases")
        print("Created DataBases folder")
    return True

def setup_exposome(interactive=False):
    """
    Sets up the files relative to the Exposome Explorer database in the "Databases" folder,
    splitting them for easier processing later on.

    Args:
        interactive (bool): Whether the session is set to be interactive or not

    .. TODO:: What if ``interactive=False``?
    """
    print("First, please put the CSVs pertaining the exposome-explorer database on the ./DataBases/ExposomeExplorer path")
    sleep(1)
    print("This database is confidential, and CANNOT be found online. Please ask IARC for the files in case you need them.")
    sleep(1)
    if not os.path.exists("./DataBases/ExposomeExplorer"): print("Let me create said folder for you..."); os.mkdir("./DataBases/ExposomeExplorer")
    print("Once you are ready, press [ENTER]", end=""); response = input()

    check_file("./DataBases/ExposomeExplorer/units.csv");                check_file("./DataBases/ExposomeExplorer/subjects.csv")
    check_file("./DataBases/ExposomeExplorer/specimens.csv");            check_file("./DataBases/ExposomeExplorer/samples.csv")
    check_file("./DataBases/ExposomeExplorer/reproducibilities.csv");    check_file("./DataBases/ExposomeExplorer/publications.csv")
    check_file("./DataBases/ExposomeExplorer/microbial_metabolite_identifications.csv")
    check_file("./DataBases/ExposomeExplorer/metabolomic_associations.csv"); check_file("./DataBases/ExposomeExplorer/cancer_associations.csv")
    check_file("./DataBases/ExposomeExplorer/measurements.csv");         check_file("./DataBases/ExposomeExplorer/experimental_methods.csv")
    check_file("./DataBases/ExposomeExplorer/correlations.csv");         check_file("./DataBases/ExposomeExplorer/components.csv")
    check_file("./DataBases/ExposomeExplorer/cohorts.csv");              check_file("./DataBases/ExposomeExplorer/cancers.csv")

    if interactive == True: print("All checks OK")

    print("Now, I will split the some files into a number of one-liner CSVs, simplifying import")
    print("Please, bear in mind that I will remove some of the the original files to avoid problems; is that ok? [y/N]:", end="")
    response = input()
    if response != "y" and response != "Y" and response != "Yes" and response != "yes" and response != "YES":
        print("Exiting..."); sleep(1); exit(1)

    for filename in os.listdir("./DataBases/ExposomeExplorer/"):
        #if filename not in ["cancer_associations.csv", "metabolomic_associations.csv", "correlations.csv"]:
        if "components" in filename:
            misc.split_csv(filename.split("/")[-1], "./DataBases/ExposomeExplorer/")

def setup_hmdb(interactive=False):
    """
    Sets up the files relative to the HMDB database in the "Databases" folder,
    splitting them for easier processing later on.

    Args:
        interactive (bool): Whether the session is set to be interactive or not
    """
    if interactive == True:
        print("""Now, lets go on with the Human Metabolome DataBase.
                 I will download all the needed folders and store them in ./DataBases/HMDB""")
    if not os.path.exists("./DataBases/HMDB"): os.mkdir("./DataBases/HMDB")
    if interactive == True: print("Depending on your internet connection, this may take some time...")

    hmdb_urls = ["https://hmdb.ca/system/downloads/current/hmdb_proteins.zip",
                "https://hmdb.ca/system/downloads/current/urine_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/serum_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/csf_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/saliva_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/feces_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/sweat_metabolites.zip"
                ]

    for url in hmdb_urls:
        filename = f"{url.split('/')[-1].split('.')[0]}.xml"
        print(f"Downloading and Unzipping: {filename}...")
        misc.download_and_unzip(url, "./DataBases/HMDB")
        if filename == "hmdb_proteins.xml":
            print("File downloaded, now splitting its contents...")
            misc.split_xml(os.path.abspath(f"./DataBases/HMDB/{filename}"), "protein", "hmdb")
            os.remove(os.path.abspath(f"./DataBases/HMDB/{filename}"))
            sleep(1) # Some waiting time to avoid overheating
        else:
            print("File downloaded, now splitting its contents...")
            misc.split_xml(os.path.abspath(f"./DataBases/HMDB/{filename}"), "metabolite", "hmdb")
            os.remove(os.path.abspath(f"./DataBases/HMDB/{filename}"))
            sleep(1)  # Some waiting time to avoid overheating

    if interactive == True: print("Everything OK")

def setup_smpdb(interactive = False):
    """
    Sets up the files relative to the SMPDB database in the "Databases" folder,
    splitting them for easier processing later on.

    Args:
        interactive (bool): Whether the session is set to be interactive or not
    """
    if interactive == True:
        print("Super! Now, lets download and unzip the files related to the Small Molecule Pathway DataBase!")
        print("This may also take some time")
    if not os.path.exists("./DataBases/SMPDB"): os.mkdir("./DataBases/SMPDB")

    smpdb_urls = ["http://smpdb.ca/downloads/smpdb_pathways.csv.zip",
                "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip",
                "http://smpdb.ca/downloads/smpdb_proteins.csv.zip",
                "http://smpdb.ca/downloads/sequences/smpdb_protein.fasta.zip",
                "http://smpdb.ca/downloads/sequences/smpdb_gene.fasta.zip",
                ]

    for url in smpdb_urls:
        split = url.split('/')[-1].split('.')[0]
        print(f"Downloading and Unzipping: {split}...")
        if split in ["smpdb_metabolites", "smpdb_proteins"]:
            misc.download_and_unzip(url, f"./DataBases/SMPDB/{split}")
        else:
            misc.download_and_unzip(url, "./DataBases/SMPDB/")

    check_file("./DataBases/SMPDB/smpdb_pathways.csv");
    check_file("./DataBases/SMPDB/smpdb_proteins/");            check_file("./DataBases/SMPDB/smpdb_metabolites/")
    check_file("./DataBases/SMPDB/smpdb_gene.fasta");           check_file("./DataBases/SMPDB/smpdb_protein.fasta")

    if interactive == True: print("All checks OK")

def setup_drugbank(interactive = False):
    """
    Sets up the files relative to the SMPDB database in the "Databases" folder,
    splitting them for easier processing later on.

    Args:
        interactive (bool): Whether the session is set to be interactive or not

    .. TODO:: What if ``interactive=False``?
    """
    print("We are almost there! Now, we have to prepare the DrugBank database. If you have a DrugBank account, I can help you automate the process!")
    if not os.path.exists("./DataBases/DrugBank"): os.mkdir("./DataBases/DrugBank")
    print("Do you have an account [Y/n]", end=""); response = input()
    if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
        print("Please, ask for a DrugBank account and run the script back. The process must be approved by the DrugBank team, and may take a week or more."); exit(1)

    print("Please introduce your DrugBank username:", end=""); username = input()
    print("Please introduce your DrugBank password:", end=""); password = input()

    print("Downloading the full database using curl. Since its 1.4 GB, this may take some time...")
    if os.path.exists("./DataBases/./DataBases/DrugBank/full_database.zip"): os.remove("./DataBases/DrugBank/full_database.zip")
    if os.path.exists("./DataBases/./DataBases/DrugBank/full_database.xml"): os.remove("./DataBases/DrugBank/full_database.xml")
    try:
        os.system(f"curl -Lf --progress-bar -o ./DataBases/DrugBank/full_database.zip -u {username}:{password} https://go.drugbank.com/releases/5-1-9/downloads/all-full-database")
    except:
        print("Something went wrong with CURL. Was your username/password OK?"); exit(1)

    print("Extracting files from downloaded zip..")

    zf = ZipFile("./DataBases/DrugBank/full_database.zip")
    zf.extractall(path = "./DataBases/DrugBank/")
    zf.close()
    print("File downloaded, now splitting its contents...")
    misc.split_xml(os.path.abspath(f"./DataBases/DrugBank/full database.xml"), "drug", "drugbank")
    os.remove(os.path.abspath("./DataBases/DrugBank/full_database.zip"))
    os.remove(os.path.abspath("./DataBases/DrugBank/full database.xml"))

    print("Everything OK")

def final_message(interactive = False):
    """
    Prompts the user with a final message.

    Args:
        interactive (bool): Whether the session is set to be interactive or not
    """
    if interactive == True:
        print("This is all! Since the WikiData and MetaNetX databases are auto-consulted using RDF, there is no need to set up anything!")
    print("Setup Ready; you can proceed to run the main script now :p")
    exit(0)

def main():
    """
    The function that executes the code
    """

    initial_message()
    install_packages(requirements_file = "requirements.txt")
    setup_folders()
    setup_exposome()
    setup_hmdb()
    setup_smpdb()
    setup_drugbank()
    final_message()


if __name__ == '__main__':

    main()
