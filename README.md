<!--
SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>

SPDX-License-Identifier: MIT
-->

<div align="center"> <img src="_static/main_schema.png" width="75%"> </div>
<br>

This Git Project, created as part of my Master's Intenship at IARC, contains a series of scripts that pulls information from a series of **five** databases from their native format (XML, CSV, etc) into a common, GraphML format, using a shared schema that has been defined to minimize the number of repeated nodes and properties. This databases are:

* Exposome-Explorer: A hand-curated, high-quality database of associations between metabolites, food intakes and outakes and different diseases, specially cancers.
* Human Metabolome DataBase: An detailed, electronic database containing detailed information about small molecule metabolites found in the human body.
* DrugBank: A unique bioinformatics and cheminformatics resource that combines detailed drug data with comprehensive drug target information.
* Small Molecule Pathway Database: An interactive database containing more than 618 small molecule pathways found in humans, More than 70% of which are unique to this DB
* WikiData: The world's largest collaboratively generated collection of Open Data worldwide.

Each of then have their unique advantages and disadvantages (size, quality, etc) but they have been chosen to work together and help in identifying metabolites and their potential cancer associations at IARC.

With regards to the schema, it can be consulted in detail in the ```new-schema.graphml``` file, which can itself be opened in Neo4J by calling: ```CALL apoc.import.graphml("new-schema.graphml", {useTypes:true, storeNodeIds:false, readLabels:True})``` after placing it in your Neo4J's import directory (you can find it in the settings shown after starting the server with ```sudo neo4j start```). It consists of a simplification of all the nodes present on the ```old-schema.graphml``` file (which itself represents the five different schemas that our five databases natively presented), arrived at by merging nodes and changing relationship names so that they are unique (and, thus, more actionable). One property, ```LabelName``` has been added as a dummy name to generate the image you can see in the header.

This repo contains two kind of scripts: first, some ```build_database.py``` scripts, which contain the information to re-build the databases in the common format from scratch, and are located in subsequent subfolders named after the database they come from (more info can be consulted on them on their respective READMEs) and a common ```main.py``` script, which can be used to query for sub-networks based solely on info presented on a **sample_input.csv** database of identified compounds which we would like to annotate.

## Intallation

To use this script, you should first clone it into your personal computer. The easiest way to do this is to [git clone](https://docs.codeberg.org/git/clone-commit-via-cli/) the repo:

1. Install git (if not already installed). On linux: ```sudo apt install git```
2. Clone the repo: ```git clone https://codeberg.org/FlyingFlamingo/graphify-databases```
3. Step into the directory ```cd graphify-databases```

Once the project has been installed, you **must** run ```setup.py```, a preparation script that guides you through the process of installing all five databases on your computer, so that then we can correctly process them and generate the sub-networks. You should also install the required python modules and run the setup script:

4. PIP install all dependencies: ```pip install -r requirements.txt```
5. Run the setup script: ```python3 setup.py```

Once this has been done, you are ready to start using the main script!

NOTE: If you do not wish to use git, you can manually download the repo by clicking [here](https://codeberg.org/FlyingFlamingo/CanGraph/archive/main.zip)

## Usage

To generate this sub-networks (the original idea of the project) you should run:

`python3 main.py neo4jadress databaseusername databasepassword databasefolder inputfile`

where:

* **neo4jadress**: is the URL of the database, in neo4j:// or bolt:// format
* **databasename**: the name of the database in use. If using the free version, there will only be one database per project (neo4j being the default name); if using the pro version, you can specify an alternate name here
* **databasepassword**: the passowrd for the **databasename** DataBase. Since the arguments are passed by BaSH onto python3, you might need to escape special characters
* **databasefolder**: The folder indicated to ```setup.py``` as the one where your databases will be stored
* **inputfile**: The location of the CSV file in which the program will search for metabolites. This file should be a Comma-Separated file, with the following format: ```MonoisotopicMass, SMILES, InChIKey, Name, InChI, Identifier, ChEBI```

All images in this repository are [CC-BY-SA-4.0 International](https://creativecommons.org/licenses/by-sa/4.0/) Licensed.
