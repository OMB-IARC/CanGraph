<!--
SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>

SPDX-License-Identifier: MIT
-->

<div align="center"> <img src="_static/hmdb_schema.png" width="50%"> </div>
<br>

This script, created as part of my Master's Intenship at IARC, imports nodes from the [Human Metabolome Database](https://hmdb.ca/) (a high quality, database containing a list of metabolites and proteins associated to different diseases) to Neo4J format in an automated way, providing an export in GraphML format.

To run, it uses `alive_progress` to generate an interactive progress bar (that shows the script is still running through its most time-consuming parts) and the `neo4j` python driver. This requirements can be installed using: `pip install -r requirements.txt`.

To run the script itself, use:

`python3 main.py neo4jadress username databasepassword`

where:

* **neo4jadress**: is the URL of the database, in neo4j:// or bolt:// format
* **username**: the username for your neo4j instance. Remember, the default is neo4j
* **password**: the passowrd for your database. Since the arguments are passed by BaSH onto python3, you might need to escape special characters

Please note that there are two kinds of functions in the associated code: those that use python f-strings, which themselves contain text that *cannot* be directly copied into Neo4J (for instance, double brackets have to be turned into simple brackets) and normal multi-line strings, which can. This is because f-strings allow for variable customization, while normal strings dont.

An archived version of this repository that takes into account the gitignored files can be created using: `git archive HEAD -o ${PWD##*/}.zip`

## Important Notices

* Please ensure you have internet access, enough espace in your hard drive (around 5 GB) and read-write access in ```./xmlfolder```. The files needed to build the database will be stored there.

* There are two kinds of high-level nodes stored in this database: "Metabolites", which are individual compounds present in the Human Metabolome; and "Proteins", which are normally enzimes and are related to one or multiple metabolites. There are different types of metabolites, but they were all imported in the same way; their origin can be differenced by the "<biospecimen>" field on the corresponding "Concentration" nodes. You could run a query such as: ```MATCH (n:Metabolite)-[r:MEASURED_AT]-(c:Concentration) RETURN DISTINCT c.Biospecimen```

* Some XML tags have been intentionally not processed; for example, the <substituents> tag seemed like too much info unrelated to our project, or the <spectra> tags, which could be useful but seemed to only link to external DBs
