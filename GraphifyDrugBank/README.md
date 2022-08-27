<!--
SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>

SPDX-License-Identifier: MIT
-->

<div align="center"> <img src="_static/drugbank_schema.png" width="50%"> </div>
<br>

This package, created as part of my Master's Intenship at IARC, imports nodes from the [DrugBank Database](https://www.drugbank.com/) (a high quality database containing a drugs and proteins with some characteristics) to Neo4J format in an automated way, providing an export in GraphML format.

To run, it uses `alive_progress` to generate an interactive progress bar (that shows the script is still running through its most time-consuming parts) and the `neo4j` python driver. This requirements can be installed using: `pip install -r requirements.txt`.

To run the script itself, use:

`python3 main.py neo4jadress username databasepassword filename`

where:

* **neo4jadress**: is the URL of the database, in neo4j:// or bolt:// format
* **username**: the username for your neo4j instance. Remember, the default is neo4j
* **password**: the passowrd for your database. Since the arguments are passed by BaSH onto python3, you might need to escape special characters
* **filename**: the database's file name. For practical purpouses, it **must** be present in ./xmlfolder

Please note that there are two kinds of functions in the associated code: those that use python f-strings, which themselves contain text that *cannot* be directly copied into Neo4J (for instance, double brackets have to be turned into simple brackets) and normal multi-line strings, which can. This is because f-strings allow for variable customization, while normal strings dont.

An archived version of this repository that takes into account the gitignored files can be created using: `git archive HEAD -o ${PWD##*/}.zip`

## Important Notices

* To access the Drugbank database, you have to previously request access at: https://go.drugbank.com/public_users/sign_up

* Some XML tags have been intentionally not processed; for example, the <international-brands> tag didn't seem to dd anything new, and the <price> and <patents> tags are likely not relevant for our project. The <ahfs-codes> tag has no info (check ```cat full\ database.xml | grep "<ahfs-codes>"```) and the <pdb-entries> tag seemed like too much work (same for <salts>).
nOT USE <price/>  <patents/>because not of use to us

* Some :Proteins might also be :Drugs and viceversa: this makes the schema difficult to understand, but provides more meaningful data. This is also why some relations (```RELATED_WITH```, for instance) dont have a relation: this way, we avoid duplicates.
