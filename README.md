<!--
SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>

SPDX-License-Identifier: GPL-3.0-or-later
-->

# transition-database

<div align="center"> <img src="header.png" width="50%"> </div>
<br>

This script, created as part of my Master's Intenship at IARC, transitions the [exposome-explorer database](http://exposome-explorer.iarc.fr/) (a high quality, hand-curated database containing associations of foods and chemical compounds with cancer) to Neo4J format in an automated way, providing an export in GraphML format.

To run, it uses `alive_progress` to generate an interactive progress bar (that shows the script is still running through its most time-consuming parts) and the `neo4j` python driver. This requirements can be installed using: `pip install -r requirements.txt`.

To run the script itself, use:

`python3 main.py neo4jadress databasename databasepassword csvfolder`

where:

* **neo4jadress**: is the URL of the database, in neo4j:// or bolt:// format
* **databasename**: the name of the database in use. If using the free version, there will only be one database per project (neo4j being the default name); if using the pro version, you can specify an alternate name here
* **databasepassword**: the passowrd for the **databasename** DataBase. Since the arguments are passed by BaSH onto python3, you might need to escape special characters
* **csvfolder**: The folder where the CSV files for the Exposome Explorer database are stored. This CSVs have to be manually exported from the (confidential) database itself, and are NOT equivalent to those find in [exposome-explorer download's page](http://exposome-explorer.iarc.fr/downloads)
