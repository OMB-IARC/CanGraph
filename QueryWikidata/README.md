<!--
SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>

SPDX-License-Identifier: MIT
-->

<div align="center"> <img src="_static/wikidata_schema.png" width="50%"> </div>
<br>

This script, created as part of my Master's Intenship at IARC, imports nodes from the [WikiData SPARQL Service](https://query.wikidata.org), creating a high-quality representation of the data therein. Although wikidata is manually curated using [the Wiki principles](https://en.wikipedia.org/wiki/Wiki), [some publications have found](https://pubmed.ncbi.nlm.nih.gov/32180547/) it might be a good source of information for life sciences, specially due to the breadth of information it contains. It also provides an export in GraphML format.

To run, it uses `alive_progress` to generate an interactive progress bar (that shows the script is still running through its most time-consuming parts) and the `neo4j` python driver. This requirements can be installed using: `pip install -r requirements.txt`.

To run the script itself, use:

`python3 build_database.py neo4jadress username databasepassword`

where:

* **neo4jadress**: is the URL of the database, in neo4j:// or bolt:// format
* **username**: the username for your neo4j instance. Remember, the default is neo4j
* **password**: the passowrd for your database. Since the arguments are passed by BaSH onto python3, you might need to escape special characters

Please note that there are two kinds of functions in the associated code: those that use python f-strings, which themselves contain text that *cannot* be directly copied into Neo4J (for instance, double brackets have to be turned into simple brackets) and normal multi-line strings, which can. This is because f-strings allow for variable customization, while normal strings dont.

An archived version of this repository that takes into account the gitignored files can be created using: `git archive HEAD -o ${PWD##*/}.zip`

Finally, please node that the general philosophy and approach of the queries have been taken from [Towards Data Science](https://towardsdatascience.com/lord-of-the-wiki-ring-importing-wikidata-into-neo4j-and-analyzing-family-trees-da27f64d675e), a genuinely useful web site.

## Important Notices

* Please ensure you have internet access, which will be used to connect to Wikidata's SPAQL endpoint and gather the necessary info.

* As Neo4J can run out of "Java Heap Space" if the number of nodes/properties to add is too high, the script has been divided in order to minimize said number: for instance, only nodes with a ```wikidata_id``` ending in a given number from 0 to 9 are processed at a time. This does not decrease performance, since these nodes would have been processed nontheless, but makes the script more reliable.

* What does impact performance, however, is having different functions for adding cancers, drugs, metabolites, etc, instead of having just one match for each created cancer node. This makes WikiData have to process more queries that are less heavy, which makes it less likely to time-out, but causes the script to run more slowly.

* The Neo4J server presents a somewhat unstable connection that is sometimes difficult to keep alive, as it tends to be killed by the system when you so much as look at it wrong. To prevent this from happening, you are encouraged to assign a high-priority to the server's process by using the ```nice``` or ```renice``` commands in Linux (note that the process will be called "Java", not "Neo4J")

* Another measure taken to prevent Neo4J's unreliability from stopping the script is the ```misc.repeat_transaction``` function, which insists a given number of times until either the problem is fixed or the error persists. This is because Neo4J tends to: random disconnects, run out of java heap space, explode... and WikiData tends to give server errors, have downtimes during the **14+ hours** the script takes to run, etc.

* The data present in the "graph.graphml" file comes from WikiData, and was provided by this service free of charge and of royalties under the permissive CC-0 license.
