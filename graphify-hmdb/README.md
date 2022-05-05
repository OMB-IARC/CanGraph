<!--
SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>

SPDX-License-Identifier: GPL-3.0-or-later
-->

# query-wikidata

This script, created as part of my Master's Intenship at IARC, imports nodes from the [WikiData SPARQL Service](https://query.wikidata.org), creating a high-quality representation of the data therein. Although wikidata is manually curated using [the Wiki principles](https://en.wikipedia.org/wiki/Wiki), [some publications have found](https://pubmed.ncbi.nlm.nih.gov/32180547/) it might be a good source of information for life sciences, specially due to the breadth of information it contains. It also provides an export in GraphML format.

To run, it uses `alive_progress` to generate an interactive progress bar (that shows the script is still running through its most time-consuming parts) and the `neo4j` python driver. This requirements can be installed using: `pip install -r requirements.txt`.

To run the script itself, use:

`python3 main.py neo4jadress username databasepassword`

where:

* **neo4jadress**: is the URL of the database, in neo4j:// or bolt:// format
* **username**: the username for your neo4j instance. Remember, the default is neo4j
* **password**: the passowrd for your database. Since the arguments are passed by BaSH onto python3, you might need to escape special characters

Please note that there are two kinds of functions in the associated code: those that use python f-strings, which themselves contain text that *cannot* be directly copied into Neo4J (for instance, double brackets have to be turned into simple brackets) and normal multi-line strings, which can. This is because f-strings allow for variable customization, while normal strings dont.

An archived version of this repository that takes into account the gitignored files can be created using: `git archive HEAD -o ${PWD##*/}.zip`

Finally, please node that the general philosophy and approach of the queries have been taken from [Towards Data Science](https://towardsdatascience.com/lord-of-the-wiki-ring-importing-wikidata-into-neo4j-and-analyzing-family-trees-da27f64d675e), a genuinely useful web site.


TEST con otros que no sean orina
ONTOLOGY No hace na
I found substituent to be too much useless info
propiedades añado las que quiero
espectra pasando que ni se que es ni parece muy util
AÑADIR PATHWAYS

WE HAVE
/metabolite_associations -> Just a relation
go_classifications
gene_properties
protein_properties
general_references
metabolite_references -> Duplicated, are just references

