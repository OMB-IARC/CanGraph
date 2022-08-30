---
sd_hide_title: true
---

<!--
SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>

SPDX-License-Identifier: MIT

CanGraph documentation master file, created by sphinx-quickstart on
Sat Aug 20 00:26:02 2022. You can adapt this file completely to your
liking, but it should at least contain the root `toctree` directive.
-->

# Main Page

```{toctree}
:hidden: true
:maxdepth: 3

ToDoList

CanGraph
```

::::{grid}
:reverse:
:gutter: 2 3 3 3
:margin: 4 4 1 2

:::{grid-item}
:columns: 12 4 4 4

```{image} ./_static/cangraph-mini.png
:width: 200px
:class: sd-m-auto sd-animate-grow50-rot20
```
:::

:::{grid-item}
:columns: 12 8 8 8
:child-align: justify
:class: sd-fs-3

A utility to study and analyse cancer-associated metabolites using knowledge graphs

```{button-ref} CanGraph
:ref-type: doc
:color: warning
:class: sd-rounded-pill float-left


Get started
```

<!-- The SVG rendering breaks latex builds for the GitHub badge, so only include in HTML -->
```{only} html
[![GitHub Stars](https://img.shields.io/github/stars/OMB-IARC/CanGraph?style=social)](https://github.com/OMB-IARC/CanGraph)
[![License](https://img.shields.io/github/license/OMB-IARC/CanGraph)](https://mit-license.org/)
```

:::
::::

CanGraph is a python program that allows you to extract information about a newly discovered or existing metabolite from the following databases:

<!-- Since emojis are not availaible in LaTeX, lets include them only in the Web build,
     replacing things for the LaTeX build-->

```{only} latex
* Human Metabolome DB
* DrugBank
* Exposome Explorer
* WikiData
* SMPDB
* MesH and MetaNetX
```

```{only} html
::::{grid} 1 1 2 3
:class-container: text-center
:gutter: 3

:::{grid-item-card}
:link: CanGraph.GraphifyHMDB
:link-type: doc
:class-header: bg-light

Human Metabolome DB 👨🏽
^^^^^^^^^^^^^^^^^^^^^

Detailed information on small molecule metabolites found in the human body
:::

:::{grid-item-card}
:link: CanGraph.GraphifyDrugBank
:link-type: doc
:class-header: bg-light

DrugBank 💉
^^^^^^^^^^

Detailed drug data with comprehensive drug target information

:::

:::{grid-item-card}
:link: CanGraph.ExposomeExplorer
:link-type: doc
:class-header: bg-light

Exposome Explorer 🔎
^^^^^^^^^^^^^^^^^^^

Dedicated to biomarkers of exposure to environmental risk factors for disease, specially cancers

:::

:::{grid-item-card}
:link: CanGraph.QueryWikidata
:link-type: doc
:class-header: bg-light

WikiData 📊
^^^^^^^^^^

Connect your book with Binder, JupyterHub, and other live environments
:::

:::{grid-item-card}
:link: CanGraph.GraphifySMPDB
:link-type: doc
:class-header: bg-light

SMPDB 🛤️
^^^^^^^^

Small Molecule Pathways, over 70% of which not found in any other db
:::

:::{grid-item-card}
:link: CanGraph.MeSHandMetaNetX
:link-type: doc
:class-header: bg-light

MeSH and MetaNetX 📍
^^^^^^^^^^^^^^^^^^^

Unified NameSpace for metabolites and Medical Subject Headings.
:::

::::
```

For this purpose, CanGraph accepts any of the following inputs:

* **InChI:** The International Chemical Identifier for the metabolite, which can be calculated using Open Source Software and identifies a metabolite with 99.99% accuracy based on its structure
* **InChIKey:** The Hashed, shortened version of the InChI, sometimes used for efficiency
* **Name:** A commonly accepted name for the metabolite; preferably, IUPAC's standarized name
* **HMDB_ID:** The Metabolite's Identifier in the Human Metabolome Database
* **ChEBI_ID:** The Metabolite's Identifier in the ChEBI Database

which must be provided as explained in [CanGraph's README](CanGraph.rst)

# Acknowledgements

::::{grid} 2 2 2 2

:::{grid-item}
:columns: 7
CanGraph has been possible thanks to the [International Agency for Research on Cancer's](https://www.iarc.who.int/) [OncoMetaBolomics Team](https://github.com/OMB-IARC), which provided funding for the project
:::

:::{grid-item}
:columns: 4

```{image} ./_static/IARC-logo.png
:class: m-auto
:width: 100px
```
:::

::::
