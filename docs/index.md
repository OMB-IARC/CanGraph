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

Downloads

Tutorials

Presentations <https://www.pablomarcos.me/es/posts/master-en-biolog%c3%ada-computacional/tfm/>

CanGraph

ToDoList
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

```{button-ref} Downloads
:ref-type: doc
:color: warning
:class: sd-rounded-pill float-left


Get started
```

<!-- The SVG rendering breaks latex builds for the GitHub badge, so only include in HTML -->
```{only} html
[![GitHub Stars](https://img.shields.io/github/stars/OMB-IARC/CanGraph?style=social)](https://github.com/OMB-IARC/CanGraph)
[![Lines of Code](https://img.shields.io/tokei/lines/codeberg.org/FlyingFlamingo/CanGraph?color=77b5fe&label=Lines%20of%20Code)](https://codeberg.org/FlyingFlamingo/CanGraph)
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

The free knowledge base with 99,190,630 data items that anyone can edit
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

::::{grid} 1 1 2 2

:::{grid-item}
:columns: 11
This research was funded by the [Agence Nationale de la Recherche](https://anr.fr/en/) Project number [ANR-19-CE45-0021](https://anr.fr/Project-ANR-19-CE45-0021) (New approaches to bridge the gap between genome-scale metabolic networks and untargeted metabolomics – **MetClassNet**) and the [Deutsche ForschungsGemeinschaft](https://www.dfg.de/) Project number [431572533](https://gepris.dfg.de/gepris/projekt/431572533?language=en) (**MetClassNet**: new approaches to bridge the gap between genome-scale metabolic networks and untargeted metabolomics )

:::

:::{grid-item}
:columns: 5

```{image} ./_static/ANR-logo.png
:align: center
:target: https://anr.fr/Project-ANR-19-CE45-0021
:width: 100px
```
:::

:::{grid-item}
:columns: 5

```{image} ./_static/DFG-logo.jpg
:align: center
:target: https://gepris.dfg.de/gepris/projekt/431572533?language=en
:width: 100px
```
:::

:::{grid-item}
:columns: 11
You can check MetClassNet's Official Website [here](http://www.metclassnet.org/)
:::

:::{grid-item}
:columns: 11
Most of the work for this project has been carried out at the [International Agency for Research on Cancer's](https://www.iarc.who.int/) [OncoMetaBolomics Team](https://github.com/OMB-IARC)
:::

:::{grid-item}
:columns: 11
:margin: 3
<!-- Leave some space for IARC's logo -->

```{image} ./_static/IARC-logo-long.jpg
:align: center
:target: https://www.iarc.who.int/
:class: m-auto
:width: 50%
```
:::


::::
