<!--
SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>

SPDX-License-Identifier: MIT
-->

# Installing CanGraph

CanGraph *a python utility to study and analyse cancer-associated metabolites using knowledge graphs*, presents itself in two versions: a command-line interface that enables you to run the program on any Linux-powered machine, and a pre-packaged [Apptainer image](https://apptainer.org/).

This section of the tutorial concerns only the installation of the CanGraph **command-line interface** in a new machine. If you would rather use {ref}`one of the pre-packaged versions of the software <Downloads>`, you can skip directly to the {ref}`usage section <Running the Software>`

---

## Installing Dependencies

CanGraph is "a python utility to study and analyse cancer-associated metabolites using knowledge graphs"; so, of course, ``python3`` needs to be installed in your system; preferably, ``python3.9``. A bunch of other requirements, detailed in the table below, should also be installed. Since the program is designed to be installed in a Linux system, and has been developed to run in Ubuntu, orders to install the programs using the APT package manager are provided.

| Package     | Description                                            | Order to Install             |
|-------------|--------------------------------------------------------|------------------------------|
| python3.9  | Programming language used                              | sudo apt install python3.9  |
| pip         | A package manager for the python package index         | sudo apt install python3-pip |
| default-jdk | The java programming language, which Neo4j needs       | sudo apt install default-jdk |
| default-jre | The java runtime environment, which Neo4j also needs   | sudo apt install default-jre |
| neo4j       | The neo4j DBMS, which is the backbone of the program   | sudo apt install neo4j       |
| git         | A version control system, to get and manage software   | sudo apt install git         |
| curl        | cURL, a C utility to download things from the internet | sudo apt install curl        |

All this packages can also be installed using a one-line command, which just abbreviates all the ones present above. This should be preceded by an update of the system, like so:

* ``sudo apt update``

* ``sudo apt install python3.9 python3-pip default-jdk default-jre neo4j git curl``

* ``apt -y autoremove; apt -y clean``

## Cloning the Git repo

Having all the required software installed, and using Git as our version control management system, we can now proceed to download the software from its Web Repository.

On any directory, open a terminal window and run:

* ``git clone https://github.com/OMB-IARC/CanGraph``

The program will be now downloaded inside the "CanGraph" folder under the directory where you run git from.

## Installing Python Requirements

CanGraph uses a lot of python packages (collections of related modules) to run in a predictable and standardized way. They are listed in the ``requirements.txt`` file which comes bundled with the CanGraph package; you can either install all of them manually or run the following commands from the CanGraph folder (``cd CanGraph``):

* First, update PIP: ``python3 -m pip --no-cache-dir install --upgrade pip``

* And then, install: ``python3 -m pip install -r requirements.txt``

## Installing Apptainer

````{eval-rst}
.. note:: You need to install Apptainer to use the :ref:`pre-packaged files <Downloads>`, but you can use the :ref:`CLI <Running the Software>` without it.
````

To make the results from the CanGraph program easily reproducible and shareable, an Apptainer image is provided by us in the {ref}`downloads page <Downloads>`. If you want to re-generate the image yourself, either to make use of the simpler Private Image which we do not normally provide ourselves, or just to replicate our results, you will need to install Apptainer first. There are [a lot of tutorials availaible online](http://apptainer.org/docs/user/main/quick_start.html), which may or may not work on your computer, and which are sometimes split between the (old) Singularity and the (new) Apptainer brands.

To simplify the process for you, we have designed a pre-built script that installs a compatible version of [Apptainer](https://apptainer.org/) (v1.1.2) and [Go](https://go.dev/) (v1.43.0, the required programming language for Apptainer to run under). You can find it {ref}`here <Install Apptainer>`

## Building the Apptainer file

Once you have the Apptainer program installed and properly registered in your system (you can check this by running: ``apptainer --version``), you can build the SIF image file (the one you can use to later on {ref}`run the software <Running the Software>`) by running ``apptainer build container_name definition_file``, where:

* ``container_name`` is the final name of the container image file; for example, ``CanGraph.yaml``

* ``definition_file`` is the location of a file that follows the [Apptainer Definition Format](http://apptainer.org/docs/user/main/definition_files.html). In our software, we have two versions of this file, and thus two possible image files: one Private, which is built with the DataBases pre-bundled so that is is easier to use, although a bit more heavy to transmit; and one Public, in which DataBases need to be manually set up at least at first use and is thus less heavy, but more resource intensive.´
