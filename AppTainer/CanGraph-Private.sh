#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

# *********************** #
# CanGraph Private
# *********************** #

# This is the runscript for the CanGraph-Private APPTainer container, of which you can
# know more by reading the ``CanGraph-Private.yaml`` file, or by accessing the online docs

# *********************** #

# Then, map the supplied arguments (which follow ``main.py`` format) for ``setup.py``
# This is specially necessary to prevent empty arguments from messing up the python calls
# To do this, we use the ``remap_args`` function
function remap_args {
  if [ -z "$1" ]; then   # If no arguments were supplied
      echo $2 # Return the second value of the function
    else
      echo $1
  fi
}

databasefolder=$(remap_args $2 "$CANGRAPH_HOME/DataBases")
restoftheargs="${@:3}"

# First, we test that the script has the correct arguments and in the correct format
if [ $# -lt 1 ] || [[ " $*" == *" -h"* ]] || [[ " $*" == *" --h"* ]] || [[ " $*" == *" -help"* ]] || [[ " $*" == *" --help"* ]]; then
  python3 $CANGRAPH_HOME/main.py -h
  exit 1
# Here, we already know there are more than 1 cmd args
elif [[ "$( python3 -u $CANGRAPH_HOME/main.py --check_arguments $1 $databasefolder $restoftheargs  2>&1 )" =~ .*"error".* ]]; then
  python3 $CANGRAPH_HOME/main.py --check_arguments $1 $databasefolder $restoftheargs
  exit 1
fi

# Of course, if all of [neo4j_home, neo4j_username, neo4j_password] are set,
# this means that we should just use these values to connect to the database; i.e.
# that a valid neo4j instance already exists and that the user just wants to connect to it.
if [ -z "$3" ]; then

  # If any of [neo4j_username, neo4j_password] is missing, we should use the default or the old ones
  neo4j_username=$(remap_args $4 "neo4j")
  if  test -f ".neo4jpassword";
    then neo4j_password=$(remap_args $5 "$(head -n 1 .neo4jpassword)");
    else neo4j_password=$(remap_args $5 "neo4j");
  fi
  neo4j_home=$(remap_args $3 "./neo4j")

  # Then, we run the Neo4J setup script to check for already installed versions or install it
  python3 $CANGRAPH_HOME/setup.py --neo4j $neo4j_home --neo4j_username $neo4j_username --neo4j_password $neo4j_password
  # Please NOTE how this will also restart any neo4j sessions in existance; thus, the bolt port will be the default
  neo4j_adress="bolt://localhost:7687"

  python3 -c "import sys; sys.path.append('$CANGRAPH_HOME'); import miscelaneous as misc; \
              misc.kill_neo4j('$neo4j_home')"

  $neo4j_home/bin/neo4j stop > /dev/null; sleep 3

  # Start the Neo4J database where info will be stored
  $neo4j_home/bin/neo4j restart > /dev/null

  # Now, wait for neo4j to start
  python3 -c "import sys; sys.path.append('$CANGRAPH_HOME'); import miscelaneous as misc; \
              misc.sleep_with_counter(8, message = 'Waiting for neo4j to start...')"

else
  neo4j_adress=$3
fi

# And run the main executable with all the argumets that were passed on to CanGraph
# Note that, if the last 3 args where missing, they will be auto-substituted here by python
python3 $CANGRAPH_HOME/main.py $1 $databasefolder $neo4j_adress $neo4j_username $neo4j_password
