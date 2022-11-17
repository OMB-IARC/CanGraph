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
# To do this, we use the ``remap_args`` function, which checks an arg and returns its value
function remap_args {
  : '
  A function that searches for an argument inside a text string.
  If it exists, it returns its value; else, it returns a given default
  Args:
    $arg_query (str): The option we are searching for;
      for instance, --dbfolder
    $default (str): The default argument, in case none is provided
    $original_args (str): The str containing the $1 string and
      its value, inmediately following it and separated by an space

  Returns:
    $arg_value (str): The value for $arg_query in $original_args
  '
  arg_query=$1; default=$2; original_args="${@:3}"

  function return_default_if_no_args {
    if [ -z "$2" ]; then   # If there was not a non-default value
        echo $1 # Return the default
      else
        echo $2 # Else return the non-default, provided value
    fi
  }

  if [[ "$original_args" == *"$arg_query"* ]]; then
    split_args_by_query=(${original_args//$arg_query/ })
    split_args_by_space=(${split_args_by_query// / })
    arg_value=${split_args_by_space[0]}

    echo $(return_default_if_no_args $default $arg_value)
  fi
}

# First, we test that the script has the correct arguments and in the correct format
if [ $# -lt 1 ] || [[ " $*" == *" -h"* ]] || [[ " $*" == *" --h"* ]] || [[ " $*" == *" -help"* ]] || [[ " $*" == *" --help"* ]]; then
  python3 $CANGRAPH_HOME/main.py -h
  exit 1
# Here, we already know there are more than 1 cmd args
elif [[ "$( python3 -u $CANGRAPH_HOME/main.py --check_arguments $@  2>&1 )" =~ .*"error".* ]]; then
  python3 $CANGRAPH_HOME/main.py --check_arguments $@
  exit 1
fi

databasefolder=$(remap_args "--dbfolder" "$CANGRAPH_HOME/DataBases" $@)
# If any of [neo4j_username, neo4j_password] is missing, we should use the default or the old ones
neo4j_username=$(remap_args "--username" "neo4j" $@)
if  test -f ".neo4jpassword";
  then neo4j_password=$(remap_args "--password" "$(head -n 1 .neo4jpassword)" $@)
  else neo4j_password=$(remap_args "--password" "neo4j" $@)
fi
neo4j_home="neo4j"

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

# else
#   neo4j_adress=$3
# fi

# And run the main executable with all the argumets that were passed on to CanGraph
# Note that, if the last 3 args where missing, they will be auto-substituted here by python
python3 $CANGRAPH_HOME/main.py $@
