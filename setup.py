#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that prepares the local environment, to be able to run the :obj:`~CanGraph.main`
and :obj:`~CanGraph.deploy` functions. This can be either run in an interactive way, requiring
user input; or in a automatic way, in order to pre-configure things, for example, if you are using
the singularity package

CanGraph.setup Usage
---------------------

To use this module:

.. argparse::
   :module: CanGraph.setup
   :func: args_parser
   :prog: python3 setup.py
   :nodefault:

CanGraph.setup Functions
-------------------------

This module is comprised of:
"""

# Import external modules necessary for the script
import neo4j                         # Import the whole package for error management
from neo4j import GraphDatabase      # The Neo4J python driver
import os, sys, shutil               # Vital modules to interact with the filesystem
import re                            # Work with regular expressions
import time                          # Manage the time, and wait times, in python
from zipfile import ZipFile          # Work with ZIP files
import subprocess                    # Manage python sub-processes
import logging                       # Make ``verbose`` messages easier to show
import random, string                # For now, used to generate passwords
from git import Repo                 # Manage Git directly from Python
import psutil                        # Kill the burden of the neo4j process
import argparse                      # Arguments pàrser for Python
from texttable import Texttable      # Draw cute tables in python
import pandas as pd                  # Analysis of tabular data
import json                          # Read JSON files from Python
from alive_progress import alive_bar # A cute progress bar that shows the script is still running

import miscelaneous as misc          # A collection of useful functions

# ********* Miscelaneous functions ********* #

def args_parser():
    """
    Parses the command line arguments into a more usable form, providing help and more

    Returns:
        argparse.ArgumentParser:
            A dictionary of the different possible options for the program as keys, specifying their set value.
            If no command-line arguments are provided, the help message is shown and the program exits.

    .. NOTE:: Note that, in Google Docstrings, if you want a multi-line ``Returns`` comment,
        you have to start it in a different line :(
    .. NOTE:: The return **must** be of type :obj:`argparse.ArgumentParser` for the ``argparse``
        directive to work and auto-gen docs
    .. NOTE:: The ``--all```option has to be adressed outside of this function in order to not
        mess up the ``argparse`` directive in sphinx
    .. NOTE:: By using :obj:`argparse.const` instead of :obj:`argparse.default`, the check_file function will check ""
        (the current dir, always exists) if the arg is not provided, not breaking the function; if it is, it checks it.
    """
    parser = argparse.ArgumentParser(
        description = "A python module that prepares the local environment, to be able to run the CanGraph.main "
                      "and CanGraph.deploy functions."
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    parser.add_argument("-i", "--interactive", action="store_true",
                        help=("tells the script if it wants interaction from the user "
                              "and information shown to them; similar to --verbose"))
    parser.add_argument("-a", "--all", action="store_true",
                        help=("runs all the options below at once; equivalent to -dgnr; "
                              "it DOES NOT activate the interactive mode"))

    parser.add_argument("--databases", nargs='?', const="DataBases", type=misc.check_file,
                        help="set up the databases from which the program will pull its info using the provided folder")
    parser.add_argument("--git", nargs='?', const=".git", type=misc.check_file,
                        help="prepare the git environment for the deploy script using the provided git folder")
    parser.add_argument("--requirements", nargs='?', const="requirements.txt", type=misc.check_file,
                        help=("installs all the requirements needed for all the possible options from the given requirements file"))
    parser.add_argument("-n", "--neo4j", nargs='?', const="neo4j",
                        help="set up  the neo4j local environment, to run from the provided folder")

    parser.add_argument('--neo4j_username', nargs='?', default="neo4j",
                        help="the username for the neo4j database")
    parser.add_argument('--neo4j_password', nargs='?', default="neo4j",
                        help="the password for the neo4j database")

    # If no args are provided, show the help message
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser

# ********* Add some verbose messages ********* #

def initial_message():
    """
    Prompts the user with an initial message if the session is set to be interactive.

    Args:
        interactive (bool): Whether the session is set to be interactive or not
    """
    logging.info("Hi! This setup script will guide you through the necessary steps to prepare ")
    logging.info("for the running of the main script")

    time.sleep(1)

    logging.info("This is because most of the databases we are using are confidential,")
    logging.info("so we cannot bundle them in the module. Also, we didn't include neo4j")
    logging.info("in the main repo to keep it light and due to copyright concerns")

    time.sleep(1)

    logging.info("Fortunately, you will only need to run the setup once if you use the --all ")
    logging.info("option; and, if everything is ready, you can just directly run main.py ")
    logging.info("(please read its README for more instructions on usage)\n")
    time.sleep(2)

def final_message(interactive = False):
    """
    Prompts the user with a final message.

    Args:
        interactive (bool): Whether the session is set to be interactive or not
    """

    logging.info("\nIf you selected the --all option, you may now proceed to run the main script;")
    logging.info("in any other case, please make sure the parts relevant to the code you will ")
    logging.info("run are propperly configure (you may call ``python3 setup.py --help`` for more ")
    logging.info("info on available options)")

    exit(0)

# ********* Install the requirements ********* #

def install_packages(requirements_file = None, package_name = None, interactive = False):
    """
    Automates installing packages using PIP

    Args:
        requirements_file (str): The path to a "requirements.txt" file, containing one requirement per line
        package_name (str): A package to be installed
        interactive (bool): Whether the session is set to be interactive or not

    Raises:
        ValueError: If neither a ``requirements_file`` nor a ``package_name`` is provided
    """
    logging.info("Checking for and installing packages using PIP...")

    quiet = "--quiet" if interactive == False else ""

    if requirements_file != None:
        sp = subprocess.run(["pip", "install", "-r", requirements_file, f"{quiet}"], stdout=subprocess.PIPE)
    elif package != None:
        sp = subprocess.run(["pip", "install", package_name, f"{quiet}"], stdout=subprocess.PIPE)
    else:
        raise ValueError("You must give at least one requirement to install")

    sp.check_returncode() # Raise error if error is found

    print("Installed Requirements")

# ********* Set Up the Git environment for the deploy script ********* #

def setup_git(path_to_repo = ".git"):
    """
    Set Up the Git environment for the :obj:`~CanGraph.deploy` script. It does so by removing
    any existing remotes and setting two new ones: github and codeberg, with their respective branches

    Args:
        path_to_repo (str): The path to the Git repo; by default, ``.git``
    """
    # Remove any existing remotes so that they dont interfere
    logging.info("Now, I will set up the git environment so that the deploy functions may work")
    logging.info("First, let me remove any existing remotes")
    repo = Repo(path_to_repo)
    remotes = repo.git.remote().split("\n")
    for remote in remotes:
        if len(remote) > 0:
            repo.git.remote("rm", remote)

    # Replace them with the new, canonic ones
    logging.info("And replace them with the new ones")
    repo.git.remote("add", "github", "https://github.com/OMB-IARC/CanGraph")
    repo.git.remote("add", "codeberg", "https://codeberg.org/FlyingFlamingo/CanGraph")

    # Finally, set the remotes to track the correct branches
    logging.info("Finally, set each branch to track the correct remotes...")
    repo.git.fetch("--all")

    print("The git is now configured :p")

# ********* Manage the Neo4J Database location, status and credentials ********* #

def find_neo4j_installation_status(neo4j_home = "neo4j", neo4j_username = "neo4j", neo4j_password="neo4j"):
    """
    Finds the installation status of Neo4J by trying to use it normally, and analyzing any thrown exceptions

    Args:
        neo4j_home (str): the installation directory for the ``neo4j`` program; by default, ``neo4j``
        neo4j_username (str): the username for the neo4j database; by default ``neo4j``
        neo4j_password (str): the password for the neo4j database; by default ``neo4j``

    Returns:
        list:
            A list of two booleans: whether neo4j exists at ``neo4j_home``,
            and whether the supplied credentials are valid or not
    """
    # First, declare the return variables as False by default
    neo4j_exists = False; default_credentials = False

    # And set neo4j_home as an absolute path, in order to standardize
    neo4j_home = os.path.abspath(neo4j_home)

    # To find out neo4j's installation status, we will try to run it and use possible exceptions
    # as a hint as to whether there are problems or not.
    try:

        # Fist, we try to start neo4j using the Linux executable
        subprocess.run([f"{os.path.abspath(neo4j_home)}/bin/neo4j", "start"],
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        logging.info("I see you do have Neo4J already installed and working!")
        logging.info("Now, lets start it and check if it presents a known user/password pair")

        # If no exceptions are thrown, then it exists! Lets try to start it then
        neo4j_exists = True; misc.sleep_with_counter(8, message = "Waiting for neo4j to start...")

        driver = misc.connect_to_neo4j(port = "bolt://localhost:7687", username = neo4j_username, password = neo4j_password)

        # We will try to get the Import Path just as a test to see if the default auths are present
        Neo4JImportPath = misc.get_import_path(driver)

        # If there are no exceptions, obviously, then the credentials you supplied us with are working!
        logging.info("Installation with known credentials found! I'll work with that :p")
        default_credentials = True

    # If it throws an ``AuthError``, its because the credentials are not the default
    except neo4j.exceptions.AuthError as error:
        logging.info("Neo4J credentials are not known. We will need to re-install...")
        logging.info("Let me stop neo4j first...")
        default_credentials = False

    # If they the default (neo4j/neo4j), it throws a ``ClientError``, because it requires a password change first
    except neo4j.exceptions.ClientError as error:
        logging.info("Installation with known credentials found! I'll work with that :p")
        default_credentials = True

    # If the executable is not found in ``neo4j_home``. we will need to install it
    except FileNotFoundError as error:
        logging.info("No Neo4J installation has been found. Defaulting to installing from website...")
        neo4j_exists = False

    # In any other case, re-install it, too
    except Exception as error:
        logging.info("An unknown error has been raised. Defaulting to installing from website...")
        logging.info(f"Error was: {error}")
        neo4j_exists = False

    misc.kill_neo4j(neo4j_home)

    return neo4j_exists, default_credentials

def install_neo4j(neo4j_home = "neo4j", interactive = False, version = "4.4.0"):
    """
    Installs the neo4j database program in the ``neo4j_home`` folder, by getting it from the internet
    according to the Operating System the script is been run in (aims for multi-platform!)

    Args:
        neo4j_home (str): the installation directory for the ``neo4j`` program; by default, ``neo4j``
        interactive (str): tells the script if it wants interaction from the user and information shown to them
        version (str): the version of the neo4j software that we wish to install
    """

    neo4j_home = os.path.abspath(neo4j_home)
    neo4j_relt = os.path.relpath(neo4j_home)

    # And install the neo4j executable
    if os.path.exists(neo4j_home):
        if interactive == True:
            print(f"I need to use folder {neo4j_relt} to store the neo4j program, but path {neo4j_relt} already exists.")
            print("Allow me to remove it? [Y/n]", end=""); response = input()
            if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
                print("Please fix the folder situation. Exiting..."); exit(1)
            print("", end="")
        shutil.rmtree(neo4j_home)

    logging.info(f"Installing neo4j in {neo4j_home}...")
    if  sys.platform == "linux" or sys.platform == "linux1" or sys.platform == "linux2" or sys.platform == "darwin":
        misc.download_and_untargz(f"https://neo4j.com/artifact.php?name=neo4j-community-{version}-unix.tar.gz", "/tmp/")
    elif  sys.platform == "win32":
        misc.download_and_untargz(f"https://neo4j.com/artifact.php?name=neo4j-community-{version}-windows.zip", "/tmp/")
    shutil.move(f"/tmp/neo4j-community-{version}/", neo4j_home)
    logging.info("Installation done!")

def update_neo4j_confs(key, value, conf_file = "neo4j/conf/neo4j.conf"):
    """
    Updates a preference on neo4j's ``conf_file``, given its name (``key``) and its expected ``value``
    If a preference is set with a value other than ``key``, said value will be overwritten; if it is commented,
    it will be uncommented (thanks to regex!)

    Args:
        conf_file (str): The location for ne4j's configuration file; usually, it should be ``neo4j/conf/neo4j.conf``
        key (str): The key for neo4j's configuration parameter that is being set up
        value (str): The value for said parameter
    """
    old_conf = f"{os.path.splitext(conf_file)[0]}.old"
    shutil.move(conf_file, old_conf)
    conf_text = ""
    with open(old_conf, "r") as conf:
        conf_text = conf.read()
        search_term = key.replace(".", "\.");
        if key in conf_text:
            conf_text = re.sub(f"\#?{search_term}=.*", f"{key}={value}", conf_text)
        else:
            conf_text = conf_text + "\n" + f"{key}={value}"
    with open(conf_file, "w") as conf:
        conf.write(conf_text)

    os.remove(old_conf)

def configure_neo4j(neo4j_home = "neo4j"):
    """
    Modifies the Neo4J conf file according to some recommendations provided by ``memrec``, neo4j's
    memory recommendator. It also enables the Awesome Procedures On Cypher (APOC) plugin from
    Neo4j Labs, and enables other basic confs such as file export and import or bigger timeouts

    Args:
        neo4j_home (str): the installation directory for the ``neo4j`` program; by default, ``neo4j``

    .. NOTE:: In order to make the setup more consistent, this function also forces the Neo4JImportPath
        (``dbms.directories.import``) to be presented in an absolute way, instead of being relative
        to ``neo4j_home``
    """
    neo4j_home = os.path.abspath(neo4j_home)
    logging.info("Now, lets configure the neo4j environment")

    # Enable the APOC plugin
    logging.info("First, let me enable the APOC plugin, which we will use heavily...")

    apoc_regex = re.compile('apoc\-.*\-core\.jar')
    files_list = [x.name for x in os.scandir(f"{neo4j_home}/labs/")]

    for line in files_list:
        if apoc_regex.match(line):
            apoc_location = apoc_regex.match(line).group(0)

    if not os.path.exists(f"{neo4j_home}/plugins/{apoc_location}"):
        shutil.copyfile(f"{neo4j_home}/labs/{apoc_location}", f"{neo4j_home}/plugins/{apoc_location}")

    # And modify the ``conf`` file
    logging.info("And, now, I will modify the conf file to update the preferences ")
    logging.info("based on sugestions provided by neo4j-admin's memrec feature")

    sp = subprocess.run([f"{neo4j_home}/bin/neo4j-admin", "memrec"], stdout=subprocess.PIPE)
    sp.check_returncode() # Raise error if error is found

    for line in str(sp.stdout).split("\\n"):
        if line.startswith("dbms."):
            key = line.split("=")[0]; val = line.split("=")[1]
            update_neo4j_confs(key, val)

    update_neo4j_confs("apoc.import.file.enabled", True)
    update_neo4j_confs("apoc.export.file.enabled", True)
    update_neo4j_confs( "apoc.http.timeout.connect", 60000)
    update_neo4j_confs("apoc.http.timeout.read", 180000)

    # Dont do this! It will kill transactions that want to happen!
    # update_neo4j_confs("dbms.memory.transaction.global_max_size", "1024m")

    # Force the "import" path to be absolute:
    update_neo4j_confs("dbms.directories.import", f"{neo4j_home}/import")

def change_neo4j_password(new_password, old_password = "neo4j", user = "neo4j", database = "system", neo4j_home = "neo4j"):
    """
    Changes the neo4j password for user ``user``, from ``old_password`` to ``new_password``, by using a
    simple query in cypher-shell

    Args:
        neo4j_home (str): the installation directory for the ``neo4j`` program; by default, ``neo4j``
        new_password (str): the new password for the database
        old_password (str): the old password for the database, needed for identification.
        user (str): the user for which the password is being changed.
        database (str): the name of the database for which we want to modify the password.
            By default, it is ``system``, since Neo4J's community edition  only allows for one database
    """

    neo4j_home = os.path.abspath(neo4j_home) # Make the path absolute, to standardize

    # Restart the Neo4J database where we want to operate
    subprocess.run([f"{neo4j_home}/bin/neo4j", "restart"],
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    try:
        driver = misc.connect_to_neo4j(port = "bolt://localhost:7687", username = "neo4j", password = "neo4j")
        Neo4JImportPath = misc.get_import_path(driver)
    except:
        pass

    misc.sleep_with_counter(8, message = "Waiting for neo4j to start...")

    # Alter the password using cypher-shell
    sp = subprocess.run([f"{neo4j_home}/bin/cypher-shell", "-d", f"{database}", "-u", f"{user}",
                        "-p", f"{old_password}",
                        f"ALTER CURRENT USER SET PASSWORD FROM '{old_password}' TO '{new_password}'"],
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    sp.check_returncode() # Raise error if error is found

    misc.kill_neo4j(neo4j_home) # Kill the Neo4J process

def setup_neo4j(neo4j_home = "neo4j", neo4j_username = "neo4j", neo4j_password = "neo4j", interactive = False):
    """
    Sets ups the neo4j environment in ``neo4j_home``, so that the functions in :obj:`~CanGraph.main` can propperly
    function. Using the functions present in this module, it finds if neo4j is installed with default credentials,
    and, if not, it installs it, changing the default password to a new one, and returning its value

    Args:
        neo4j_home (str): the installation directory for the ``neo4j`` program; by default, ``neo4j``
        neo4j_username (str): the username for the neo4j database; by default ``neo4j``
        neo4j_password (str): the password for the neo4j database; by default ``neo4j``
        interactive (str): tells the script if it wants interaction from the user and information shown to them

    Returns:
        str: The password that was set up for the new neo4j database. This is also written to .neo4jpassword

    .. NOTE:: This has been designed to be used with a ``neo4j_home`` located in the WorkDir, but can be used
        in any other location with read/write access, or even with ``apt install`` installed versions! Just
        find its ``neo4j_home``, make sure it has r/w access, and provide it to the program!
    .. NOTE:: If no neo4j_password is provided or if neo4j_password = "neo4j", the function will check for a previously
        created ".neo4jpassword" file, signalling a possible pre-existing database with known credentials
    """

    # Fist, we need to kill any existing neo4j instances. As explained in :obj:`~CanGraph.miscelaneous.kill_neo4j`,
    #this is because the process may be lingering on otherwise, making the program crash
    if os.path.exists(neo4j_home): misc.kill_neo4j(neo4j_home)

    # Display an initial message
    print("Now, lets set up the Neo4J environment")
    logging.info("First, I will check for an existing Neo4J installation with the default")
    logging.info("user/passwd, as this might save time.")

    # Set the variables that will help us find the status of neo4j's installation to ``False`` by default
    neo4j_exists = False; default_credentials = False

    # Check for existing credentials
    if neo4j_password == "neo4j" and os.path.exists(".neo4jpassword"):
        with open('.neo4jpassword', "r") as password_file:
            old_password = password_file.readline().rstrip()
    else: old_password = neo4j_password

    # Find out the status of neo4j installation
    neo4j_exists, default_credentials = find_neo4j_installation_status(neo4j_home, neo4j_username, old_password)

    # If needed, install neo4j
    if neo4j_exists == False or default_credentials == False:
        install_neo4j(neo4j_home, interactive)

        # In which case, we will need to re-generate the password (or change it if its the default and thus errors)
        logging.info("Now, let me generate a new password for the neo4j database and establish it")
        characters = string.ascii_letters + string.digits
        new_password = ''.join(random.choice(characters) for i in range(12))

        change_neo4j_password(new_password, "neo4j", "neo4j", "system", "neo4j")
        if interactive == True: print(f"Neo4J's password has been set to: {new_password}")
        with open('.neo4jpassword', "w") as password_file:
            password_file.write(new_password)
        if interactive == True: print(f"This password has been written to ``.neo4jpassword``")

    # If the credentials provided are OK, then no need to change anything of course
    else: new_password = old_password

    # Finally, we configure the Neo4J environment. This might be duplicated, but checking
    # if the configs are already there will almost take the same time, so easier this way
    configure_neo4j(neo4j_home)

    # We Kill Neo4J one last time:
    misc.kill_neo4j(neo4j_home)

    print("The neo4j environment is ready for use") # And done!

    return new_password

# ********* Set Up the DataBases from which the ``main`` function will get its data ********* #

def setup_folders(databasefolder = "./DataBases", interactive=False):
    """
    Creates the ``databasefolder`` if it does not exist. If it does, it either asks
    before overwriting in ``interactive`` mode, or directly overwrites in auto mode.

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found
        interactive (bool): Whether the session is set to be interactive or not

    Raises:
        ValueError: If the Databases folder already exists (so as not to overwrite)

    Returns:
        bool: True if successful, False otherwise.
    """

    # Set the databasefolder to be an absolute path
    databasefolder = os.path.abspath(databasefolder)

    # If the path exists, remove it and create anew; else, pass
    if os.path.exists(databasefolder):
        if interactive == True:
            print(f"I see you already have a {databasefolder} folder in this directory. "
                   "May we use it for the project? (We may ovewrite its contents!)")
            print("Use existing folder? [Y/n]", end=""); response = input()
        else:
            response = "y"

        if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
            print("Please fix the folder situation. Exiting..."); time.sleep(1); exit(1)
    else:
        os.mkdir(databasefolder)
        logging.info("Created DataBases folder")

    return True

def check_exposome_files(databasefolder = "./DataBases"):
    """
    Checks for the presence of all the files that should be in "``databasefolder``/ExposomeExplorer"`
    for the ExposomeExplorer part of the script to run

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found

    Returns:
        bool:
            True if everything went okay; False otherwise. If False,
            Exposome-Explorer should not be used as a data source
    """

    # Set the databasefolder to be an absolute path
    databasefolder = os.path.abspath(databasefolder)

    # Check for the presence of all files; if :obj:`CanGraph.miscelaneous.check_file` throws an error, set exposome_explorer_ok to ``False``
    try:
        misc.check_file(f"{databasefolder}/ExposomeExplorer/units.csv");                misc.check_file(f"{databasefolder}/ExposomeExplorer/subjects.csv")
        misc.check_file(f"{databasefolder}/ExposomeExplorer/specimens.csv");            misc.check_file(f"{databasefolder}/ExposomeExplorer/samples.csv")
        misc.check_file(f"{databasefolder}/ExposomeExplorer/reproducibilities.csv");    misc.check_file(f"{databasefolder}/ExposomeExplorer/publications.csv")
        misc.check_file(f"{databasefolder}/ExposomeExplorer/microbial_metabolite_identifications.csv")
        misc.check_file(f"{databasefolder}/ExposomeExplorer/metabolomic_associations.csv"); misc.check_file(f"{databasefolder}/ExposomeExplorer/cancer_associations.csv")
        misc.check_file(f"{databasefolder}/ExposomeExplorer/measurements.csv");         misc.check_file(f"{databasefolder}/ExposomeExplorer/experimental_methods.csv")
        misc.check_file(f"{databasefolder}/ExposomeExplorer/correlations.csv");          misc.check_file(f"{databasefolder}/ExposomeExplorer/components.csv")
        misc.check_file(f"{databasefolder}/ExposomeExplorer/cohorts.csv");              misc.check_file(f"{databasefolder}/ExposomeExplorer/cancers.csv")
        exposome_explorer_ok = True
    except:
        exposome_explorer_ok = False

    return exposome_explorer_ok

def setup_exposome(databasefolder = "./DataBases", interactive=False):
    """
    Sets up the files relative to the Exposome Explorer database in the ``databasefolder``,
    splitting them for easier processing later on. If the session is set to be interactive,
    the user will be given time to add the files themselves; if not, the full suite of necessary
    files will be checked for their presence in ``databasefolder``

    Then, the "components" file will be splitted into one record oer line,
    as :obj:`~CanGraph.main` requires

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found
        interactive (bool): Whether the session is set to be interactive or not

    Returns:
        bool:
            True if everything went okay; False otherwise. If False,
            Exposome-Explorer should not be used as a data source
    """

    # Set the databasefolder to be an absolute path
    databasefolder = os.path.abspath(databasefolder)
    reldatabasedir = os.path.relpath(databasefolder)

    # If interactive, ask the user to add the E-E CSVs and whether they approve of file split.
    # If not, just check for the presence of the files and split components without problem
    exposome_explorer_ok = check_exposome_files()

    print("Setting up the Exposome Explorer database...")

    if interactive == True and exposome_explorer_ok == False:
        logging.info("First, please put the appropriate CSVs "
                     f"on the {reldatabasedir}/ExposomeExplorer path")
        time.sleep(1)
        logging.info("This database is confidential, and CANNOT be found online. "
                     "Please ask IARC for the files in case you need them.")
        time.sleep(1)
        if not os.path.exists(f"{databasefolder}/ExposomeExplorer"):
            logging.info("Let me create said folder for you..."); os.mkdir(f"{databasefolder}/ExposomeExplorer")
        print("Once you are ready, press [ENTER]", end=""); response = input()

        # If the user didn't add all the necessary files, exit
        exposome_explorer_ok = check_exposome_files()
        if exposome_explorer_ok == False:
            raise ValueError(f"Some files where not found in the {reldatabasedir}/ExposomeExplorer "
                             f"folder. Please, revise that it is correctly setup")
        else: logging.info("All checks OK")

        # If everything is ok, ask the user for split permission
        logging.info("Now, I will split the some files into a number of one-liner CSVs, simplifying import")
        print("Please, bear in mind that I will remove some of the the original files to avoid problems; is that ok? [Y/n]:", end="")
        response = input()
        if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
            print("You need to approve splitting for the setup to work. Exiting..."); time.sleep(1); exit(1)
    elif interactive == True and exposome_explorer_ok == True:
        logging.info(f"All files found on the {reldatabasedir}/ExposomeExplorer path. Moving on...")

    # Split the "components" files
    if exposome_explorer_ok:
        for filename in os.listdir(f"{databasefolder}/ExposomeExplorer/"):
            if "components" in filename:
                misc.split_csv(filename.split("/")[-1], f"{databasefolder}/ExposomeExplorer/")
        logging.info("Done!")

    return exposome_explorer_ok

def setup_hmdb(databasefolder = "./DataBases"):
    """
    Sets up the files relative to the HMDB database in the ``databasefolder``,
    splitting them for easier processing later on.

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found

    Returns:
        bool:
            True if everything went okay; False otherwise. If False,
            DrugBank should not be used as a data source
    """

    # Set the databasefolder to be an absolute path
    databasefolder = os.path.abspath(databasefolder)

    print("Setting up the Human Metabolome DataBase...")
    logging.info(f"I will download all the needed files and store them in {os.path.relpath(databasefolder)}/HMDB")
    if not os.path.exists(f"{databasefolder}/HMDB"): os.mkdir(f"{databasefolder}/HMDB")
    logging.info("Depending on your internet connection, this may take some time...")

    # The files in the HMDB database which we will download
    hmdb_urls = ["https://hmdb.ca/system/downloads/current/hmdb_proteins.zip",
                "https://hmdb.ca/system/downloads/current/urine_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/serum_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/csf_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/saliva_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/feces_metabolites.zip",
                "https://hmdb.ca/system/downloads/current/sweat_metabolites.zip"
                ]

    # For each file in the list, download and unzip them
    for url in hmdb_urls:
        filename = f"{url.split('/')[-1].split('.')[0]}.xml"
        splittag = "protein" if filename == "hmdb_proteins.xml" else "metabolite"

        try:
            misc.check_file(f"{databasefolder}/HMDB/{filename}")
            logging.info(f"File: {filename} was found on {databasefolder} and will not be re-downloaded")
        except:
            logging.info(f"Downloading file: {filename}...")
            misc.download_and_unzip(url, f"{databasefolder}/HMDB")

        logging.info(f"Unzipping file: {filename}...")
        misc.split_xml(os.path.abspath(f"{databasefolder}/HMDB/{filename}"), splittag, "hmdb")
        time.sleep(1) # Some waiting time to avoid overheating

    logging.info("Everything OK")

    return True

def setup_smpdb(databasefolder = "./DataBases"):
    """
    Sets up the files relative to the SMPDB database in the ``databasefolder``,
    splitting them for easier processing later on.

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found

    Returns:
        bool:
            True if everything went okay; False otherwise. If False,
            DrugBank should not be used as a data source
    """
    print("Setting up the Small Molecule Pathway DataBase...")
    logging.info("Since we have to download and unzip some files related, this may also take some time! :p")

    # Create the database folder if it does not previously exist
    if not os.path.exists(f"{databasefolder}/SMPDB"): os.mkdir(f"{databasefolder}/SMPDB")

    smpdb_urls = ["http://smpdb.ca/downloads/smpdb_pathways.csv.zip",
                "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip",
                "http://smpdb.ca/downloads/smpdb_proteins.csv.zip",
                "http://smpdb.ca/downloads/sequences/smpdb_protein.fasta.zip",
                "http://smpdb.ca/downloads/sequences/smpdb_gene.fasta.zip",
                ]

    # For each URL in the DataBase
    for url in smpdb_urls:

        # Calculate the basename of the files from SMPDB
        split = url.split('/')[-1].split('.')
        filename = ".".join(split[:2]); basename = split[0]

        # Assign if we want to check
        filetocheck = basename if basename in ["smpdb_metabolites", "smpdb_proteins"] else filename

        try:
            misc.check_file(f"{databasefolder}/SMPDB/{filetocheck}")
            logging.info(f"File: {filetocheck} was found on {databasefolder} and will not be re-downloaded")
        except:
            logging.info(f"Downloading and Unzipping: {filetocheck} ...")
            folder = basename if basename in ["smpdb_metabolites", "smpdb_proteins"] else ""
            misc.download_and_unzip(url, f"{databasefolder}/SMPDB/{folder}")

    # Since SMPDB already gives splitted files, there is no need for us to re-split! :p
    logging.info("All checks OK")

    return True

def setup_drugbank(databasefolder = "./DataBases", interactive = False):
    """
    Sets up the files relative to the SMPDB database in the ``databasefolder``,
    splitting them for easier processing later on.

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found
        interactive (bool): Whether the session is set to be interactive or not

    Returns:
        bool:
            True if everything went okay; False otherwise. If False,
            DrugBank should not be used as a data source
    """

    # Set the databasefolder to be an absolute path
    databasefolder = os.path.abspath(databasefolder)

    drugbank_ok = False # Initialize some databases

    print("Setting up the DrugBank database...")

    # Create the databasefolder. If created, just ignore its contents and use it (no problem as
    #overwrite is not likely, and, if it happens, we dont care much)
    if not os.path.exists(f"{databasefolder}/DrugBank"): os.mkdir(f"{databasefolder}/DrugBank")

    # If interactive, ask for DrugBank data to help the user download the data.
    # else, just check that the full_database file exists
    try:
        misc.check_file(f"{databasefolder}/DrugBank/full database.xml"); drugbank_ok = True
        logging.info(f"File: full database.xml was found on {databasefolder} and will not be re-downloaded")
    except:
        if interactive == True:
            logging.info("If you have a DrugBank account, I can help you automate the process!")
            print("Do you have an account [Y/n]", end=""); response = input()
            if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
                logging.error("Please, ask for a DrugBank account and run the script back. "
                            "The process must be approved by the DrugBank team, and may take a week or more.")
                exit(1)

            print("Please introduce your DrugBank username: ", end=""); username = input()
            print("Please introduce your DrugBank password: ", end=""); password = input()
            sys.stdout.write("\033[F"); print("Please introduce your DrugBank password: ********** HIDDEN **********")

            logging.info("Downloading the full database using curl. ")
            logging.info("Since its 1.4 GB, this may take some time...")
            if os.path.exists(f"{databasefolder}/DrugBank/full database.zip"):
                    os.remove(f"{databasefolder}/DrugBank/full database.zip")
            if os.path.exists(f"{databasefolder}/DrugBank/full database.xml"):
                    os.remove(f"{databasefolder}/DrugBank/full database.xml")
            try:
                os.system(f"curl -Lf --progress-bar -o \"{databasefolder}/DrugBank/full database.zip\" "
                          f"-u {username}:{password} https://go.drugbank.com/releases/5-1-9/downloads/all-full-database")
            except Exception as error:
                logging.error(f"cURL exited with error: {error}. "
                              f"Was your username/password OK?"); exit(1)

            logging.info("Extracting files from the full zip...")

            misc.unzip(f"{databasefolder}/DrugBank/full database.zip", f"{databasefolder}/DrugBank")
            os.remove (f"{databasefolder}/DrugBank/full database.zip")
            misc.check_file(f"{databasefolder}/DrugBank/full database.xml"); drugbank_ok = True

        else:
            drugbank_ok = False

    # If everything went OK, split the ``full_database`` file into its components ( 1 record per line )
    if drugbank_ok:
        logging.info("Splitting the contents of the full zip......")
        misc.split_xml(os.path.abspath(f"{databasefolder}/DrugBank/full database.xml"), "drug", "drugbank")
        logging.info("Everything OK")

    return drugbank_ok

def setup_database_index(databasefolder = "./DataBases"):
    """
    Prepares the index file for all the databases present in the ``databasefolder`` folder,
    which will helpfully reduce processing time a lot

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found

    Returns:
        dict:
            A dictionary containing the index for all the databases in ``databasefolder``.
            This index will be written as JSON in ``databasefolder``/index.json
    """
    print(f"Generating search index for all the files in {databasefolder}...")

    # Set the databasefolder to be an absolute path
    databasefolder = os.path.abspath(databasefolder)

    # First, we prepare a scan of all the files available on our "DataBases" folder
    # We will cycle through them later on to try and find the index keys
    all_files = misc.scan_folder(databasefolder)
    all_files = [ x for x in all_files if "index.json" not in x]

    # Prepare the dictionary which will hold the index
    index_dict = {"ChEBI_ID":{}, "HMDB_ID":{}, "InChI":{}, "MeSH_ID":{}, "Name":{}}
    inchi_regexp = "InChI\=1S?\/[A-Za-z0-9\.]+(\+[0-9]+)?(\/[cnpqbtmsih][A-Za-z0-9\-\+\(\)\,\/\?\;\.]+)*(\"|\<)"

    # Declare the progress bar we will use:
    with  alive_bar(len(all_files), title = "Generating Search Index...") as bar:

        # And run the indexer!
        for filepath in all_files:
            relpath = os.path.relpath(filepath, databasefolder)
            hmdb_ids, mesh_ids, inchis, chebi_ids, names = [], [], [], [], []
            with open(f'{filepath}', "r") as f:
                text = f.read()

            if "InChI" in text:
                results =  re.finditer(inchi_regexp, text)
                inchis = [x.group(0) for x in results]

            if "ExposomeExplorer/components" in relpath:
                component = pd.read_csv(os.path.abspath(filepath), dtype = str)
                chebi_ids = str(component["chebi_id"][0])
                names.append(component["name"][0])
                if str(component["alternative_names"][0]) != "nan":
                    names.extend(component["alternative_names"][0].split(";"))
            elif "SMPDB/smpdb_metabolites" in relpath:
                component = pd.read_csv(os.path.abspath(filepath), dtype = str)
                chebi_ids = list(component["ChEBI ID"])
                names.extend(list(component["Metabolite Name"]))
            elif "SMPDB/smpdb_proteins" in relpath:
                component = pd.read_csv(os.path.abspath(filepath), dtype = str)
                names.extend(list(component["Protein Name"]))
            elif "<chebi_id>" in text:
                chebi_ids = [ str(x.replace("<chebi_id>","").replace("</chebi_id>","")) for x in
                              list( re.findall("<chebi_id>.*</chebi_id>", text) ) ]

            if "<name>" in text:
                names.extend( [ x.replace("<name>","").replace("</name>","") for x in
                             list( re.findall("<name>.*</name>", text) ) ] )
            if "<synonym>" in text:
                names.extend( [ x.replace("<synonym>","").replace("</synonym>","") for x in
                             list( re.findall("<synonym>.*</synonym>", text) ) ] )
            if "DrugBank" in relpath:
                if "<kind>" in text:
                    names.extend( [ x.replace("Name</kind>\n      <value>","").replace("</value>","") for x in
                                list( re.findall("Name</kind>\n      <value>.*</value>", text) ) ] )
            elif "HMDB" in relpath:
                if "<iupac_name>" in text:
                    names.extend( [ x.replace("<iupac_name>","").replace("</iupac_name>","") for x in
                                list( re.findall("<iupac_name>.*</iupac_name>", text) ) ] )
                if "<traditional_iupac>" in text:
                    names.extend( [ x.replace("<traditional_iupac>","").replace("</traditional_iupac>","") for x in
                                list( re.findall("<traditional_iupac>.*</traditional_iupac>", text) ) ] )
                if "<uniprot_name>" in text:
                    names.extend( [ x.replace("<uniprot_name>","").replace("</uniprot_name>","") for x in
                                list( re.findall("<uniprot_name>.*</uniprot_name>", text) ) ] )

            if "<mesh-id>" in text:
                mesh_ids = [ x.replace("<mesh-id>","").replace("</mesh-id>","") for x in
                             list( re.findall("<mesh-id>.*</mesh-id>", text) ) ]

            if "HMDB" in text:
                hmdb_ids = list( re.findall("HMDB\d\d\d\d\d", text) )

            def create_and_append_to_index(item_list, item_type, relpath):
                for item in [x for x in item_list if x and str(x) != "nan"]:
                    index_dict[item_type].setdefault(item, [])
                    if relpath not in index_dict[item_type][item]:
                        index_dict[item_type][item].append(relpath)

            create_and_append_to_index(chebi_ids, "ChEBI_ID", relpath)
            create_and_append_to_index(hmdb_ids, "HMDB_ID", relpath)
            create_and_append_to_index(inchis, "InChI", relpath)
            create_and_append_to_index(mesh_ids, "MeSH_ID", relpath)
            create_and_append_to_index(names, "Name", relpath)

            bar()

    logging.info("Saving the index file...")
    with open(f"{databasefolder}/index.json", "w") as outfile:
        json.dump(index_dict, outfile)

    return index_dict

def setup_databases(databasefolder = "./DataBases", interactive = False):
    """
    Set Up the ``databasefolder`` from where the :obj:`~CanGraph.main` script will take its data. It
    does so by creating or removing and re-creating the ``databasefolder``, and putting inside it, or
    asking/checking if the user has put inside, the necessary files

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found
        interactive (bool): Whether the session is set to be interactive or not
    """
    setup_folders(databasefolder, interactive)
    exposome_explorer_ok = setup_exposome(databasefolder, interactive)
    hmdb_ok = setup_hmdb(databasefolder)
    smpdb_ok = setup_smpdb(databasefolder)
    drugbank_ok = setup_drugbank(databasefolder, interactive)

    print(f"The following databases have been set up in the {databasefolder} folder:")
    t = Texttable()
    t.add_rows([
                ['Exposome Explorer', 'HMDB', 'SMPDB', 'DrugBank', ],
                ["False" if item == 0 else "True" for item in [exposome_explorer_ok, hmdb_ok, smpdb_ok, drugbank_ok]]
                ])
    print(t.draw())

    logging.info("This is all! Since the WikiData and MetaNetX databases are auto-consulted ")
    logging.info("using RDF, the only thing we have left is generating the search index!")

    index_dict = setup_database_index(databasefolder, interactive)

    logging.info(f"The following number of identifiers have been indexed:")
    t = Texttable()
    t.add_rows([
                ['ChEBI_ID', 'HMDB_ID', 'InChI', 'MeSH_ID', "Name"],
                [ len(index_dict["ChEBI_ID"]), len(index_dict["HMDB_ID"]),
                  len(index_dict["InChI"]), len(index_dict["MeSH_ID"]),
                  len(index_dict["Name"])]
                ])
    logging.info(t.draw())

# ********* And, finally, the main function ********* #

def main():
    """
    The function that executes the code
    """

    # Parse the command line arguments
    # Done first in order to show errors if bad commands are issued
    parser = args_parser(); args = parser.parse_args()

    # Enable the --all argument
    if args.all: args.requirements = args.git = args.neo4j = args.databases = True

    # If the session is set to be interactive, display the logging messages
    logging.basicConfig(format='%(message)s')
    if args.interactive: logging.getLogger().setLevel(logging.INFO)

    # In any case, disable excessive verbosity on the neo4j driver
    logging.getLogger("neo4j").setLevel(logging.CRITICAL)

    # Only in interactive mode to avoid unnecessary waits
    if args.interactive:
        initial_message()
    else:
        print("Starting setup script")

    if args.requirements:
        install_packages(requirements_file = args.requirements)

    if args.git:
        setup_git(args.git)

    if args.neo4j:
        password = setup_neo4j(args.neo4j, args.neo4j_username, args.neo4j_password)
    else:
        password = ""

    if args.databases:
        setup_databases(args.databases, args.interactive)

    if args.interactive:
        final_message()
    else:
        print("The setup script has ended successfully")

    return password


if __name__ == '__main__':

    main()
