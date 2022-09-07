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
import neo4j                        # Import the whole package for error management
from neo4j import GraphDatabase     # The Neo4J python driver
import os, sys, shutil              # Vital modules to interact with the filesystem
import re                           # Work with regular expressions
from time import sleep              # Cute go slowly
from zipfile import ZipFile         # Work with ZIP files
import subprocess                   # Manage python sub-processes
import logging                      # Make ``verbose`` messages easier to show
import random, string               # For now, used to generate passwords
from git import Repo                # Manage Git directly from Python
import pip                          # Interact with the Python Package Index directly
import psutil                       # Kill the burden of the neo4j process
import argparse                     # Arguments pàrser for Python

import miscelaneous as misc         # A collection of useful functions

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
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--databases", action="store_true",
                        help="set up the databases from which the program will pull its info")
    parser.add_argument("-g", "--git", action="store_true",
                        help="prepare the git environment for the deploy script")
    parser.add_argument("-n", "--neo4j", action="store_true",
                        help="set up  the neo4j local environment, to run from ``./neo4j``")
    parser.add_argument("-r", "--requirements", action="store_true",
                        help=("installs all the requirements needed for all the possible options "
                              "of the program to work"))
    parser.add_argument("-i", "--interactive", action="store_true",
                        help=("tells the script if it wants interaction from the user "
                              "and information shown to them; similar to --verbose"))
    parser.add_argument("-a", "--all", action="store_true",
                        help=("runs all the options above at once; equivalent to -dgnr; "
                              "it DOES NOT activate the interactive mode"))

    parser.add_argument('databasefolder', nargs='?', default="DataBases",
                        help="the folder where the databases will be stored; by default ``./DataBases``")
    parser.add_argument('neo4j_home', nargs='?', default="neo4j",
                        help="the installation directory for the ``neo4j`` program; by default, ``./neo4j``")
    parser.add_argument('path_to_repo', nargs='?', default=".git",
                        help="the path to the Git repo; by default, ``.git``")

    # If no args are provided, show the help message
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser

def check_file(filepath):
    """
    Checks for the presence of a file, and exits the program

    Args:
        filepath (str): The path of the file its existence is being checked

    Raises:
        ValueError: If the file does not exist

    Returns:
        bool: True if successful, False otherwise.
    """
    if not os.path.exists(filepath):
        raise ValueError(f"Missing file: {filepath}. Please add the file and run the script back")
    return True

# ********* Add some verbose messages ********* #

def initial_message():
    """
    Prompts the user with an initial message if the session is set to be interactive.

    Args:
        interactive (bool): Whether the session is set to be interactive or not
    """
    logging.info("Hi! This setup script will guide you through the necessary steps to prepare ")
    logging.info("for the running of the main script \n")

    sleep(1)

    logging.info("This is because most of the databases we are using are confidential,")
    logging.info("so we cannot bundle them in the module. Also, we didn't include neo4j")
    logging.info("in the main repo to keep it light and due to copyright concerns \n")

    sleep(1)

    logging.info("Fortunately, you will only need to run the setup once if you use the --all ")
    logging.info("option; and, if everything is ready, you can just directly run main.py ")
    logging.info("(please read its README for more instructions on usage)\n")
    sleep(2)

def final_message(interactive = False):
    """
    Prompts the user with a final message.

    Args:
        interactive (bool): Whether the session is set to be interactive or not
    """

    logging.info("\n If you selected the --all option, you may now proceed to run the main script;")
    logging.info("in any other case, please make sure the parts relevant to the code you will ")
    logging.info("run are propperly configure (you may call ``python3 setup.py --help`` for more ")
    logging.info("info on available options)")

    exit(0)

# ********* Install the requirements ********* #

def install_packages(requirements_file = None, package_name = None):
    """
    Automates installing packages using PIP

    Args:
        requirements_file (str): The path to a "requirements.txt" file, containing one requirement per line
        package_name (str): A package to be installed

    Raises:
        ValueError: If neither a ``requirements_file`` nor a ``package_name`` is provided
    """
    logging.info("Checking for and installing packages using PIP...")
    if requirements_file != None:
        pip.main(["install", "-r", requirements_file])
    elif package != None:
        pip.main(["install", package_name])
    else:
        raise ValueError("You must give at least one requirement to install")

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

    logging.info("Done! Git is now configured :p")

# ********* Manage the Neo4J Database location, status and credentials ********* #

def kill_neo4j():
    """
    A simple function that kills any process that was started using a cmd argument including "neo4j"

    .. WARNING:: This function may unintendedly kill any command run from the ``neo4j`` folder.
        This is unfortunate, but the creation of this function was essential given that ``neo4j stop``
        does not work properly; instead of dying, the process lingers on, interfering
        with :obj:`~CanGraph.setup.find_neo4j_installation_status` and hindering the main program
    """
    for proc in psutil.process_iter():
        if "neo4j" in " ".join(proc.cmdline()):
            logging.info("Killing existing Neo4j sessions...")
            proc.kill()

def find_neo4j_installation_status(neo4j_home = "neo4j"):
    """
    Finds the installation status of Neo4J by trying to use it normally, and analyzing any thrown exceptions

    Args:
        neo4j_home (str): the installation directory for the ``neo4j`` program; by default, ``neo4j``
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
        logging.info("Now, lets start it and check if it presents the default user/passwd")

        # If no exceptions are thrown, then it exists! Lets try to start it then
        neo4j_exists = True; sleep(8)

        instance = "bolt://localhost:7687"; user = "neo4j"; passwd = "neo4j"
        driver = GraphDatabase.driver(instance, auth=(user, passwd))

        # We will try to get the Import Path just as a test to see if the default auths are present
        Neo4JImportPath = misc.get_import_path(driver)

    # If it throws an ``AuthError``, its because the credentials are not the default
    except neo4j.exceptions.AuthError as error:
        logging.info("Neo4J credentials are not default. We will need to re-install...")
        logging.info("Let me stop neo4j first...")
        default_credentials = False

    # If they are, it throws a ``ClientError``, because it requires a password change first
    except neo4j.exceptions.ClientError as error:
        logging.info("Installation with default credentials found! I'll work with that :p")
        default_credentials = True

    # If the executable is not found in ``neo4j_home``. we will need to install it
    except FileNotFoundError as error:
        logging.info("No Neo4J installation has been found. Defaulting to installing from website...")
        neo4j_exists = False

    # In any other case, re-install it, too
    except Exception as error:
        logging.info("An unknown error has been raised. Defaulting to installing from website...")
        neo4j_exists = False

    kill_neo4j()

    return neo4j_exists, default_credentials

def install_neo4j(neo4j_home = "neo4j", interactive = False):
    """
    Installs the neo4j database program in the ``neo4j_home`` folder, by getting it from the internet
    according to the Operating System the script is been run in (aims for multi-platform!)

    Args:
        neo4j_home (str): the installation directory for the ``neo4j`` program; by default, ``neo4j``
        interactive (str): tells the script if it wants interaction from the user and information shown to them
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

    logging.info("Installing neo4j in local folder...")
    if  sys.platform == "linux" or sys.platform == "linux1" or sys.platform == "linux2" or sys.platform == "darwin":
        misc.download_and_untargz("https://neo4j.com/artifact.php?name=neo4j-community-4.4.10-unix.tar.gz", "/tmp/")
    elif  sys.platform == "win32":
        misc.download_and_untargz("https://neo4j.com/artifact.php?name=neo4j-community-4.4.10-windows.zip", "/tmp/")
    shutil.move("/tmp/neo4j-community-4.4.10/", neo4j_home)
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

def configure_neo4j(neo4j_home):
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

    sp = subprocess.run(["neo4j/bin/neo4j-admin", "memrec"], stdout=subprocess.PIPE)
    sp.check_returncode() # Raise error if error is found

    for line in str(sp.stdout).split("\\n"):
        if line.startswith("dbms."):
            key = line.split("=")[0]; val = line.split("=")[1]
            update_neo4j_confs(key, val)

    # Enable the APOC plugin
    logging.info("Now, lets configure the neo4j environment")
    logging.info("First, let me enable the APOC plugin, which we will use heavily...")
    if not os.path.exists(f"{neo4j_home}/plugins/apoc-4.4.0.7-core.jar"):
        shutil.copyfile(f"{neo4j_home}/labs/apoc-4.4.0.7-core.jar", f"{neo4j_home}/plugins/apoc-4.4.0.7-core.jar")

    # And modify the ``conf`` file
    logging.info("And, now, I will modify the conf file to update the preferences")
    update_neo4j_confs("apoc.import.file.enabled", True)
    update_neo4j_confs("apoc.export.file.enabled", True)
    update_neo4j_confs( "apoc.http.timeout.connect", 60000)
    update_neo4j_confs("apoc.http.timeout.read", 180000)
    update_neo4j_confs("dbms.memory.transaction.global_max_size", "1024m")

    # Force the "import" path to be absolute:
    update_neo4j_confs("dbms.directories.import", f"{neo4j_home}/import")

def change_neo4j_password(neo4j_home, new_password, old_password = "neo4j", user = "neo4j", database = "system"):
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
    sleep(8) # Enable the DBMS to start

    # Alter the password using cypher-shell
    sp = subprocess.run([f"{neo4j_home}/bin/cypher-shell", "-d", f"{database}", "-u", f"{user}",
                        "-p", f"{old_password}",
                        f"ALTER CURRENT USER SET PASSWORD FROM '{old_password}' TO '{new_password}'"],
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    sp.check_returncode() # Raise error if error is found

    kill_neo4j() # Kill the Neo4J process

def setup_neo4j(neo4j_home = "neo4j", interactive = False):
    """
    Sets ups the neo4j environment in ``neo4j_home``, so that the functions in :obj:`~CanGraph.main` can propperly
    function. Using the functions present in this module, it finds if neo4j is installed with default credentials,
    and, if not, it installs it, changing the default password to a new one, and returning its value

    Args:
        neo4j_home (str): the installation directory for the ``neo4j`` program; by default, ``neo4j``
        interactive (str): tells the script if it wants interaction from the user and information shown to them

    Returns:
        str: The password that was set up for the new neo4j database

    .. NOTE:: This has been designed to be used with a ``neo4j_home`` located in the WorkDir, but can be used
        in any other location with read/write access, or even with ``apt install`` installed versions! Just
        find its ``neo4j_home``, make sure it has r/w access, and provide it to the program!
    """

    # Fist, we need to kill any existing neo4j instances. As explained in :obj:`~CanGraph.setup.kill_neo4j`,
    #this is because the process may be lingering on otherwise, making the program crash
    kill_neo4j()

    # Display an initial message
    logging.info("Now, lets set up the Neo4J environment")
    logging.info("First, I will check for an existing Neo4J installation with the default user/passwd, as this might save time.")

    # Set the variables that will help us find the status of neo4j's installation to ``False`` by default
    neo4j_exists = False; default_credentials = False

    # Find out the status of neo4j installation
    neo4j_exists, default_credentials = find_neo4j_installation_status(neo4j_home)

    # If needed, install neo4j
    if neo4j_exists == False or default_credentials == False:
        install_neo4j(neo4j_home, interactive)

    configure_neo4j(neo4j_home) # And configure it

    # Generate a random password for the neo4j database
    logging.info("Finally, let me generate a password for the neo4j database and establish it")
    characters = string.ascii_letters + string.digits
    password = ''.join(random.choice(characters) for i in range(12))
    if interactive == True: print(f"Neo4J's password has been set to: {password}")

    # And change it, making the program finish :p
    change_neo4j_password(neo4j_home, password)

    logging.info("Changed password. The neo4j environment is ready for use")

    return password

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
                   "May we use it for the project? (We may ovewrite its contents!")
            print("Use existing folder? [Y/n]", end=""); response = input()
        else:
            response = "n"

        if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
            print("Please fix the folder situation. Exiting..."); sleep(1); exit(1)
        shutil.rmtree(databasefolder)
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

    # Check for the presence of all files; if :obj:`CanGraph.setup.check_file` throws an error, set exposome_explorer_ok to ``False``
    try:
        check_file(f"{databasefolder}/ExposomeExplorer/units.csv");                check_file(f"{databasefolder}/ExposomeExplorer/subjects.csv")
        check_file(f"{databasefolder}/ExposomeExplorer/specimens.csv");            check_file(f"{databasefolder}/ExposomeExplorer/samples.csv")
        check_file(f"{databasefolder}/ExposomeExplorer/reproducibilities.csv");    check_file(f"{databasefolder}/ExposomeExplorer/publications.csv")
        check_file(f"{databasefolder}/ExposomeExplorer/microbial_metabolite_identifications.csv")
        check_file(f"{databasefolder}/ExposomeExplorer/metabolomic_associations.csv"); check_file(f"{databasefolder}/ExposomeExplorer/cancer_associations.csv")
        check_file(f"{databasefolder}/ExposomeExplorer/measurements.csv");         check_file(f"{databasefolder}/ExposomeExplorer/experimental_methods.csv")
        check_file(f"{databasefolder}/ExposomeExplorer/correlations.csv");          check_file(f"{databasefolder}/ExposomeExplorer/components.csv")
        check_file(f"{databasefolder}/ExposomeExplorer/cohorts.csv");              check_file(f"{databasefolder}/ExposomeExplorer/cancers.csv")
        exposome_explorer_ok = True
    except:
        logging.error(f"Some files where not found in the {databasefolder}/ExposomeExplorer folder. Please, revise that it is correctly setup")
        exposome_explorer_ok = False

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

    exposome_explorer_ok = False # Initialize variable

    # If interactive, ask the user to add the E-E CSVs and whether they approve of file split.
    # If not, just check for the presence of the files and split components without problem
    if interactive == True:

        logging.info("First, please put the CSVs pertaining the exposome-explorer database on the ./DataBases/ExposomeExplorer path")
        sleep(1)
        logging.info("This database is confidential, and CANNOT be found online. Please ask IARC for the files in case you need them.")
        sleep(1)
        if not os.path.exists(f"{databasefolder}/ExposomeExplorer"): logging.info("Let me create said folder for you..."); os.mkdir(f"{databasefolder}/ExposomeExplorer")
        print("Once you are ready, press [ENTER]", end=""); response = input()

        exposome_explorer_ok = check_exposome_files()

        logging.info("All checks OK")

        logging.info("Now, I will split the some files into a number of one-liner CSVs, simplifying import")
        print("Please, bear in mind that I will remove some of the the original files to avoid problems; is that ok? [y/N]:", end="")
        response = input()
        if response != "y" and response != "Y" and response != "Yes" and response != "yes" and response != "YES":
            print("You need to approve splitting for the setup to work. Exiting..."); sleep(1); exit(1)
    else:
        exposome_explorer_ok = check_exposome_files()

    # Split the "components" files
    if exposome_explorer_ok:
        for filename in os.listdir(f"{databasefolder}/ExposomeExplorer/"):
            if "components" in filename:
                misc.split_csv(filename.split("/")[-1], f"{databasefolder}/ExposomeExplorer/")

    return exposome_explorer_ok

def setup_hmdb(databasefolder = "./DataBases"):
    """
    Sets up the files relative to the HMDB database in the ``databasefolder``,
    splitting them for easier processing later on.

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found
    """

    # Set the databasefolder to be an absolute path
    databasefolder = os.path.abspath(databasefolder)

    logging.info("Now, lets go on with the Human Metabolome DataBase. I will download all "
                 f"the needed folders and store them in {os.path.relpath(databasefolder)}/HMDB")
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
        logging.info(f"Downloading and Unzipping: {filename}...")
        misc.download_and_unzip(url, f"{databasefolder}/HMDB")
        logging.info("File downloaded, now splitting its contents...")
        if filename == "hmdb_proteins.xml":
            misc.split_xml(os.path.abspath(f"{databasefolder}/HMDB/{filename}"), "protein", "hmdb")
        else:
            misc.split_xml(os.path.abspath(f"{databasefolder}/HMDB/{filename}"), "metabolite", "hmdb")
        sleep(1) # Some waiting time to avoid overheating

    logging.info("Everything OK")

def setup_smpdb(databasefolder = "./DataBases"):
    """
    Sets up the files relative to the SMPDB database in the ``databasefolder``,
    splitting them for easier processing later on.

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found
    """
    logging.info("Super! Now, lets download and unzip the files related to the Small Molecule Pathway DataBase!")
    logging.info("This may also take some time")

    if not os.path.exists(f"{databasefolder}/SMPDB"): os.mkdir(f"{databasefolder}/SMPDB")

    smpdb_urls = ["http://smpdb.ca/downloads/smpdb_pathways.csv.zip",
                "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip",
                "http://smpdb.ca/downloads/smpdb_proteins.csv.zip",
                "http://smpdb.ca/downloads/sequences/smpdb_protein.fasta.zip",
                "http://smpdb.ca/downloads/sequences/smpdb_gene.fasta.zip",
                ]

    for url in smpdb_urls:
        split = url.split('/')[-1].split('.')[0]
        logging.info(f"Downloading and Unzipping: {split}...")
        if split in ["smpdb_metabolites", "smpdb_proteins"]:
            misc.download_and_unzip(url, f"{databasefolder}/SMPDB/{split}")
        else:
            misc.download_and_unzip(url, f"{databasefolder}/SMPDB/")

    check_file(f"{databasefolder}/SMPDB/smpdb_pathways.csv");
    check_file(f"{databasefolder}/SMPDB/smpdb_proteins/");            check_file(f"{databasefolder}/SMPDB/smpdb_metabolites/")
    check_file(f"{databasefolder}/SMPDB/smpdb_gene.fasta");           check_file(f"{databasefolder}/SMPDB/smpdb_protein.fasta")

    # Since SMPDB already gives splitted files, there is no need for us to re-split! :p

    logging.info("All checks OK")

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

    logging.info("We are almost there! Now, we have to prepare the DrugBank database. "
                 "If you have a DrugBank account, I can help you automate the process!")

    # Create the databasefolder. If created, just ignore its contents and use it (no problem as
    #overwrite is not likely, and, if it happens, we dont care much)
    if not os.path.exists(f"{databasefolder}/DrugBank"): os.mkdir(f"{databasefolder}/DrugBank")

    # If interactive, ask for DrugBank data to help the user download the data.
    # else, just check that the full_database file exists
    if interactive == True:
        print("Do you have an account [Y/n]", end=""); response = input()
        if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
            logging.error("Please, ask for a DrugBank account and run the script back. "
                          "The process must be approved by the DrugBank team, and may take a week or more.")
            exit(1)

        print("Please introduce your DrugBank username: ", end=""); username = input()
        print("Please introduce your DrugBank password: ", end=""); password = input()
        sys.stdout.write("\033[F"); print("Please introduce your DrugBank password: ********** HIDDEN **********")

        logging.info("Downloading the full database using curl. Since its 1.4 GB, this may take some time...")
        if os.path.exists(f"{databasefolder}/DrugBank/full_database.zip"): os.remove(f"{databasefolder}/DrugBank/full_database.zip")
        if os.path.exists(f"{databasefolder}/DrugBank/full_database.xml"): os.remove(f"{databasefolder}/DrugBank/full_database.xml")
        try:
            os.system(f"curl -Lf --progress-bar -o \"{databasefolder}/DrugBank/full_database.zip\" -u {username}:{password} https://go.drugbank.com/releases/5-1-9/downloads/all-full-database")
            check_file(f"{databasefolder}/DrugBank/full_database.xml"); drugbank_ok = True
        except:
            logging.error("Something went wrong with CURL. Was your username/password OK?"); exit(1)
    else:
        try:
            check_file(f"{databasefolder}/DrugBank/full_database.xml"); drugbank_ok = True
        except:
            drugbank_ok = False

    # If everything went OK, split the ``full_database`` file into its components ( 1 record per line )
    if drugbank_ok:
        logging.info("Extracting files from the full zip..")

        misc.unzip(f"{databasefolder}/DrugBank/full_database.zip")
        logging.info("File extracted, now splitting its contents...")
        misc.split_xml(os.path.abspath(f"{databasefolder}/DrugBank/full database.xml"), "drug", "drugbank")
        os.remove(os.path.abspath(f"{databasefolder}/DrugBank/full_database.zip"))

        logging.info("Everything OK")

    return drugbank_ok

def databases(databasefolder = "./DataBases", interactive = False):
    """
    Set Up the ``databasefolder`` from where the :obj:`~CanGraph.main` script will take its data. It
    does so by creating or removing and re-creating the ``databasefolder``, and putting inside it, or
    asking/checking if the user has put inside, the necessary files

    Args:
        databasefolder (str): The main folder where all the databases we will be using are to be found
        interactive (bool): Whether the session is set to be interactive or not

    WARNING FIXEATE LOS MENSAJES INICIALES Y FINALES
    """
    setup_folders(databasefolder, interactive)
    setup_exposome(databasefolder, interactive)
    setup_hmdb(databasefolder)
    setup_smpdb(databasefolder)
    setup_drugbank(databasefolder, interactive)

    logging.info("This is all! Since the WikiData and MetaNetX databases are auto-consulted ")
    logging.info("using RDF, there is no need to set up anything more!")

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

    # Only in interactive mode to avoid unnecessary waits
    if args.interactive: initial_message()

    if args.requirements:
        install_packages(requirements_file = "requirements.txt")

    if args.git:
        setup_git(args.path_to_repo)

    if args.neo4j:
        setup_neo4j(args.neo4j_home)

    if args.databases:
        setup_databases(args.databasefolder, args.interactive)

    if args.interactive: final_message()


if __name__ == '__main__':

    main()
