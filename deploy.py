#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

"""
A python module that simplifies deploying the different formats the program can present itself as:

* A web documentation that is made using sphinx
* A PDF manual likewise made, also using LaTeX
* A git repo where the program can be accessed and version-tagged
* And, in the future... a singularity container!

CanGraph.deploy Usage
---------------------

To use this module:

.. argparse::
   :module: CanGraph.deploy
   :func: args_parser
   :prog: python3 deploy.py
   :nodefault:

.. NOTE:: For this program to work, the Git environment **has to be set up first**.
    You can ensure this by using: :obj:`CanGraph.setup.setup_git`

CanGraph.deploy Functions
-------------------------

This module is comprised of:
"""

# Import external modules necessary for the script
from git import Repo    # Manage Git directly from Python
import argparse         # Arguments pàrser for Python
import os, sys, shutil  # Vital modules to interact with the filesystem

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
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--main", help="deploy code to your dev branch", action="store_true")
    parser.add_argument("-d", "--dev", help="deploy code to your main branch", action="store_true")
    parser.add_argument("-w", "--web", help="deploys the documentation to the 'pages' web site, depending on which is activated", action="store_true")
    parser.add_argument("-p", "--pdf", help="generates a PDF version of the sphinx manual, and saves it to the repo", action="store_true")

    # If no args are provided, show the help message
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser

def git_push(path_to_repo, remote_names, commit_message, force = False):
    """
    Pushes the current repo's state and current branch to a remote git repository

    Args:
        path_to_repo (str): The path to the local ``.git`` folder
        remote_names (list or str): The names of the remote to which we want to commit, which must be previously
            configured (see :obj:`CanGraph.setup.setup_git`). e.g.: ["github", "codeberg"]
        commit_message (str): The Git Commit Message for the current repo's state
        force (bool): Whether to force the commit (necessary if you are resetting the HEAD)

    .. NOTE:: ``gitpython`` is not good at managing complex commit messages (i.e. those with a Subject and a Body).
        If you want to add one of those, please, use ``\\n`` as the separator; the function will take care of the rest

    .. seealso:: The approach taken hare was inspired
        by `StackOverflow #41836988 <https://stackoverflow.com/questions/41836988/git-push-via-gitpython>`_
    """
    repo = Repo(path_to_repo)
    repo.git.add('.')

    current_branch = repo.active_branch.name
    repo.index.commit(commit_message.replace("\\n", "\n"))

    i = 1
    while i <= 3:
        try:
            if isinstance(remote_names, list):
                for remote in remote_names:
                    origin = repo.remote(name=remote)
                    origin.push(refspec=f"{current_branch}:{current_branch}", force=force)
            elif isinstance(remote_names, str):
                origin = repo.remote(name=remote_names)
                origin.push(refspec=f"{current_branch}:{current_branch}", force=force)
            break
        except Exception as e:
            i += 1
            print(f"Error: {e} \nPlease try again ({i}/3)...")

def make_sphinx_prechecks(docs_folder = "./docs/", work_dir = ".", gen_apidocs = False):
    """
    Generates sphinx api-docs for automatic documentation and uses ``make linkcheck``` to check for broken links

    Args:
        gen_apidocs (bool): Whether to re-generate the API docs. Default is ``False`` since we use MD (we dont want RST files)
        docs_folder (str): The path to sphinx's docs folder, where the tests will be run; by default ``./docs/``
        work_dir (str): The current Working Directory; by default ``.``
    """
    cwd = os.getcwd(); os.chdir(os.path.abspath(docs_folder))
    if gen_apidocs == True:
        print("First, lets generate all the needed RST auto-doc files that don't already exist")
        os.system("sphinx-apidoc -o . ..")
        print("If everything is OK, press [ENTER]", end=""); response = input()

    print("Now, I will use sphinx's make linkcheck to list all broken links for you to fix them")
    os.system("make linkcheck")
    print("If everything is OK, press [ENTER]", end=""); response = input()
    os.chdir(cwd)

def deploy_webdocs(docs_folder = "./docs/", work_dir = ".", prechecks_done=False, custom_domain=None):
    """
    Generates the HTML web docs, and publishes it both to Github and Codeberg pages

    Args:
        docs_folder (str): The path to sphinx's docs folder, where the tests will be run;
            by default ``./docs/``
        work_dir (str): The current Working Directory; by default
        prechecks_done (bool): Whether the prechecks present
            in `~CanGraph.deploy.make_sphinx_prechecks` have already been made
        custom_domain (str): A custom domain to deploy de docs to.

    Returns:
        bool: Whether the prechecks have already been done; always True if the function is run

    .. NOTE:: For ``custom_domain`` to work, please configure your DNS records apparently
    .. NOTE:: ``modules.rst`` is not removed, but it is correctly ignored in conf.py
    """

    print("You have decided to update the project's web documentation.")

    # Avoid double pre-checking with sphinx-linkcheck and sphink-apidoc
    if prechecks_done == True:
        print("I see you have already checked sphinx's linkcheck and generated the RST auto-doc files")
        print("Let's skip to sphinx's make html")

    else:
        make_sphinx_prechecks()

        print("Finally, I will use sphinx's make html to generate the web docs")

    # Make HTML documentation
    cwd = os.getcwd(); os.chdir(os.path.abspath(docs_folder)); os.system("make html")
    print("If everything is OK, press [ENTER]", end=""); response = input()

    print("Done! Now, lets deploy the docs to the Web Repository using Git")

    os.chdir(cwd) # Remember to go back to the workdir!

    # First, we copy the "docs" folder to a higher folder so that it stops being affected by Git
    if os.path.exists(os.path.abspath("../html")):
        print("I need to use ../html to preserve files from Git, but path ../html already exists.")
        print("Allow me to remove it? [Y/n]", end=""); response = input()
        if response == "n" or response == "N" or response == "No" or response == "no" or response == "NO":
            raise ValueError("Please fix the folder situation. Exiting..."); exit(1)
        else:
            print("Saving contents from git...")
            shutil.rmtree("../html")
    shutil.copytree(os.path.abspath("./docs/_build/html/"), os.path.abspath("../html"))

    # Then, we save the git repo's current state before checkout
    repo = Repo(".git")
    repo.git.add('.')
    repo.git.stash('save')

    # Annotating the current branch name for the future
    current_branch = repo.active_branch.name

    print("Let me clean the 'pages' repo first...")

    # And checkout to the pages branch
    repo.git.checkout('pages')

    # Where we remove previous contents
    for item in os.listdir("."):
        if item not in ".git":
            if os.path.isfile(item):
                os.remove(item)
            elif os.path.isdir(item):
                shutil.rmtree(item)

    # Reset to the previous commit to make the repo smaller
    repo.git.add('.'); repo.git.reset('HEAD~1')

    # And then copy the contents of the saved "html" folder
    shutil.copytree(os.path.abspath("../html/"), os.path.abspath("./"), dirs_exist_ok=True)

    # [OPTIONAL] And we add some files to enable accessing from domains
    if custom_domain:
        # On Codeberg
        with open(".domains", "w") as domains_file:
            domains_file.write(("cangraph.pablomarcos.me \n"
                                "cangraph.flyingflamingo.codeberg.page \n"
                                "pages.cangraph.flyingflamingo.codeberg.page \n"))

        # And on GitHub, too:
        with open("CNAME", "w+") as CNAME_file:
            CNAME_file.write("""cangraph.pablomarcos.me""")

    # Add .nojekyll to force github to correctly process the pages branch
    with open(".nojekyll", "w+") as nojekyll:
            nojekyll.write("\n")

    # The repo is ready! Lets commit:
    print("Now, please provide a comment for Git ...")
    commit_message = input()

    print("Uploading to Git Repos ...")
    git_push(".git", ["github", "codeberg"], commit_message, force = True)

    # And leave everything as it was
    print("Site has been deployed")
    shutil.rmtree(os.path.abspath("../html/"))

    repo.git.checkout(current_branch)
    repo.git.stash('apply')

    # Prechecks have been done
    return True

def deploy_pdf_manual(docs_folder = "./docs/", work_dir = ".", manual_location = "./CanGraph_Manual.pdf",
                      prechecks_done=False, custom_domain=None):
    """
    Parses the command line arguments into a more usable form, providing help and more

    Generates the PDF docs guide, and publishes it in ``manual_location``

    Args:
        docs_folder (str): The path to sphinx's docs folder, where the tests will be run;
            by default ``./docs/``
        work_dir (str): The current Working Directory; by default
        prechecks_done (bool): Whether the prechecks present
            in `~CanGraph.deploy.make_sphinx_prechecks` have already been made
        custom_domain (str): A custom domain to deploy de docs to.
        manual_location (str): The location (**including filename**) of the finalised PDF manual,
            relative to the location of the script

    Returns:
        bool: Whether the prechecks have already been done; always True if the function is run
    """

    # Avoid double pre-checking with sphinx-linkcheck and sphink-apidoc
    if prechecks_done == True:
        print("I see you have already checked sphinx's linkcheck and generated the RST auto-doc files")
        print("Let's skip to sphinx's make latexpdf")

    else:
        make_sphinx_prechecks()

        print("Finally, I will use sphinx's make latexpdf to generate the web docs")

    # Generate the PDF docs
    print("Now, let me use sphinx's latexpdf to generate a PDF doc")
    cwd = os.getcwd(); os.chdir(os.path.abspath(docs_folder)); os.system("make latexpdf > /dev/null")
    print("If everything is OK, press [ENTER]", end=""); response = input()

    # And copy them in the appropriate place
    os.chdir(os.path.abspath(cwd))
    if os.path.exists(manual_location):
        os.remove(manual_location)
    shutil.copyfile(os.path.abspath("./docs/_build/latex/cangraph.pdf"), os.path.abspath(manual_location))
    print("The PDF Manual has been copied in the Repo's Root Folder (CanGraph_Manual.pdf)")

    # Prechecks have been done
    return True

def deploy_code(branch = "dev"):
    """
    Deploys code from a given branch to the corresponding remote.

    Args:
        branch (str): The name of the branch of the **local** git repo that we want to deploy

    .. NOTE:: Normally, the ``pages`` branch should be published using :obj:`~CanGraph.deploy.deploy_webdocs`,
        which in theory would be weird to publish without updating the docs first
    """

    # Ask for a Git Commit Message
    print(f"You have decided to push the updated code to your {branch} branch.")
    print("First, please provide a comment for Git ...")
    commit_message = input()

    repo = Repo('.git')

    # Annotating the current branch name for the future
    current_branch = repo.active_branch.name

    if branch == "main":
        repo.git.add('.')
        repo.git.stash ("save")
        repo.git.checkout("main")
        repo.git.merge("dev")

        commit_message = "🔀 Merged advances from dev into main"

        git_push(".git", ["github", "codeberg"], commit_message, force = True)

        # Finally, we reset the dev branch
        repo.git.branch("-D", current_branch)
        repo.git.checkout("-b", current_branch)

    if branch == "dev":
        # And git add / commit / push
        print("Uploading to Git Repos ...")
        git_push(".git", ["github", "codeberg"], commit_message, force = True)


def main():
    """
    The function that executes the code
    """

    # Parse the command line arguments
    # Done first in order to show errors if bad commands are issued
    parser = args_parser(); args = parser.parse_args()

    print("Hi! This is a helper program intended to ease the burden of maintaining Codeberg \n"
          "and GitHub mirror repos for the CanGraph project, as well as some other things, \n"
          "such as the Singularity Container or the PDF Manuals")

    prechecks_done = False # No prechecks are done at first, of course

    if os.path.exists("./docs/_build/"):
        print("Let me remove any cached files in the '_build' folder...")
        shutil.rmtree("./docs/_build/")

    # Do what we were asked to do
    if args.web:
        prechecks_done = deploy_webdocs(prechecks_done=prechecks_done)

    if args.pdf:
        prechecks_done = deploy_pdf_manual(prechecks_done=prechecks_done)

    if args.dev:
        deploy_code(branch = "dev")

    if args.main:
        deploy_code(branch = "main")

    print("That was it! Thanks for using this helper script :p")


if __name__ == '__main__':

    main()
