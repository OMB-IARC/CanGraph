# Configuration file for the Sphinx documentation builder.

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath('../ExposomeExplorer'))
sys.path.insert(0, os.path.abspath('../GraphifySMPDB'))
sys.path.insert(0, os.path.abspath('../GraphifyHMDB'))
sys.path.insert(0, os.path.abspath('../GraphifyDrugBank'))
sys.path.insert(0, os.path.abspath('../QueryWikidata'))
sys.path.insert(0, os.path.abspath('../MeSHandMetaNetX'))


# -- Project information -----------------------------------------------------

project = 'CanGraph'
copyright = '2022, Pablo Marcos <software @ loreak . org>'

version = '0.9'  # Short release number
release = '0.9'  # Full release number


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings.
extensions = [
   'sphinx.ext.napoleon',         # Enable parsing of NumPy/Google style docstrings
   'sphinx.ext.autodoc',          # Auto-document the code
   'sphinx.ext.autosummary',      # Auto-create the summaries
   'sphinx.ext.todo',             # Enable the TODO and TODOLIST directives
   'sphinx.ext.intersphinx',      # Link to other projects’ documentation
   'sphinx_design',               # Enable cute docs layouts
   'myst_parser',                 # Enable MarkDown Documentation
]

# Enable correct parsing by myst and remove superflous MD warnings
myst_enable_extensions = ["colon_fence"]
suppress_warnings = ["myst.header"]

# Autogenerate summaries
autosummary_generate = True

# Remove bloat by ignoring module names
add_module_names = False

# Include TO-DOs
todo_include_todos = False

# Set the documents which Sphinx will process
source_suffix = ['.rst', '.md']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'modules.rst']

# Check internal links to modules
nitpicky = True

# Enable linking to external project's documentation
intersphinx_mapping = {'neo4j': ('https://neo4j.com/docs/api/python-driver/current/', None),
                       #"python": ("https://docs.python.org/3.8", None),
                       #"sphinx": ("https://www.sphinx-doc.org/en/master", None)
                       }


# -- Options for HTML output -------------------------------------------------

#html_title =
html_theme = 'sphinx_book_theme' # .. NOTE:: REPLACE "-" WITH "_" HERE
html_logo = "_static/cangraph-logo.png"
html_favicon = "_static/cangraph-mini.png"
html_last_updated_fmt = ""

html_theme_options = {
  "logo_only": True,
  "home_page_in_toc": True,
  "use_download_button": True,
  "repository_url": "https://github.com/OMB-IARC/CanGraph/",
  "use_repository_button": True,
  "use_issues_button": True,
  "use_edit_page_button": True,
  "repository_branch": "main",
  "path_to_docs": "docs",
  "toc_title": "In this page",
  "extra_navbar": '''Logo by <a href="https://www.instagram.com/danonino.caducado/"> Daniel Marcos </a> </br>
                     Website by <a href="https://www.pablomarcos.me/"> Pablo Marcos </a>''',
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
