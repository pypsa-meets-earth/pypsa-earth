# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import datetime
import os
import sys

# sys.path.insert(0, os.path.abspath('scripts'))
# sys.path.insert(0, os.path.abspath('../scripts'))
# sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath("../scripts"))
for p in sys.path:
    print(p)

# -- Project information -----------------------------------------------------

project = "PyPSA meets Africa"
author = "Max Parzen"
copyright = f"{datetime.datetime.today().year}, {author}"

# The full version, including alpha/beta/rc tags
release = "0.0.1"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    'sphinx.ext.graphviz',
    # "sphinx.ext.pngmath",
    # "sphinxcontrib.tikz",
    # "rinoh.frontend.sphinx",
    "sphinx.ext.imgconverter",  # for SVG conversion
]

numpydoc_show_class_members = False

# autodoc_mock_imports leave the package out and does not require for building
# the documentation. If not mocked out errors can appear i.e. not automated
# documentation
autodoc_mock_imports = [
    "setuptools",
    # "esy",
    # "pypsa",
    # "numpy",
    # "pandas",
    # "geopandas",
    # "shapely",
    # "geoplot",
    # "matplotlib",
]
autodoc_default_options = {}  # {"members": None}
autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = [".rst", ".md"]
source_suffix = ".rst"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The master toctree document.
master_doc = "index"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "pypsa_meets_africa",
        "pypsa-meets-africa Documentation",
        author,
        "pypsa-meets-africa",
        "One line description of project.",
        "Miscellaneous",
    ),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "pypsa_meets_africa",
              "pypsa_meets_africa Documentation", [author], 1)]

# If true, show URL addresses after external links.
# man_show_urls = False

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {"https://docs.python.org/": None}
