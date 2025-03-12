# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
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
import shutil
import sys

from git import Repo

sys.path.insert(0, os.path.abspath("../scripts"))
for p in sys.path:
    print(p)

print(
    shutil.rmtree(
        "../../documentation",
        ignore_errors=True,
    )
)

print(
    Repo.clone_from(
        "https://github.com/pypsa-meets-earth/documentation", "../../documentation"
    )
)

print(
    shutil.rmtree(
        "img/",
        ignore_errors=True,
    )
)

print(
    shutil.copytree(
        "../../documentation/doc/img/",
        "img/",
        symlinks=False,
        ignore=None,
        copy_function=shutil.copy2,
        ignore_dangling_symlinks=False,
        dirs_exist_ok=False,
    )
)


# -- Project information -----------------------------------------------------

project = "PyPSA-Earth"
author = "Max Parzen"
copyright = f"{datetime.datetime.today().year}, {author}"

# The full version, including alpha/beta/rc tags
release = "0.6.0"

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.graphviz",
    "sphinx_copybutton",
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
    # "esy-osmfilter",
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
html_theme = "sphinx_book_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "repository_url": "https://github.com/pypsa-meets-earth/pypsa-earth",
    "use_repository_button": True,
    "show_navbar_depth": 1,
}

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "PyPSA-Earth"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "PyPSA-Earth"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "https://github.com/pypsa-meets-earth/pypsa-meets-earth.github.io/raw/main/assets/img/logo.png"

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "PyPSA-Earth",
        "PyPSA-Earth Documentation",
        author,
        "PyPSA-Earth",
        "One line description of project.",
        "Miscellaneous",
    ),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "pypsa_earth", "pypsa_earth Documentation", [author], 1)]

# If true, show URL addresses after external links.
# man_show_urls = False

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

copybutton_prompt_text = r"\.{3}/[^$]*\$ "  # Matches the pattern .../something $
copybutton_prompt_is_regexp = True
copybutton_only_copy_prompt_lines = True
