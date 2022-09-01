#!/usr/bin/env python3

# Copyright 2017-2020 Andrew Myers, Axel Huebl, Burlen Loring
# Maxence Thevenet, Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# -*- coding: utf-8 -*-
#
# WarpX documentation build configuration file, created by
# sphinx-quickstart on Tue Mar  7 22:08:26 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import subprocess
import sys
import urllib.request

import sphinx_rtd_theme

sys.path.insert(0, os.path.join( os.path.abspath(__file__), '../Python') )

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx_design',
    'breathe'
 ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'WarpX'
copyright = '2017-2021, WarpX collaboration'
author = 'WarpX collaboration'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = u'22.08'
# The full version, including alpha/beta/rc tags.
release = u'22.08'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

numfig = True

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'WarpXdoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'WarpX.tex', 'WarpX Documentation',
     'WarpX collaboration', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'warpx', 'WarpX Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'WarpX', 'WarpX Documentation',
     author, 'WarpX', 'WarpX is an advanced electromagnetic Particle-In-Cell code.',
     'Miscellaneous'),
]




# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://amrex-codes.github.io/': None}

# Setup the breathe extension
breathe_projects = {
    "WarpX": "../doxyxml/"
}
breathe_default_project = "WarpX"

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'
# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

# Download AMReX & openPMD-api Doxygen Tagfile to interlink Doxygen docs
url = 'https://amrex-codes.github.io/amrex/docs_xml/doxygen/amrex-doxygen-web.tag.xml'
urllib.request.urlretrieve(url, '../amrex-doxygen-web.tag.xml')

url = 'https://openpmd-api.readthedocs.io/en/latest/_static/doxyhtml/openpmd-api-doxygen-web.tag.xml'
urllib.request.urlretrieve(url, '../openpmd-api-doxygen-web.tag.xml')

# Build Doxygen
subprocess.call('cd ../; doxygen;'
                'mkdir -p source/_static;'
                'cp -r doxyhtml source/_static/;'
                'cp warpx-doxygen-web.tag.xml source/_static/doxyhtml/',
                shell=True)
