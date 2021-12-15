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
import os
import sys
from pygments.lexer import RegexLexer
from pygments import token
from sphinx.highlighting import lexers

sys.path.insert(0, os.path.abspath('../'))

# -- Pygmentize GMAK input file ----------------------------------------------

class GmiLexer(RegexLexer):
    name = 'gmi'
    tokens = {
        'root': [
            (r'^\$\S+$', token.Keyword.Namespace),
            (r'^\S+', token.Name.Variable),
            (r'\S+', token.Text),
            (r'\s+', token.Text),
        ]
    }

class TreeLexer(RegexLexer):
    name = 'tree'
    tokens = {
        'root': [
            (r'%[a-zA-Z0-9-]+', token.Keyword.Variable),
            (r'\s+', token.Text),
            (r'\S+_', token.Text),
            (r'\S+', token.Text),
        ]
    }

# -- Pygmentize file trees ---------------------------------------------------

lexers['gmi'] = GmiLexer()
lexers['tree'] = TreeLexer()


# -- Project information -----------------------------------------------------

project = 'gridmaker'
copyright = '2021, MSSM/LabMMol'
author = 'MSSM/LabMMol'

# The full version, including alpha/beta/rc tags
release = '1.0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# Add the extension
extensions = [
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosectionlabel',
    'sphinxcontrib.bibtex',
]

bibtex_bibfiles = ['/home/yan/WORK/BIB/main.bib']


intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://docs.scipy.org/doc/numpy', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
    'sklearn': ('https://scikit-learn.org/stable', None),
    'pandas': ('http://pandas.pydata.org/pandas-docs/stable', None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# Make sure the target is unique
autosectionlabel_prefix_document = True


# Activate figure numbers
numfig = True
numfig_format = {'figure': "Fig. %s:"}
