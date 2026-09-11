# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

project = 'AstroEBSD_v2'
copyright = '2025, T.B Britton, A. Foden, T. McAuliffe, T. Zhang, L. Berners.'
author = 'T.B Britton, A. Foden, T. McAuliffe, T. Zhang, L. Berners.'
release = 'v2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
extensions = ['sphinxcontrib.matlab', 'sphinx.ext.autodoc']#, 'sphinx_iframes']
primary_domain = "mat"
templates_path = ['_templates']

sys.path.insert(0, os.path.abspath(os.path.join("..", "..", "..")))

matlab_src_dir = os.path.abspath(os.path.join("..", ".."))
matlab_show_property_default_value = True
matlab_show_property_specs = True



exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'agogo'
html_static_path = ['_static']
html_theme = 'sphinx_rtd_theme'