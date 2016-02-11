#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import os
import shutil
import sys

import shlex

sys.path.insert(0, os.path.abspath('.'))
# Use metadata from setup.py for author, version, etc.
from setup import config

# -- Generate stubs for automodule documentation --------------------------

def generate_stubs():
    """Generate an autodoc stub for each module in ./src"""
    template = '{m}\n{u}\n\n.. automodule:: src.{m}\n   :members:\n'
    shutil.rmtree('docs/src/', ignore_errors=True)
    os.mkdir('docs/src/')
    for f in os.listdir('src'):
        fname, ext = os.path.splitext(f)
        if ext == '.py' and fname != '__init__':
            with open('docs/src/' + fname + '.rst', 'w') as f:
                f.write(template.format(m=fname, u='#'*len(fname)))

generate_stubs()

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
# napoleon included from 1.3; edit if using earlier version with the addon
needs_sphinx = '1.3'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

# Automatic class documentation should also include the __init__ docstring
autoclass_content = 'both'

# Inline return types look nicer.
napoleon_use_rtype = False

# Add any paths that contain templates here, relative to this directory.
templates_path = []

# The suffix(es) of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = '3D-tools'
author = config['author']
copyright = '{}, {}'.format(datetime.date.today().year, author)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
release = config['version']
version = '{ver[0]}.{ver[1]}'.format(ver=release.split('.'))

# strftime format for |today| and update tags in html
today_fmt = html_last_updated_fmt = '%Y-%m-%d'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'README*']

# The reST default role (used for this markup: `text`) for all documents.
default_role = 'any'

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'alabaster'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    #'logo': 'logo.png',
    'github_user': 'borevitzlab',
    'github_repo': '3D-tools',
    'github_button': False,
    'travis_button': False,
}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# -- Options for LaTeX output ---------------------------------------------
latex_documents = [(master_doc, '3D-tools.tex', '3D-tools Documentation',
                    'Zac Hatfield Dodds', 'manual')]
