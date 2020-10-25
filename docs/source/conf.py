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
import glob
sys.path.insert(0, os.path.abspath('../../example_scripts/'))
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'SugarPy'
copyright = '2020, Stefan Schulze, Anne Oltmanns, Christian Fufezan, Julia Krägenbring, Mechthild Pohlschröder, Michael Hippler'
author = 'Stefan Schulze, Anne Oltmanns, Christian Fufezan, Julia Krägenbring, Mechthild Pohlschröder, Michael Hippler'

# The version is stored in a simple txt file
version_path = os.path.join(
    os.path.dirname(__file__),
    os.pardir, os.pardir,
    'sugarpy', 'version.txt'
)
with open(version_path, 'r') as version_file:
    sugarpy_version = version_file.read().strip()

version = sugarpy_version
release = sugarpy_version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = ['.txt', '.rst']

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

#
# Example script parsing
#
print('''
    Formatting example scripts into rst files for the docs
''')
# input()
for example_script in glob.glob('../../example_scripts/*.py'):
    if os.path.exists(example_script) is False:
        continue
    basename= os.path.basename(example_script)
    print('Reading: {0}'.format(example_script))
    with open('code_inc/{0}'.format(basename.replace('.py','.inc')), 'w') as o:
        print('''.. code-block:: python\n''', file=o)
        with open( example_script ) as infile:
            for line in infile:
                print('\t{0}'.format( line.rstrip()), file=o)

