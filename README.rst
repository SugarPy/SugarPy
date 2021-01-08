SugarPy
#######

**Universal, discovery-driven analysis of intact glycopeptides**

|doc-status|

.. |doc-status| image:: https://readthedocs.org/projects/sugarpy-ms/badge/?version=latest
    :target: https://sugarpy-ms.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Summary
*******

SugarPy facilitates the glycan database independent, discovery-driven identification of intact glycopeptides from in-source collision-induced dissociation mass spectrometry data.

Abstract
********

Protein glycosylation is a complex post-translational modification with crucial cellular functions in all domains of life. Currently, large-scale glycoproteomics approaches rely on glycan database dependent algorithms and are thus unsuitable for discovery-driven analyses of glycoproteomes. Therefore, we devised SugarPy, a glycan database independent Python module, and validated it on the glycoproteome of human breast milk. We further demonstrated its applicability by analyzing glycoproteomes with uncommon glycans stemming from the green algae Chlamydomonas reinhardtii and the archaeon Haloferax volcanii. SugarPy also facilitated the novel characterization of glycoproteins from the red alga Cyanidioschyzon merolae.

Download and Installation
*************************

Get the latest version via GitHub: https://github.com/SugarPy/SugarPy

as a .zip package: https://github.com/SugarPy/SugarPy/archive/master.zip

or via git clone::

    git clone https://github.com/SugarPy/SugarPy.git

Navidate into the downloaded/cloned folder and install the requirements::

    user@localhost:~$ cd SugarPy
    user@localhost:~/SugarPy$ pip install -r requirements.txt

After that, SugarPy can be installed into the Python sitepackages using::

    user@localhost:~/SugarPy$ python setup.py install

.. note::

    You might need administrator privileges to write in the Python site-package folder.
    On Linux or OS X, use ```sudo python setup.py install``` or write into a user folder
    by using this command ```python setup.py install --user```. On Windows, you have to
    start the command line with administrator privileges.

Usage
*****

SugarPy can be used as a standalone engine and some example scripts are provided
in the example_scripts folder.
However, we recommend to use it within the Python framework `Ursgal`_ ,
which combines various proteomics tools and therefore allows to easily set up a complete workflow.
Example scripts for Ursgal workflows that include SugarPy are given in Ursgals example_scripts folder.

.. _Ursgal:
    https://github.com/ursgal/ursgal
    

Documentation
*************

For a more detailed documentation of SugarPy's structure, please refer to
the documentation folder or https://sugarpy-ms.readthedocs.io


Questions and Participation
***************************

If you encounter any problems you can open up issues at GitHub, or write an email to sugarpy.team@gmail.com.

For any contributions, fork us at https://github.com/SugarPy/SugarPy and open up pull requests!


Disclaimer
**********

While automatic filtering of SugarPy results leads to the reliable identification of glycopeptides
in the cases we have tested, we recommend manual inspection of results and to alter filter criteria
as needed.

Copyrights
**********

Copyright 2019-today by authors and contributors in alphabetical order

* Christian Fufezan
* Michael Hippler
* Julia Krägenbring
* Michael Mormann
* Anne Oltmanns
* Mecky Pohlschröder
* Stefan Schulze

Citation
********


Schulze, S.; Oltmanns, A.; Fufezan, C.; Krägenbring, J.; Mormann, M.; Pohlschröder, M.; Hippler, M. (2020). SugarPy facilitates the universal, discovery-driven analysis of intact glycopeptides. `Bioinformatics`_ 

.. _Bioinformatics:
    https://doi.org/10.1093/bioinformatics/btaa1042

