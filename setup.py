#!/usr/bin/env python3

from setuptools import setup
import os

# We store our version number in a simple text file:
version_path = os.path.join(
    os.path.dirname(__file__),
    'sugarpy', 'version.txt'
)
with open(version_path, 'r') as version_file:
    sugary_version = version_file.read().strip()

with open('requirements.txt') as req:
    requirements = req.readlines()

setup(
    name='sugarpy',
    version=sugary_version,
    packages=['sugarpy'],
    package_dir={'sugarpy': 'sugarpy'},
    description='sugarpy',
    package_data={
        'sugarpy' : [
            'version.txt',
        ]
    },
    python_requires='>=3.5.0',
    build_requires=[
        'numpy',
    ],
    install_requires=requirements,
    long_description='SugarPy facilitates the universal, discovery-driven analysis of intact glycopeptides',
    author='Stefan Schulze, Anne Oltmanns, Christian Fufezan, Julia Krägenbring, Mechthild Pohlschröder and Michael Hippler',
    author_email='sugarpy.team@gmail.com',
    url='https://github.com/SugarPy/SugarPy',
    license='Lesser GNU General Public License (LGPL)',
    platforms='Any that supports python 3.4',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: SunOS/Solaris',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
)
