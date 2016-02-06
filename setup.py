#!/usr/bin/env python3
"""
Metadata and setup info for installation.

See  https://packaging.python.org/en/latest/distributing.html
"""

from setuptools import setup, find_packages

with open('README.rst') as f:
    long_description = f.read()

config = dict(
    name='forestutils',
    version='0.1.0',
    description='Tools to analyse 3D scans of a forest.',
    long_description=long_description,

    author='Zac Hatfield-Dodds',
    author_email='zac.hatfield.dodds@gmail.com',
    license='GPL3+',
    url='https://github.com/borevitzlab/3D-tools',

    keywords='forest LIDAR photogrammetery environment remote-sensing',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    install_requires=[],
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example "pip install -e .[test]" installs from path "." and
    # includes the 'test' dependency group.
    extras_require={'test': ['hypothesis', 'nose', 'pylint']},
)

if __name__ == '__main__':
    setup(**config)
