#!/usr/bin/env python

import os
import sys
from setuptools import find_packages

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

REQUIRES = ['numpy>=1.13', 'astropy>=2.0', 'matplotlib']
PACKAGE_PATH = os.path.abspath(os.path.join(__file__, os.pardir))

setup(
    name='pytelesto',
    version='0.0.1',
    description='Telesto reduction pipeline',
    long_description=read('README.md'),
    author='LASTRO - EPFL',
    author_email='lastro.epfl@gmail.com',
    url='https://github.com/LASTRO-EPFL/telesto-pipeline',
    download_url='https://github.com/LASTRO-EPFL/telesto-pipeline.git',
    packages=find_packages(PACKAGE_PATH),
    install_requires=REQUIRES,
    license='MIT',
    keywords='telesto',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
    ],
)
