# -*- coding: utf-8 -*-

"""
Setup file for the *PUMAS_Cube* package.

"""


# %% IMPORTS
# Built-in imports
from codecs import open
import re

# Package imports
from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup


# %% SETUP DEFINITION
# Get the requirements list
with open('requirements.txt', 'r') as f:
    requirements = f.read().splitlines()

# Read the __version__.py file
with open('pumas_cube/__version__.py', 'r') as f:
    vf = f.read()

# Obtain version from read-in __version__.py file
version = re.search(r"^_*version_* = ['\"]([^'\"]*)['\"]", vf, re.M).group(1)

# Setup function declaration
setup(name="pumas_cube",
      version=version,
      author="Ellert van der Velden",
      author_email='evandervelden@swin.edu.au',
      ext_modules=cythonize(
          [Extension(
              name='pumas_cube.cgeometry_double_cube',
              sources=["pumas_cube/cgeometry_double_cube.pyx"],
              libraries=['hdf5', 'cube', 'pumas'],
              extra_compile_args=["-O3", "-fPIC"])]),
      python_requires='>=3.6, <4',
      packages=find_packages(),
      package_dir={'pumas_cube': "pumas_cube"},
      include_package_data=True,
      install_requires=requirements,
      zip_safe=False,
      )
