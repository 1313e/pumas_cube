# -*- coding: utf-8 -*-

"""
Setup file for the *PUMAS_Cube* package.

"""


# %% IMPORTS
# Built-in imports
from codecs import open
from os import path
import re

# Package imports
from Cython.Build import cythonize
from distutils.dist import Distribution
from setuptools import Extension, find_packages, setup


# %% SETUP DEFINITION
# Obtain library information from setup.cfg
dist = Distribution()
dist.parse_config_files()
libs_cfg = {key: val[1] for key, val in dist.get_option_dict('libs').items()}

# Obtain absolute paths to cube_dir and pumas_dir
cube_dir = path.abspath(libs_cfg['cube_dir'])
pumas_dir = path.abspath(libs_cfg['pumas_dir'])

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
              name='pumas_cube.cgeometry_multi_cube',
              sources=["pumas_cube/cgeometry_multi_cube.pyx"],
              libraries=['hdf5', 'cube', 'pumas'],
              extra_compile_args=["-O3", "-fPIC"],
              library_dirs=[
                  path.join(cube_dir, "lib"),
                  path.join(pumas_dir, "lib")],
              runtime_library_dirs=[
                  path.join(cube_dir, "lib"),
                  path.join(pumas_dir, "lib")],
              include_dirs=[
                  path.join(cube_dir, "include"),
                  path.join(pumas_dir, "include")],
              language='c'
              )],
          compiler_directives={
              'embedsignature': True}
          ),
      python_requires='>=3.6, <4',
      packages=find_packages(),
      package_dir={'pumas_cube': "pumas_cube"},
      include_package_data=True,
      install_requires=requirements,
      zip_safe=False,
      )
