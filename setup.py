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

# Obtain paths to cube_dir, pumas_dir and hdf5_dir
cube_dir = libs_cfg.get('cube_dir', 'None')
pumas_dir = libs_cfg.get('pumas_dir', 'None')
hdf5_dir = libs_cfg.get('hdf5_dir', 'None')

# Get absolute paths to these dirs
cube_dir = path.abspath(cube_dir) if cube_dir != 'None' else None
pumas_dir = path.abspath(pumas_dir) if pumas_dir != 'None' else None
hdf5_dir = path.abspath(hdf5_dir) if hdf5_dir != 'None' else None

# Create lib_dirs and include_dirs lists
lib_dirs = [path.join(d, 'lib')
            for d in [cube_dir, pumas_dir, hdf5_dir] if d is not None]
include_dirs = [path.join(d, 'include')
                for d in [cube_dir, pumas_dir, hdf5_dir] if d is not None]

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
              library_dirs=lib_dirs,
              runtime_library_dirs=lib_dirs,
              include_dirs=include_dirs,
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
