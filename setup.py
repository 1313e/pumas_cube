# -*- coding: utf-8 -*-

from setuptools import Extension, setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        [Extension(
            name='cgeometry_double_cube',
            sources=["cgeometry_double_cube.pyx"],
            libraries=['hdf5', 'cube', 'pumas'],
            extra_compile_args=[
                "-O3",
                "-fPIC"])])
)
