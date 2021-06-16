Rubik's Cube Model Scripts
==========================
This repo contains the scripts that I wrote during my ADACS internship in 2021 for usage with my `Rubik's cube model`_ library.
It provides the C and accompanied Python scripts used for the simulating muons in the model of the Stawell gold mine in Victoria, Australia using the `PUMAS muon transport engine`_.

.. _Rubik's cube model: https://github.com/1313e/rubiks-cube-model
.. _PUMAS muon transport engine: https://github.com/niess/pumas

Installation
------------
Before installation, one has to make sure that both the `Rubik's cube model`_ and `PUMAS muon transport engine`_ libraries are correctly installed.
The paths towards these libraries must be provided in the ``setup.cfg`` file.
Additionally, this library requires the HDF5 library as well.

After that, this library can be easily built with the ``Makefile`` by simply using ``make`` in the root directory.
This will compile both the C scripts and the Cython extensions used by Python.
Alternatively, one can solely install the Python interface with ``pip install .`` in the root directory.

Usage
-----
After installing the package, it can be imported with ``import pumas_cube`` in any Python script.
The double Rubik's cube model together with PUMAS can then be executing by providing the ``run_double_cube`` function with a valid input parameters file (which can be found in `input <./input/input.par>`_).

