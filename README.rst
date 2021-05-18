Rubik's Cube Model Scripts
==========================
This repo contains the scripts that I wrote during my ADACS internship in 2021 for usage with my `Rubik's cube model`_ library.
It provides the C and accompanied Python scripts used for the simulating muons in the model of the Stawell gold mine in Victoria, Australia using the `PUMAS muon transport engine`_.

.. _Rubik's cube model: https://github.com/1313e/rubiks-cube-model
.. _PUMAS muon transport engine: https://github.com/niess/pumas

Installation
------------
This library can be easily built with the ``Makefile`` by simply using ``make`` in the root directory.
This will compile both the C scripts and the Cython extensions used by Python.

Usage
-----
The double Rubik's cube model in ``geometry_double_cube.c`` can be executed by either calling it directly from the command-line or using the ``run_double_cube()``-function in Python.
