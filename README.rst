PUMAS Cube
==========
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
Running the models
++++++++++++++++++
After installing the package, it can be imported with ``import pumas_cube`` in any Python script.
This gives one access to the ``run_double_cube`` function, which runs the double Rubik's cube model together with PUMAS, and the various utility functions used for analyzing the data that this process produces.

Using the ``run_double_cube`` function is pretty straight-forward.
First, we need to prepare an input parameters file, which can be found in `input <./input/input.par>`_.
This file takes several parameters, including where all the output HDF5-files should be stored; what Rubik's cube model files should be used; what the position of the muon detector is; where the MDF file is with all required materials; and what ``rock_id`` in the Rubik's cube files corresponds to what material.

After making sure that the input parameters file has the correct values in it, one can provide it to the ``run_double_cube`` function together with a few parameters stating what the arrival directions and energies are of the muons we want to explore.
For example:

.. code:: python

    # Import pumas_cube
    import pumas_cube

    # Run double cube model
    pumas_cube.run_double_cube('input.par',
                               N=100,
                               az_rng=(230, 235),
                               el_rng=(10, 20),
                               logE_rng=(-1.0, 1.5))

Here, we ran the double Rubik's cube model with PUMAS using the ``input.par`` parameters file, and a series of input parameters.
Every simulation uses a resolution of one square degree and a logarithmic energy bin-size of :math:`0.1` (although I am planning on allowing for this particular value to be changeable by the user).
That means that in this particular case, we are asking for :math:`5*10*25=1,250` different simulations (:math:`5` azimuth angles; :math:`10` elevation angles; :math:`25` logarithmic energy bins).
Therefore, given that we also asked for :math:`100` muons PER simulation, a total of :math:`125,000` muons will be simulated with this function call.

The arrival direction and energy of each individually simulated muon is randomized within the boundaries it was given.
So, for example, if a particular run has :math:`az=230; el=10; logE_rng=(-1.0, -0.9)`, then the arrival direction and energy of every muon in that simulation will be randomized in the ranges :math:`az_rng=(230, 231); el_rng=(10, 11); logE_rng=(-1.0, -0.9)`.
The actual values that were used are stored in the output HDF5-files for every simulated muon.

Analyzing the results
+++++++++++++++++++++
After a simulation run has finished, we can read the results back into memory with the ``read_cube_HDF5`` function.
This function takes the path to a directory that contains the output HDF5-files of a simulation run, and several arguments that specify what specific part of the data one wants to read in.
For example, if we assume we executed the script above; stored the data in a folder called ``test``; and we are interested in the end positions of all muons for a detected elevation of :math:`[15, 16)`, we can do this with:

.. code:: python

    # Option 1: Use only positional arguments
    data = pumas_cube.read_cube_HDF5('test', None, 15, None,
                                     'position_xf', 'position_yf', 'position_zf')

    # Option 2: Use keyword arguments for required input arguments
    data = pumas_cube.read_cube_HDF5('position_xf', 'position_yf', 'position_zf',
                                      output_dir='test',
                                      az_rng=None,
                                      elevation=15,
                                      logE_rng=None)

Here, we specified that we want to read all azimuth angles and logarithmic energy bins that are present in the output HDF5-files in the directory ``test`` for an elevation of :math:`[15, 16)`.
What we want from these simulations is the end position of every muon, which is given by ``'position_xf', 'position_yf', 'position_zf'``.
Providing no arguments to what data should be returned will return all data from every valid simulation instead.

The ``data`` variable we end up with is a Python dict, that contains an entry called ``'attrs'`` (a dict with all attributes of the HDF5-file, like what models were used or what the detector position was) and a series of keys that each describe the azimuth angle and logarithmic energy bin range for a specific simulation.
That sounds very complicated, so let me give an example.
One of the entries in ``data`` that we obtained above, will be ``(230, -1.0, -0.9)``.
This means that this entry describes the simulation that was done with the parameters :math:`az=230; el=15; logE_rng=(-1.0, -0.9)`.
We know that the elevation was :math:`15` because that is what we asked for when calling the ``read_cube_HDF5`` function, whereas the other parameters are in the key.
The dict that belongs to this specific simulation then in turn contains all the datasets that was asked for, in this case ``'position_xf', 'position_yf', 'position_zf'``.