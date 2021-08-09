PUMAS Cube
==========
*PUMAS Cube* is a Python library that I wrote during my ADACS internship in 2021 for usage with my `Rubik's cube model`_ library.
It provides the C and accompanied Python scripts used for the simulating muons in the model of the Stawell gold mine in Victoria, Australia using the `PUMAS muon transport engine`_.

.. _Rubik's cube model: https://github.com/1313e/rubiks_cube
.. _PUMAS muon transport engine: https://github.com/niess/pumas

Installation
------------
Before installation, one has to make sure that both the `Rubik's cube model`_ and the ``dev`` branch of the `PUMAS muon transport engine`_ libraries are correctly installed, in addition to the HDF5 library.
The paths towards these libraries can be provided in the ``setup.cfg`` file in the ``[libs]`` section or set to ``None`` to be auto-discovered through the ``CPATH``; ``LIBRARY_PATH``; and ``LD_LIBRARY_PATH`` (``DYLD_LIBRARY_PATH`` on MacOS-X) environment variables.
After that, this library can be easily installed with ``$ make`` in the root directory.
If one wants to also built the C-extensions for use within C, this can be done with ``$ make lib`` (make sure that the paths to all required libraries can be found by the compiler).

WARNING:
This package will not install properly on recent MacOS-X versions, like Catalina, if one sets the library paths in the ``setup.cfg`` file.
I am not entirely sure why, but for some reason, the shared file objects are not stored properly with their absolute paths on these systems, causing the package to not be able to find them.
Instead, set the paths to the libraries to ``None`` (or remove them entirely from the file) and instead make sure they can be discovered through the environment variables mentioned above.
On older MacOS-X versions and most common Linux distributions, this is not an issue.


Usage
-----
Preparing the models
++++++++++++++++++++
The *PUMAS Cube* library requires Rubik's Cube models, created with the `Rubik's cube model`_ library, in order to function.
While the `Rubik's cube model`_ library itself does not require the `PUMAS muon transport engine`_, *PUMAS Cube* does and therefore care must be taken that the models generated are compatible with PUMAS.
In particular, the positions of every cube must be given in :math:`m`; the density is given in :math:`g/cm^3`; and the ``rock_id`` is 0-indexed.
In the input parameters file, as discussed below as well, one specifies which ``rock_id`` corresponds to which material in the PUMAS MDF.

Running the models
++++++++++++++++++
After installing the package and preparing the models, we can use the library by importing it with ``import pumas_cube`` in any Python script.
This gives one access to the ``run_multi_cube`` function, which runs the multi Rubik's Cube model together with PUMAS, and the various utility functions used for analyzing the data that this process produces.

Using the ``run_multi_cube`` function is pretty straight-forward.
First, we need to prepare an input parameters file, which can be found in `input <./input/input.par>`_.
This file takes several parameters, including where all the output HDF5-files should be stored; what Rubik's Cube model files should be used; what the position of the muon detector is; where the MDF file is with all required materials; and what ``rock_id`` in the Rubik's cube files corresponds to what material.

After making sure that the input parameters file has the correct values in it, one can provide it to the ``run_multi_cube`` function together with a few parameters stating what the arrival directions and energies are of the muons we want to explore.
For example:

.. code:: python

    # Import pumas_cube
    import pumas_cube

    # Run multi cube model
    pumas_cube.run_multi_cube('input.par',
                               N=100,
                               az_rng=(230, 235),
                               el_rng=(10, 20),
                               logE_rng=(-1.0, 1.5))

Here, we ran the multi Rubik's cube model with PUMAS using the ``input.par`` parameters file, and a series of input parameters.
Every simulation uses a resolution of one square degree and a logarithmic energy bin-size of :math:`0.1` (although I am planning on allowing for this particular value to be changeable by the user).
That means that in this particular case, we are asking for :math:`5*10*25=1,250` different simulations (:math:`5` azimuth angles; :math:`10` elevation angles; :math:`25` logarithmic energy bins).
Therefore, given that we also asked for :math:`100` muons PER simulation, a total of :math:`125,000` muons will be simulated with this function call.

The arrival direction and energy of each individually simulated muon is randomized within the boundaries it was given.
So, for example, if a particular run has ``az=230; el=10; logE_rng=(-1.0, -0.9)``, then the arrival direction and energy of every muon in that simulation will be randomized in the ranges ``az_rng=(230, 231); el_rng=(10, 11); logE_rng=(-1.0, -0.9)``.
The actual values that were used are stored in the output HDF5-files for every simulated muon.

Analyzing the results
+++++++++++++++++++++
After a simulation run has finished, we can read the results back into memory with the ``read_cube_HDF5`` function.
This function takes the path to a directory that contains the output HDF5-files of a simulation run, and several arguments that specify what specific part of the data one wants to read in.
For example, if we assume we executed the script above; stored the data in a folder called ``test``; and we are interested in the end positions of all muons for a detected elevation of :math:`[15, 16)`, we can do this with:

.. code:: python

    # Read data
    data = pumas_cube.read_cube_HDF5('position_xf', 'position_yf', 'position_zf',
                                      output_dir='test',
                                      az_rng=None,
                                      elevation=15,
                                      logE_rng=None)

Here, we specified that we want to read all azimuth angles and logarithmic energy bins that are present in the output HDF5-files in the directory ``test`` for an elevation of :math:`[15, 16)`.
What we want from these simulations is the end position of every muon, which is given by ``'position_xf', 'position_yf', 'position_zf'``.
Providing no arguments to what data should be returned will return all data from every valid simulation instead.
If required, one can check the ``pumas_cube.dset_unit_dct`` dict for the names of all dataset values that this function can take.

The ``data`` variable we end up with is a Python dict, that contains an entry called ``'attrs'`` (a dict with all attributes of the HDF5-file, like what models were used or what the detector position was) and a series of keys that each describe the azimuth angle and logarithmic energy bin range for a specific simulation.
That sounds very complicated, so let me give an example.
One of the entries in ``data`` that we obtained above, will be ``(230, -1.0, -0.9)``.
This means that this entry describes the simulation that was done with the parameters ``az=230; el=15; logE_rng=(-1.0, -0.9)``.
We know that the elevation was :math:`15` because that is what we asked for when calling the ``read_cube_HDF5`` function, whereas the other parameters are in the key.
The dict that belongs to this specific simulation then in turn contains all the datasets that was asked for, in this case ``'position_xf', 'position_yf', 'position_zf'``.

Plotting the results
++++++++++++++++++++
While we can use the ``read_cube_HDF5`` function described above to analyze the results in any way we want and write our own plotting scripts, *PUMAS Cube* provides two generic plotting functions already: ``make_hist`` and ``make_scatter``.
The ``make_hist`` function can be used to create a simple histogram of a SINGLE dataset that is stored for the simulations that satisfy the specific simulation parameters.
As stated above, one can check the ``pumas_cube.dset_unit_dct`` dict for the names of all dataset values that this function can take.
For example, let's say that we want to make a histogram of the final energies of all muons in the simulation:

.. code:: python

    # Create histogram of final energies
    pumas_cube.make_hist('energy_f',
                         output_dir='test',
                         az_rng=None,
                         el_rng=(10, 20),
                         logE_rng=None,
                         savefig='hist.png')

As shown above, the requesting data to be used in this function is almost identical to the ``read_cube_HDF5`` function, except that now a range of elevations can be given.
Be warned however that providing a large range of elevations can give a figure that might be very hard to interpret, as different elevations often result in different average distances from the detector to the edge of the union of the models.

The other function, ``make_scatter``, creates a 3D scatter plot of the end positions of all simulations that satisfy the specific simulation parameters.
Its use is very similar to the ``make_hist`` function:

.. code:: python

    # Create scatter plot of final positions
    pumas_cube.make_scatter(output_dir='test',
                            az_rng=None,
                            el_rng=(10, 20),
                            logE_rng=None,
                            savefig='scatter.png')

Like with the previous plotting function, using an elevation range that is too wide might create a figure that is hard to interpret.
