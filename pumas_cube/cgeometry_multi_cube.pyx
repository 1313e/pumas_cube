# -*- coding: utf-8 -*-
# distutils: include_dirs = ../src
# cython: language_level=3, c_string_type=unicode, c_string_encoding=utf8


# %% IMPORTS
from pumas_cube cimport geometry_multi_cube as cmulti_cube


# %% DEFINITIONS
def init_structs(const char *input_par):
    """
    Initialize all structs required for running the multi Rubik's cube model.

    This function MUST be called before any other function can be used.

    Parameters
    ----------
    input_par : str
        The path towards the input parameters file to use for all runs.

    """

    cmulti_cube.init_structs(input_par)


def run_multi_cube(int N, double azimuth, double elevation, double logE_min,
                    double logE_max, int verbose):
    """
    Run the multi Rubik's cube with the provided parameters.

    Parameters
    ----------
    N : int
        The number of muons to simulate.
    azimuth : float in [0, 360)
        The lower azimuth angle to use as the detection azimuth of every muon.
        The actual detection azimuth is randomized uniformly between [0, 1]
        plus this value for every muon.
    elevation : float in [0, 90)
        The lower elevation angle to use as the detection elevation of every
        muon.
        The actual detection elevation is randomized uniformly between [0, 1]
        plus this value for every muon.
    logE_min, logE_max : float
        The minimum/maximum detection energy of a muon in GeV in log.
        The actual detection energy is randomized uniformly between these
        values for every muon.
    verbose : {0; 1}
        The verbosity of the call.

    """

    # Call model
    cmulti_cube.run_multi_cube(N, azimuth, elevation, logE_min, logE_max,
                               verbose)


def destroy_structs():
    """
    Destroy all structs initialized with :func:`~init_structs` and remove them
    from memory.

    This function must be called in order to prevent memory leaks.

    """

    cmulti_cube.destroy_structs()
