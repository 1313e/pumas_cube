# -*- coding: utf-8 -*-
# distutils: include_dirs = .
# cython: language_level=3, c_string_type=unicode, c_string_encoding=utf8

cimport geometry_double_cube as cdouble_cube


def init_structs():
    """
    Initialize all structs required for running the double Rubik's cube model.

    This function MUST be called before any other function can be used.

    """

    cdouble_cube.init_structs()


def run_double_cube(int N, double azimuth, double elevation, double logE_min,
                    double logE_max, int verbose):
    """
    Run the double Rubik's cube with the provided parameters.

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
    cdouble_cube.run_double_cube(N, azimuth, elevation, logE_min, logE_max,
                                 verbose)


def destroy_structs():
    """
    Destroy all structs initialized with :func:`~init_structs` and remove them
    from memory.

    This function must be called in order to prevent memory leaks.

    """

    cdouble_cube.destroy_structs()
