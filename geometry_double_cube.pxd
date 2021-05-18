# -*- coding: utf-8 -*-

cdef extern from "geometry_double_cube.c":
    void init_structs()

    void run_double_cube(
        int n_times, double azimuth, double elevation, double logE_min,
        double logE_max, int verbose)

    void destroy_structs()
