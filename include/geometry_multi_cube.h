#ifndef geometry_multi_cube_h
#define geometry_multi_cube_h

void init_structs(const char *input_par);

void run_multi_cube(
    int n_times, double azimuth, double elevation, double logE_min,
    double logE_max, int verbose);

void destroy_structs(void);
#endif
