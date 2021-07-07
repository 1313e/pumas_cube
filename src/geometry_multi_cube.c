/*
Copyright (C) 2021 Ellert van der Velden
All rights reserved.

This software is free software, distributed under the BSD-3 License.
You may redistribute and/or modify it without any restrictions, as long as the conditions specified in the terms of the BSD-3 license (included) are met.
*/

/* Standard library includes */
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

// HDF5-library
#include "H5Gpublic.h"
#include "H5Tpublic.h"
#include "hdf5.h"

/* The PUMAS API */
#include "pumas.h"

// Rubik's cube
#include "cube.h"

// Header file
#include "../include/geometry_multi_cube.h"

// Utility C-files
#include "read_par_file.c"


// Define pi
#ifndef M_PI
/* Define pi, if unknown */
#define M_PI 3.14159265358979323846
#endif

/* Reference ground altitude out of the Rubik's cube models */
#define OUTER_GROUND_LEVEL 0

/* Reference density for the outer ground */
#define OUTER_GROUND_DENSITY 2.67

// Macros
#define min(A, B) (A < B ? A : B)
#define max(A, B) (A > B ? A : B)
#define ffloor(x) ((long)x-(x<(long)x))

// Macros for HDF5 writing
#define CHECK_STATUS_AND_RETURN_ON_FAIL(status, return_value, ...) \
    do {                                                           \
        if(status < 0) {                                           \
            fprintf(stderr, __VA_ARGS__);                          \
            return return_value;                                         \
        }                                                          \
    } while (0)

#define CREATE_SINGLE_ATTRIBUTE(group_id, attribute_name, attribute_value, h5_dtype) { \
    hid_t macro_dataspace_id = H5Screate(H5S_SCALAR);               \
    CHECK_STATUS_AND_RETURN_ON_FAIL(macro_dataspace_id, (int32_t) macro_dataspace_id, \
                                    "Could not create an attribute dataspace.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"' and the HDF5 datatype was '" #h5_dtype".\n"); \
    if(sizeof(attribute_value) != H5Tget_size(h5_dtype)) {    \
        fprintf(stderr,"Error: attribute " #attribute_name" the C size = %zu does not match the hdf5 datatype size=%zu\n", \
                sizeof(attribute_value), H5Tget_size(h5_dtype));    \
        return -1;                                                  \
    }                                                               \
    hid_t macro_attribute_id = H5Acreate(group_id, attribute_name, h5_dtype, macro_dataspace_id, H5P_DEFAULT, H5P_DEFAULT); \
    CHECK_STATUS_AND_RETURN_ON_FAIL(macro_attribute_id, (int32_t) macro_attribute_id, \
                                    "Could not create an attribute ID.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"' and the HDF5 datatype was '" #h5_dtype".\n"); \
    herr_t status = H5Awrite(macro_attribute_id, h5_dtype, &(attribute_value)); \
    CHECK_STATUS_AND_RETURN_ON_FAIL(status, (int32_t) status,       \
                                    "Could not write an attribute.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"' and the HDF5 datatype was '" #h5_dtype".\n"); \
    status = H5Aclose(macro_attribute_id);                          \
    CHECK_STATUS_AND_RETURN_ON_FAIL(status, (int32_t) status,       \
                                    "Could not close an attribute ID.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"' and the HDF5 datatype was '" #h5_dtype".\n"); \
    status = H5Sclose(macro_dataspace_id);                          \
    CHECK_STATUS_AND_RETURN_ON_FAIL(status, (int32_t) status,       \
                                    "Could not close an attribute dataspace.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"' and the HDF5 datatype was '" #h5_dtype".\n"); \
}

#define CREATE_1D_ARRAY_ATTRIBUTE(file_id, attribute_name, dims, buffer, h5_dtype) { \
    if(sizeof(buffer[0]) != H5Tget_size(h5_dtype)) {                \
        fprintf(stderr,"Error: For attribute " #attribute_name", the C size = %zu does not match the hdf5 datatype size=%zu\n", \
                sizeof(buffer[0]),  H5Tget_size(h5_dtype));         \
        return -1;                                                  \
    }                                                               \
    hid_t macro_dataspace_id = H5Screate_simple(1, dims, NULL);     \
    CHECK_STATUS_AND_RETURN_ON_FAIL(macro_dataspace_id, (int32_t) macro_dataspace_id, \
                                    "Could not create a dataspace for attribute " #attribute_name".\n" \
                                    "The dimensions of the dataspace was %d\n", (int32_t) dims[0]); \
    hid_t macro_dataset_id = H5Acreate(file_id, attribute_name, h5_dtype, macro_dataspace_id, H5P_DEFAULT, H5P_DEFAULT); \
    CHECK_STATUS_AND_RETURN_ON_FAIL(macro_dataset_id, (int32_t) macro_dataset_id, \
                                    "Could not create a dataset for attribute " #attribute_name".\n" \
                                    "The dimensions of the dataset was %d\nThe file id was %d\n.", \
                                    (int32_t) dims[0], (int32_t) file_id); \
    herr_t dset_status = H5Awrite(macro_dataset_id, h5_dtype, buffer); \
    CHECK_STATUS_AND_RETURN_ON_FAIL(dset_status, (int32_t) dset_status, \
                                    "Failed to write a dataset for attribute " #attribute_name".\n" \
                                    "The dimensions of the dataset was %d\nThe file ID was %d\n." \
                                    "The dataset ID was %d.", (int32_t) dims[0], (int32_t) file_id, \
                                    (int32_t) macro_dataset_id);    \
    dset_status = H5Aclose(macro_dataset_id);                       \
    CHECK_STATUS_AND_RETURN_ON_FAIL(dset_status, (int32_t) dset_status, \
                                    "Failed to close the dataset for attribute " #attribute_name".\n" \
                                    "The dimensions of the dataset was %d\nThe file ID was %d\n." \
                                    "The dataset ID was %d.", (int32_t) dims[0], (int32_t) file_id, \
                                    (int32_t) macro_dataset_id);    \
    dset_status = H5Sclose(macro_dataspace_id);                     \
    CHECK_STATUS_AND_RETURN_ON_FAIL(dset_status, (int32_t) dset_status, \
                                    "Failed to close the dataspace for attribute " #attribute_name".\n" \
                                    "The dimensions of the dataset was %d\nThe file ID was %d\n." \
                                    "The dataspace ID was %d.", (int32_t) dims[0], (int32_t) file_id, \
                                    (int32_t) macro_dataspace_id);  \
}

#define CREATE_STRING_ATTRIBUTE(group_id, attribute_name, attribute_value) { \
    hid_t macro_dataspace_id = H5Screate(H5S_SCALAR);                   \
    CHECK_STATUS_AND_RETURN_ON_FAIL(macro_dataspace_id, (int32_t) macro_dataspace_id, \
                                    "Could not create an attribute dataspace for a String.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
    hid_t atype = H5Tcopy(H5T_C_S1);                                  \
    CHECK_STATUS_AND_RETURN_ON_FAIL(atype, (int32_t) atype,             \
                                    "Could not copy an existing data type when creating a String attribute.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
    herr_t attr_status = H5Tset_size(atype, strlen(attribute_value));                 \
    CHECK_STATUS_AND_RETURN_ON_FAIL(attr_status, (int32_t) attr_status, \
                                    "Could not set the total size of a datatype when creating a String attribute.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
    attr_status = H5Tset_strpad(atype, H5T_STR_NULLTERM);               \
    CHECK_STATUS_AND_RETURN_ON_FAIL(attr_status, (int32_t) attr_status, \
                                    "Could not set the padding when creating a String attribute.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
    hid_t macro_attribute_id = H5Acreate(group_id, attribute_name, atype, macro_dataspace_id, H5P_DEFAULT, H5P_DEFAULT); \
    CHECK_STATUS_AND_RETURN_ON_FAIL(macro_attribute_id, (int32_t) macro_attribute_id, \
                                    "Could not create an attribute ID for string.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
    attr_status = H5Awrite(macro_attribute_id, atype, attribute_value); \
    CHECK_STATUS_AND_RETURN_ON_FAIL(attr_status, (int32_t) attr_status, \
                                    "Could not write an attribute.\n"   \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
    attr_status = H5Aclose(macro_attribute_id);                         \
    CHECK_STATUS_AND_RETURN_ON_FAIL(attr_status, (int32_t) attr_status,                            \
                                    "Could not close an attribute ID.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
    attr_status = H5Tclose(atype);                                      \
    CHECK_STATUS_AND_RETURN_ON_FAIL(attr_status, (int32_t) attr_status, \
                                    "Could not close atype value.\n"    \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
    attr_status = H5Sclose(macro_dataspace_id);                         \
    CHECK_STATUS_AND_RETURN_ON_FAIL(attr_status, (int32_t) attr_status, \
                                    "Could not close an attribute dataspace when creating a String attribute.\n" \
                                    "The attribute we wanted to create was '" #attribute_name"'.\n"); \
}

#define CREATE_AND_WRITE_1D_ARRAY(file_id, field_name, dims, buffer, h5_dtype) { \
    if(sizeof(buffer[0]) != H5Tget_size(h5_dtype)) {                \
        fprintf(stderr,"Error: For field " #field_name", the C size = %zu does not match the hdf5 datatype size=%zu\n", \
                sizeof(buffer[0]),  H5Tget_size(h5_dtype));         \
        return -1;                                                  \
    }                                                               \
    hid_t macro_dataspace_id = H5Screate_simple(1, dims, NULL);     \
    CHECK_STATUS_AND_RETURN_ON_FAIL(macro_dataspace_id, (int32_t) macro_dataspace_id, \
                                    "Could not create a dataspace for field " #field_name".\n" \
                                    "The dimensions of the dataspace was %d\n", (int32_t) dims[0]); \
    hid_t macro_dataset_id = H5Dcreate(file_id, field_name, h5_dtype, macro_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); \
    CHECK_STATUS_AND_RETURN_ON_FAIL(macro_dataset_id, (int32_t) macro_dataset_id, \
                                    "Could not create a dataset for field " #field_name".\n" \
                                    "The dimensions of the dataset was %d\nThe file id was %d\n.", \
                                    (int32_t) dims[0], (int32_t) file_id); \
    herr_t dset_status = H5Dwrite(macro_dataset_id, h5_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer); \
    CHECK_STATUS_AND_RETURN_ON_FAIL(dset_status, (int32_t) dset_status, \
                                    "Failed to write a dataset for field " #field_name".\n" \
                                    "The dimensions of the dataset was %d\nThe file ID was %d\n." \
                                    "The dataset ID was %d.", (int32_t) dims[0], (int32_t) file_id, \
                                    (int32_t) macro_dataset_id);    \
    dset_status = H5Dclose(macro_dataset_id);                       \
    CHECK_STATUS_AND_RETURN_ON_FAIL(dset_status, (int32_t) dset_status, \
                                    "Failed to close the dataset for field " #field_name".\n" \
                                    "The dimensions of the dataset was %d\nThe file ID was %d\n." \
                                    "The dataset ID was %d.", (int32_t) dims[0], (int32_t) file_id, \
                                    (int32_t) macro_dataset_id);    \
    dset_status = H5Sclose(macro_dataspace_id);                     \
    CHECK_STATUS_AND_RETURN_ON_FAIL(dset_status, (int32_t) dset_status, \
                                    "Failed to close the dataspace for field " #field_name".\n" \
                                    "The dimensions of the dataset was %d\nThe file ID was %d\n." \
                                    "The dataspace ID was %d.", (int32_t) dims[0], (int32_t) file_id, \
                                    (int32_t) macro_dataspace_id);  \
}

/* Handles for PUMAS Physics and simulation context */
static struct pumas_physics * physics = NULL;
static struct pumas_context * context = NULL;

static struct run_params params;
static struct rubiks_cube **models;
static struct pumas_medium *media;

// Estimated location of detector in (GDA94, MGA54, AHD)
static double det_position[3];

/* Altitude at which the primary flux is sampled */
static double PRIMARY_ALTITUDE;

// Define struct for storing initial and final states of a muon
struct pumas_states{
    // Number of particles
    size_t n_particles;
    // Lower azimuth limit
    double azimuth;
    // Lower elevation limit
    double elevation;
    // Angular bin solid angle
    double solid_angle;
    // Detector position
    double det_position[3];
    // Minimum log starting energy
    double logE_min;
    // Maximum log starting energy
    double logE_max;
    // Minimum starting energy
    double energy_min;
    // Maximum starting energy
    double energy_max;
    // Energy threshold
    double energy_threshold;
    // Average flux
    double avg_flux;
    // Average flux error
    double avg_flux_err;
    // Random seed of the Monte Carlo engine
    unsigned long random_seed;
    // Charge of particles
    int *charge;
    // Monte-Carlo initial (generation) weight
    float *weight_i;
    // Monte-Carlo final weight including backward transport
    float *weight_f;
    // Initial energy
    double *energy_i;
    // Final energy
    double *energy_f;
    // Distance travelled
    double *distance_f;
    // Grammage
    double *grammage_f;
    // Proper time
    double *time_f;
    // Final X-position
    double *position_xf;
    // Final Y-position
    double *position_yf;
    // Final Z-position
    double *position_zf;
    // Initial Azimuth angle
    double *azimuth_i;
    // Initial Elevation angle
    double *elevation_i;
    // Initial X-direction
    double *direction_xi;
    // Initial Y-direction
    double *direction_yi;
    // Initial Z-direction
    double *direction_zi;
    // Final X-direction
    double *direction_xf;
    // Final Y-direction
    double *direction_yf;
    // Final Z-direction
    double *direction_zf;
    // Final rock_id
    int *rock_id_f;
    // Final medium density
    double *density_f;
    // Reason for muon leaving simulation
    enum pumas_event *event;
    // Decayed flag
    hbool_t *decayed;
    // Flux
    double *flux_f;
};

// Function that frees memory allocated to a pumas_states struct
void destroy_pumas_states(struct pumas_states **states_ptr){
    // Return if no valid pumas_states was provided
    if (states_ptr == NULL || *states_ptr == NULL) {
        return;
    }

    // Free memory of all allocated arrays
    struct pumas_states *s = *states_ptr;
    free(s->charge);
    free(s->weight_i);
    free(s->weight_f);
    free(s->energy_i);
    free(s->energy_f);
    free(s->distance_f);
    free(s->grammage_f);
    free(s->time_f);
    free(s->position_xf);
    free(s->position_yf);
    free(s->position_zf);
    free(s->azimuth_i);
    free(s->elevation_i);
    free(s->direction_xi);
    free(s->direction_yi);
    free(s->direction_zi);
    free(s->direction_xf);
    free(s->direction_yf);
    free(s->direction_zf);
    free(s->rock_id_f);
    free(s->density_f);
    free(s->event);
    free(s->decayed);
    free(s->flux_f);

    // Free memory of pumas_states itself
    free(s);

    // Assign NULL to pointer
    *states_ptr = NULL;
}

void destroy_structs(){
    free(media);
    run_params_destroy(&params);
    pumas_context_destroy(&context);
    pumas_physics_destroy(&physics);
    for (int i=0; i<params.n_models; i++) {
        rubiks_cube_destroy(&models[i]);
    }
    free(models);
}

/* Gracefully exit to the OS */
static int exit_gracefully(int rc){
    destroy_structs();
    exit(rc);
}

/* Error handler for PUMAS with a graceful exit */
static void handle_error(
    enum pumas_return rc, pumas_function_t * caller, const char * message){
    /* Dump the error summary */
    fputs("pumas: library error. See details below\n", stderr);
    fprintf(stderr, "error: %s\n", message);

    /* Exit to the OS */
    exit_gracefully(EXIT_FAILURE);
}

// Define get_cube convenience function
enum rubiks_cube_return get_cube(struct pumas_state *state, int *rock_id_ptr,
                                 double *step_ptr, double *density_ptr){
    // Declare the direction sign and the cube status
    const int sgn = (context->mode.direction == PUMAS_MODE_FORWARD) ? 1 : -1;
    enum rubiks_cube_return status;

    // Search all models in order
    int i;
    for (i=0; i<params.n_models; i++) {
        status = rubiks_cube_find_cube(state->position[0], state->position[1],
                                       state->position[2], state->direction[0]*sgn,
                                       state->direction[1]*sgn, state->direction[2]*sgn,
                                       models[i], NULL, NULL, NULL, density_ptr, rock_id_ptr,
                                       step_ptr);

        // If cube was found in this model, return status
        if (status != RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND) {
            return(status);
        }
    }

    // If cube was not found in any of the models, search the outer layers instead
    const double z = state->position[2];
    double ztop;

    // Check which outer layer particle is at right now
    if (z < params.outer_level_altitude) {
        // Below ground level
        if (rock_id_ptr != NULL)
            *rock_id_ptr = params.outer_lower_material;
        if (density_ptr != NULL)
            *density_ptr = params.outer_lower_density;
        ztop = OUTER_GROUND_LEVEL;
    }
    else if (z < PRIMARY_ALTITUDE) {
        // Above ground level, below maximum altitude
        if (rock_id_ptr != NULL)
            *rock_id_ptr = params.outer_upper_material;
        if (density_ptr != NULL)
            *density_ptr = params.outer_upper_density;
        ztop = PRIMARY_ALTITUDE;
    }
    else {
        // Above maximum altitude
        return(RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND);
    }

    // If step was requested, calculate it
    if (step_ptr != NULL) {
        double uz = state->direction[2]*sgn;
        if (uz < 1E-03){
            uz = 1E-03;
        }
        const double step = (ztop - z) / uz;
        const double step_min = 1E-06;
        *step_ptr = (step > step_min) ? step : step_min;
    }

    // Return success
    return(RUBIKS_CUBE_RETURN_SUCCESS);
}

// Declare locals function
double get_locals(struct pumas_medium *medium, struct pumas_state *state,
                  struct pumas_locals *locals){
    // Declare variables
    double density;

    // Obtain the density of the cube the muon is in right now
    get_cube(state, NULL, NULL, &density);

    // Assign density
    locals->density = density*1000;

    // Return 0
    return(0);
}

// Declare geometry function
enum pumas_step get_medium(struct pumas_context *context, struct pumas_state *state,
                           struct pumas_medium **medium_ptr, double *step_ptr){
    // Declare variables
    double step;
    int rock_id;
    enum rubiks_cube_return status;

    // Obtain the cube in which the muon is right now
    status = get_cube(state, &rock_id, &step, NULL);

    // Assign medium if requested
    if (medium_ptr != NULL) {
        // Determine if the current cube is valid
        if (status != RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND) {
            // If so, retrieve proper media
            *medium_ptr = &media[rock_id];
        }
        else {
            // If not, set it to NULL
            *medium_ptr = NULL;
        }
    }

    // Assign stepping distance if requested
    if (step_ptr != NULL) {
        *step_ptr = step;
    }

    // Return that the stepping distance is exact
    return(PUMAS_STEP_RAW);
}

// This function initializes the PUMAS physics and rubiks_cube cube
void init_geometry(const char *input_par){
    // Read input parameter file
    read_par_file(input_par, &params);

    // Load in physics
    pumas_physics_create(&physics, PUMAS_PARTICLE_MUON, params.MDF_filename,
                         params.DEDX_dir, NULL);

    // Create output directory
    mkdir(params.output_dir, 0755);

    // Allocate memory for models
    models = (struct rubiks_cube **)malloc(sizeof(struct rubiks_cube *)*params.n_models);
    if (models == NULL) {
        perror("Memory cannot be allocated for Rubik's cube models!!!\n");
        exit_gracefully(EXIT_FAILURE);
    }

    // Read in models
    int i;
    for (i=0; i<params.n_models; i++) {
        rubiks_cube_create(&models[i], params.model_filenames[i]);
    }

    // Obtain detector position
    det_position[0] = params.det_position[0];
    det_position[1] = params.det_position[1];
    det_position[2] = params.det_position[2];

    // Allocate memory for media
    media = (struct pumas_medium *)malloc(sizeof(struct pumas_medium)*params.n_materials);
    if (media == NULL) {
        perror("Memory cannot be allocated for geometry media!!!\n");
        exit_gracefully(EXIT_FAILURE);
    }

    // Map PUMAS materials indices
    for (i=0; i<params.n_materials; i++) {
        pumas_physics_material_index(physics, params.material_names[i], &media[i].material);
    }

    // Set all locals to get_locals
    for (int i=0; i<params.n_materials; i++) {
        media[i].locals = &get_locals;
    }

    // Determine primary altitude
    rubiks_cube_primary_altitude(models[params.n_models-1], &PRIMARY_ALTITUDE);
}

// Function that initializes the pumas_states struct
void init_pumas_states(struct pumas_states **states_ptr, size_t n_particles){
    // Allocate memory for the pumas_states struct
    struct pumas_states *s = (struct pumas_states *)malloc(sizeof(struct pumas_states));

    // Check if memory allocation was done correctly
    if (s == NULL){
        // If not, raise error and exit
        perror("Memory cannot be allocated for pumas_states struct!!!\n");
        exit_gracefully(EXIT_FAILURE);
    }

    /* Allocate memory for all member arrays */
    s->charge = (int *)malloc(sizeof(int)*n_particles);
    s->weight_i = (float *)malloc(sizeof(float)*n_particles);
    s->weight_f = (float *)malloc(sizeof(float)*n_particles);
    s->energy_i = (double *)malloc(sizeof(double)*n_particles);
    s->energy_f = (double *)malloc(sizeof(double)*n_particles);
    s->distance_f = (double *)malloc(sizeof(double)*n_particles);
    s->grammage_f = (double *)malloc(sizeof(double)*n_particles);
    s->time_f = (double *)malloc(sizeof(double)*n_particles);
    s->position_xf = (double *)malloc(sizeof(double)*n_particles);
    s->position_yf = (double *)malloc(sizeof(double)*n_particles);
    s->position_zf = (double *)malloc(sizeof(double)*n_particles);
    s->azimuth_i = (double *)malloc(sizeof(double)*n_particles);
    s->elevation_i = (double *)malloc(sizeof(double)*n_particles);
    s->direction_xi = (double *)malloc(sizeof(double)*n_particles);
    s->direction_yi = (double *)malloc(sizeof(double)*n_particles);
    s->direction_zi = (double *)malloc(sizeof(double)*n_particles);
    s->direction_xf = (double *)malloc(sizeof(double)*n_particles);
    s->direction_yf = (double *)malloc(sizeof(double)*n_particles);
    s->direction_zf = (double *)malloc(sizeof(double)*n_particles);
    s->rock_id_f = (int *)malloc(sizeof(int)*n_particles);
    s->density_f = (double *)malloc(sizeof(double)*n_particles);
    s->event = (enum pumas_event *)malloc(sizeof(enum pumas_event)*n_particles);
    s->decayed = (hbool_t *)malloc(sizeof(hbool_t)*n_particles);
    s->flux_f = (double *)malloc(sizeof(double)*n_particles);

    // Check if memory allocation was done correctly
    if (s->charge == NULL || s->weight_i == NULL || s->weight_f == NULL ||
        s->energy_i == NULL || s->energy_f == NULL ||
        s->distance_f == NULL || s->grammage_f == NULL || s->time_f == NULL ||
        s->position_xf == NULL || s->position_yf == NULL || s->position_zf == NULL ||
        s->azimuth_i == NULL || s->elevation_i == NULL ||
        s->direction_xi == NULL || s->direction_yi == NULL || s->direction_zi == NULL ||
        s->direction_xf == NULL || s->direction_yf == NULL || s->direction_zf == NULL ||
        s->rock_id_f == NULL || s->density_f == NULL || s->event == NULL ||
        s->decayed == NULL || s->flux_f == NULL) {
        // If not, raise error and exit
        perror("Memory cannot be allocated for pumas_states struct members!!!\n");
        exit_gracefully(EXIT_FAILURE);
    }

    // Assign states to states_ptr
    *states_ptr = s;
}

// Function that writes all data in a struct pumas_states to file
static int write_states_to_file(struct pumas_states *states, const char *filename){
    // Create array for dimension size of attributes
    hsize_t dims_pos[1] = {3};

    // Open file, or create if it does not exist
    hid_t file_id;
    if (access(filename, F_OK) == 0) {
        file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else {
        file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

        // If it did not exist yet, store some additional attributes
        int i;
        for (i=0; i<params.n_models; i++) {
            char model_filename_key[80];
            sprintf(model_filename_key, "cube_model_%i", i);
            CREATE_STRING_ATTRIBUTE(file_id, model_filename_key, params.model_filenames[i]);
        }
        CREATE_SINGLE_ATTRIBUTE(file_id, "n_particles", states->n_particles, H5T_NATIVE_ULONG);
        CREATE_SINGLE_ATTRIBUTE(file_id, "elevation", states->elevation, H5T_NATIVE_DOUBLE);
        CREATE_SINGLE_ATTRIBUTE(file_id, "solid_angle", states->solid_angle, H5T_NATIVE_DOUBLE);
        CREATE_SINGLE_ATTRIBUTE(file_id, "energy_threshold", states->energy_threshold, H5T_NATIVE_DOUBLE);
        CREATE_1D_ARRAY_ATTRIBUTE(file_id, "det_position", dims_pos, states->det_position, H5T_NATIVE_DOUBLE);
    }

    // Check if azimuth group exists
    char group_name[80];
    hid_t group_id;
    sprintf(group_name, "Az%03g", states->azimuth);
    herr_t status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    status = H5Gget_objinfo(file_id, group_name, 0, NULL);

    // Open group for this azimuth, or create if it does not exist
    if (status == 0) {
        group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
    }
    else {
        group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // If it did not exist yet, write azimuth as an attribute
        CREATE_SINGLE_ATTRIBUTE(group_id, "azimuth", states->azimuth, H5T_NATIVE_DOUBLE);
    }

    // Create group for this particular run
    char subgroup_name[80];
    sprintf(subgroup_name, "logE%+04.1f_%+04.1f", states->logE_min, states->logE_max);
    hid_t subgroup_id = H5Gcreate(group_id, subgroup_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write attributes
    CREATE_SINGLE_ATTRIBUTE(subgroup_id, "logE_min", states->logE_min, H5T_NATIVE_DOUBLE);
    CREATE_SINGLE_ATTRIBUTE(subgroup_id, "logE_max", states->logE_max, H5T_NATIVE_DOUBLE);
    CREATE_SINGLE_ATTRIBUTE(subgroup_id, "energy_min", states->energy_min, H5T_NATIVE_DOUBLE);
    CREATE_SINGLE_ATTRIBUTE(subgroup_id, "energy_max", states->energy_max, H5T_NATIVE_DOUBLE);
    CREATE_SINGLE_ATTRIBUTE(subgroup_id, "avg_flux", states->avg_flux, H5T_NATIVE_DOUBLE);
    CREATE_SINGLE_ATTRIBUTE(subgroup_id, "avg_flux_err", states->avg_flux_err, H5T_NATIVE_DOUBLE);
    CREATE_SINGLE_ATTRIBUTE(subgroup_id, "random_seed", states->random_seed, H5T_NATIVE_ULONG);

    // Write all arrays to file
    hsize_t dims[1] = {states->n_particles};
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "charge", dims, states->charge, H5T_NATIVE_INT);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "weight_i", dims, states->weight_i, H5T_NATIVE_FLOAT);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "weight_f", dims, states->weight_f, H5T_NATIVE_FLOAT);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "energy_i", dims, states->energy_i, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "energy_f", dims, states->energy_f, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "distance_f", dims, states->distance_f, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "grammage_f", dims, states->grammage_f, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "time_f", dims, states->time_f, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "position_xf", dims, states->position_xf, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "position_yf", dims, states->position_yf, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "position_zf", dims, states->position_zf, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "azimuth_i", dims, states->azimuth_i, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "elevation_i", dims, states->elevation_i, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "direction_xi", dims, states->direction_xi, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "direction_yi", dims, states->direction_yi, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "direction_zi", dims, states->direction_zi, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "direction_xf", dims, states->direction_xf, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "direction_yf", dims, states->direction_yf, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "direction_zf", dims, states->direction_zf, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "rock_id_f", dims, states->rock_id_f, H5T_NATIVE_INT);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "density_f", dims, states->density_f, H5T_NATIVE_DOUBLE);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "event", dims, states->event, H5T_NATIVE_INT);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "decayed", dims, states->decayed, H5T_NATIVE_HBOOL);
    CREATE_AND_WRITE_1D_ARRAY(subgroup_id, "flux_f", dims, states->flux_f, H5T_NATIVE_DOUBLE);

    // Close file
    H5Gclose(subgroup_id);
    H5Gclose(group_id);
    H5Fclose(file_id);

    // Return
    return(EXIT_SUCCESS);
}

/* Wrap the default PUMAS PRNG (a Mersenne Twister)
 *
 * By default the seed is initialized from /dev/urandom. It can be retrieved
 * with the `pumas_context_random_seed_get` function. You might want to add the
 * used random seed to the HDF5 file?
 *
 * Alternativelly one can set a specific random seed with
 * `pumas_context_random_seed_set`.
 */
static inline double uniform01(struct pumas_context * context){
    return context->random(context);
}

/* Fraction of the muon flux for a given charge */
static double charge_fraction(double charge)
{
        /* Use a constant charge ratio.
         * Ref: CMS (https://arxiv.org/abs/1005.5332)
         */
        const double charge_ratio = 1.2766;

        if (charge < 0.)
                return 1. / (1. + charge_ratio);
        else if (charge > 0.)
                return charge_ratio / (1. + charge_ratio);
        else
                return 1.;
}

/* Gaisser's flux model
 * Ref: see e.g. the ch.30 of the PDG (https://pdglive.lbl.gov)
 */
double flux_gaisser(double cos_theta, double kinetic_energy, double charge)
{
        const double Emu = kinetic_energy + 0.10566;
        const double ec = 1.1 * Emu * cos_theta;
        const double rpi = 1. + ec / 115.;
        const double rK = 1. + ec / 850.;
        return 1.4E+03 * pow(Emu, -2.7) * (1. / rpi + 0.054 / rK) *
            charge_fraction(charge);
}

/* Volkova's parameterization of cos(theta*) */
static double cos_theta_star(double cos_theta)
{
        const double p[] = { 0.102573, -0.068287, 0.958633, 0.0407253,
                0.817285 };
        const double cs2 =
            (cos_theta * cos_theta + p[0] * p[0] + p[1] * pow(cos_theta, p[2]) +
                p[3] * pow(cos_theta, p[4])) /
            (1. + p[0] * p[0] + p[1] + p[3]);
        return cs2 > 0. ? sqrt(cs2) : 0.;
}

/*
 * Guan et al. parameterization of the sea level flux of atmospheric muons
 * Ref: https://arxiv.org/abs/1509.06176
 */
double flux_gccly(double cos_theta, double kinetic_energy, double charge)
{
        const double Emu = kinetic_energy + 0.10566;
        const double cs = cos_theta_star(cos_theta);
        return pow(1. + 3.64 / (Emu * pow(cs, 1.29)), -2.7) *
            flux_gaisser(cs, kinetic_energy, charge);
}

void init_structs(const char *input_par){
    // Set the error handler callback
    pumas_error_handler_set(&handle_error);

    // Load in physics and geometry data
    init_geometry(input_par);

    // Create new simulation context
    pumas_context_create(&context, physics, 0);

    // Configure the context for a backward transport
    context->mode.direction = PUMAS_MODE_BACKWARD;

    // Set the medium callback
    context->medium = &get_medium;

    // Enable external limit on the kinetic energy
    context->event |= PUMAS_EVENT_LIMIT_ENERGY;
}

void run_multi_cube(int n_times, double azimuth, double elevation, double logE_min, double logE_max, int verbose){
    // Initialize variables for Monte Carlo
    const double deg = M_PI/180;
    const double solid_angle = deg*fabs(sin((elevation+1)*deg)-(sin(elevation*deg))); /* The angular bin solid angle */
    double energy_min = pow(10, logE_min);
    double energy_max = pow(10, logE_max);
    const double rk = log(energy_max / energy_min);
    double w = 0., w2 = 0.;
    double dist = 0;
    int i;
    struct timespec time_1, time_2;
    double energy_threshold;

    // Create variable to store all final muon states in and add all variables that apply to the whole simulation
    static struct pumas_states *states;
    init_pumas_states(&states, n_times);
    pumas_physics_table_value(physics, PUMAS_PROPERTY_KINETIC_ENERGY, 0, 0,
                              pumas_physics_table_length(physics)-1,
                              &energy_threshold);
    pumas_context_random_seed_get(context, &states->random_seed);
    states->logE_min = logE_min;
    states->logE_max = logE_max;
    states->energy_min = energy_min;
    states->energy_max = energy_max;
    states->energy_threshold = energy_threshold;
    states->azimuth = azimuth;
    states->elevation = elevation;
    states->solid_angle = solid_angle;
    states->n_particles = n_times;
    states->det_position[0] = det_position[0];
    states->det_position[1] = det_position[1];
    states->det_position[2] = det_position[2];

    // Perform Monte Carlo
    clock_gettime(CLOCK_MONOTONIC_RAW, &time_1);
    for (i = 0; i < n_times; i++) {
        /* Set the muon final state */
        double kf, wf;
        if (rk) {
            /* The final state kinetic energy is randomised over
                * a log-uniform distribution. The Monte-Carlo weight is
                * initialised according to this generating bias PDF,
                * i.e. wf = 1 / PDF(kf).
                */
            kf = energy_min*exp(rk*uniform01(context));
            wf = kf*rk;
        } else {
            /* A point estimate is computed, for a fixed final
                * state energy.
                */
            kf = energy_min;
            wf = 1;
        }

        /* Randomize the muon charge */
        const int charge = (uniform01(context) > 0.5) ? -1 : 1;
        wf *= 2; /* Update the Monte Carlo weight accordingly */

        // Calculate random azimuth and elevation
        double azimuth_i = azimuth+uniform01(context);
        double elevation_i;
        {
                /* Generate uniformly over the bin solid angle, i.e. in
                 * sin(elevation) x azimuth
                 */
                const double smin = sin(elevation*deg);
                const double smax = sin((elevation+1.)*deg);
                const double sin_el = smin+(smax - smin)*uniform01(context);
                elevation_i = asin(sin_el)/deg;

                /* Update the weight with the bin solid angle (i.e. 1 / PDF) */
                wf *= solid_angle;
        }
        double dx = cos(elevation_i*deg)*cos(azimuth_i*deg);
        double dy = cos(elevation_i*deg)*sin(azimuth_i*deg);
        double dz = sin(elevation_i*deg);

        // Create muon state
        struct pumas_state state = {
            .charge = charge,
            .energy = kf,
            .weight = wf,
            .decayed = 0,
            .direction = {-dx, -dy, -dz}
        };
        state.position[0] = det_position[0];
        state.position[1] = det_position[1];
        state.position[2] = det_position[2];

        // Add this state to states
        states->charge[i] = charge;
        states->weight_i[i] = state.weight;
        states->energy_i[i] = state.energy;
        states->flux_f[i] = 0;
        states->azimuth_i[i] = azimuth_i;
        states->elevation_i[i] = elevation_i;
        states->direction_xi[i] = state.direction[0];
        states->direction_yi[i] = state.direction[1];
        states->direction_zi[i] = state.direction[2];

        // Initialize initial/final medium and event
        enum pumas_event event;
        struct pumas_medium *medium[2];

        /* Transport the muon backwards */
        while (state.energy < energy_threshold - FLT_EPSILON) {
            if (state.energy < 1E+02 - FLT_EPSILON) {
                /* Below 100 GeV do a detailed simulation
                    * à la Geant4, including transverse transport
                    */
                context->mode.energy_loss = PUMAS_MODE_DETAILED;
                context->mode.scattering = PUMAS_MODE_FULL_SPACE;
                context->limit.energy = 1E+02;
            } else {
                /* Do a fast simulation à la MUM */
                context->mode.energy_loss = PUMAS_MODE_HYBRID;
                context->mode.scattering = PUMAS_MODE_LONGITUDINAL;
                context->limit.energy = energy_threshold;
            }
            pumas_context_transport(context, &state, &event, medium);

            /* Check if the muon has exit the simulation area */
            if (event == PUMAS_EVENT_MEDIUM) {
                if (medium[1] == NULL) {
                    if (state.position[2] >= PRIMARY_ALTITUDE - FLT_EPSILON) {
                        /* Update the integrated flux */
                        const double fi = flux_gccly(-state.direction[2],
                                                     state.energy, state.charge);
                        const double wi = state.weight*fi;
                        w += wi;
                        w2 += wi*wi;
                        states->flux_f[i] = fi;
                    }
                    break;
                }
            } else if (event != PUMAS_EVENT_LIMIT_ENERGY) {
                /* This should not happen */
                fprintf(stderr, "error: unexpected PUMAS event `%d`\n", event);
                exit_gracefully(EXIT_FAILURE);
            }
        }

        // Add this final muon state to states
        dist += state.distance;
        states->weight_f[i] = state.weight;
        states->energy_f[i] = state.energy;
        states->distance_f[i] = state.distance;
        states->grammage_f[i] = state.grammage;
        states->time_f[i] = state.time;
        states->position_xf[i] = state.position[0];
        states->position_yf[i] = state.position[1];
        states->position_zf[i] = state.position[2];
        states->direction_xf[i] = state.direction[0];
        states->direction_yf[i] = state.direction[1];
        states->direction_zf[i] = state.direction[2];
        states->event[i] = event;
        states->decayed[i] = state.decayed;

        // Obtain the rock_id and density of the cube the muon ended up in
        if (event == PUMAS_EVENT_MEDIUM) {
            // If muon reached the muon, set both to -1
            states->rock_id_f[i] = -1;
            states->density_f[i] = -1;
        }
        else {
            // Else, obtain the values from the actual final cube
            get_cube(&state, &states->rock_id_f[i], NULL, &states->density_f[i]);
        }
    }

    /* Print the average flux over the angular bin */
    w /= n_times;
    const double sigma = sqrt(((w2/n_times)-w*w)/n_times) / solid_angle;
    w /= solid_angle;
    const char * unit = rk ? "" : "GeV^{-1} ";
    if (verbose) {
        printf("Flux : %.5lE \\pm %.5lE %sm^{-2} s^{-2} sr^{-1}\n", w, sigma, unit);
    }

    // Write all the data in states to a file
    states->avg_flux = w;
    states->avg_flux_err = sigma;
    char HDF5_filename[80];
    sprintf(HDF5_filename, "%s/multi_El%02g.hdf5", params.output_dir, elevation);
    write_states_to_file(states, HDF5_filename);
    clock_gettime(CLOCK_MONOTONIC_RAW, &time_2);

    if (verbose) {
        printf("Average distance travelled per muon in meters: %lf\n", dist/n_times);
        printf("Time taken in seconds for %i muons: %lf\n", n_times,
               (time_2.tv_sec-time_1.tv_sec)+(time_2.tv_nsec-time_1.tv_nsec)/1e9);
    }
    destroy_pumas_states(&states);
}

/* The executable main entry point */
int main(int narg, char * argv[]){
    /* Check the number of arguments */
    if (narg < 5) {
        fprintf(stderr,
            "Usage: %s INPUT_PAR N_TIMES AZIMUTH ELEVATION "
            "KINETIC_ENERGY[_MIN] [KINETIC_ENERGY_MAX]\n",
            argv[0]);
        exit_gracefully(EXIT_FAILURE);
    }

    /* Parse the arguments */
    const char *input_par = argv[1];
    const int n_times = strtod(argv[2], NULL);
    const double azimuth = strtod(argv[3], NULL);
    const double elevation = strtod(argv[4], NULL);
    const double energy_min = strtod(argv[5], NULL);
    const double energy_max =
        (narg >= 7) ? strtod(argv[6], NULL) : energy_min;

    // Init structs
    init_structs(input_par);

    // Call the model
    run_multi_cube(n_times, azimuth, elevation, energy_min, energy_max, 1);

    /* Exit to the OS */
    exit_gracefully(EXIT_SUCCESS);
}
