// Standard libraries
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Define different datatypes
enum data_types {
    DOUBLE = 1,
    INT = 2,
    STRING = 3
};

// Global variables
#define MAX_PARS        50      // Maximum number of parameters
#define MAX_PAR_LEN     256     // Maximum number of characters in parameter value

// Define struct that holds all read-in parameters
struct run_params{
    // Path to directory where output must be written to
    char output_dir[MAX_PAR_LEN];
    // Number of Rubik's cube models to use
    int n_models;
    // Model filenames
    char **model_filenames;
    // Detector position
    double det_position[3];
    // MDF filename
    char MDF_filename[MAX_PAR_LEN];
    // Path to directory with energy loss tables
    char DEDX_dir[MAX_PAR_LEN];
    // Number of end materials to use from the MDF
    int n_materials;
    // Assigned material names
    char **material_names;
    // Altitude of split between outer layers
    double outer_level_altitude;
    // Material index of outer lower layer
    int outer_lower_material;
    // Density of outer lower layer
    double outer_lower_density;
    // Material index of outer upper layer
    int outer_upper_material;
    // Density of outer upper layer
    double outer_upper_density;
};

// This function reads in the parameters from file
int read_par_file(const char *filename, struct run_params *params){
    // Define parameters
    int par_types[MAX_PARS];
    char par_keys[MAX_PARS][MAX_PAR_LEN+1];
    void *par_ptrs[MAX_PARS];

    // Set final character of all parameter values to be the null terminator
    int i;
    for (i=0; i<MAX_PARS; i++) {
        par_keys[i][MAX_PAR_LEN] = '\0';
    }

    // Define all parameters
    int n_par = 0;

    // OUTPUT_DIR
    strncpy(par_keys[n_par], "output_dir", MAX_PAR_LEN);
    par_ptrs[n_par] = params->output_dir;
    par_types[n_par++] = STRING;

    // N_MODELS
    strncpy(par_keys[n_par], "n_models", MAX_PAR_LEN);
    par_ptrs[n_par] = &(params->n_models);
    par_types[n_par++] = INT;

    // DET_POSITION
    strncpy(par_keys[n_par], "det_position_x", MAX_PAR_LEN);
    par_ptrs[n_par] = &(params->det_position[0]);
    par_types[n_par++] = DOUBLE;
    strncpy(par_keys[n_par], "det_position_y", MAX_PAR_LEN);
    par_ptrs[n_par] = &(params->det_position[1]);
    par_types[n_par++] = DOUBLE;
    strncpy(par_keys[n_par], "det_position_z", MAX_PAR_LEN);
    par_ptrs[n_par] = &(params->det_position[2]);
    par_types[n_par++] = DOUBLE;

    // MDF_FILENAME
    strncpy(par_keys[n_par], "MDF_filename", MAX_PAR_LEN);
    par_ptrs[n_par] = params->MDF_filename;
    par_types[n_par++] = STRING;

    // DEDX_DIR
    strncpy(par_keys[n_par], "DEDX_dir", MAX_PAR_LEN);
    par_ptrs[n_par] = params->DEDX_dir;
    par_types[n_par++] = STRING;

    // N_MATERIALS
    strncpy(par_keys[n_par], "n_materials", MAX_PAR_LEN);
    par_ptrs[n_par] = &(params->n_materials);
    par_types[n_par++] = INT;

    // OUTER_LEVEL_ALTITUDE
    strncpy(par_keys[n_par], "outer_level_altitude", MAX_PAR_LEN);
    par_ptrs[n_par] = &(params->outer_level_altitude);
    par_types[n_par++] = DOUBLE;

    // OUTER_LOWER_MATERIAL
    char outer_lower_material_name[MAX_PAR_LEN];
    strncpy(par_keys[n_par], "outer_lower_material", MAX_PAR_LEN);
    par_ptrs[n_par] = &outer_lower_material_name;
    par_types[n_par++] = STRING;

    // OUTER_LOWER_DENSITY
    strncpy(par_keys[n_par], "outer_lower_density", MAX_PAR_LEN);
    par_ptrs[n_par] = &(params->outer_lower_density);
    par_types[n_par++] = DOUBLE;

    // OUTER_UPPER_MATERIAL
    char outer_upper_material_name[MAX_PAR_LEN];
    strncpy(par_keys[n_par], "outer_upper_material", MAX_PAR_LEN);
    par_ptrs[n_par] = &outer_upper_material_name;
    par_types[n_par++] = STRING;

    // OUTER_UPPER_DENSITY
    strncpy(par_keys[n_par], "outer_upper_density", MAX_PAR_LEN);
    par_ptrs[n_par] = &(params->outer_upper_density);
    par_types[n_par++] = DOUBLE;

    // Open parameter file
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        // If file cannot be opened, it does not exist
        fprintf(stderr, "ERROR: Parameter file %s not found!!!\n", filename);
        return(EXIT_FAILURE);
    }

    // Start reading file
    char buffer[MAX_PAR_LEN+1];
    while (fgets(&(buffer[0]), MAX_PAR_LEN, file) != NULL) {
        // Create buffers for key and value
        char key_buf[MAX_PAR_LEN+1], value_buf[MAX_PAR_LEN+1];

        // Read in two fields from this line and skip if not possible
        if(sscanf(buffer, "%1024s %1024[^\n]", key_buf, value_buf) < 2) {
            continue;
        }

        // Skip if this line is a comment
        if (key_buf[0] == '%') {
            continue;
        }

        // Check which parameter key is on this line
        int i, j = -1, k = -1, m = -1;
        for (i=0; i<n_par; i++) {
            // Check if this line's key matches any of the normal keys
            if (strncasecmp(key_buf, par_keys[i], MAX_PAR_LEN) == 0) {
                // If the value of par_keys matches the key, then this is the correct index
                j = i;
                par_keys[i][0] = 0;
                break;
            }

            // If not, check if it is a model name
            if (strncasecmp(key_buf, "cube_model_", 10) == 0) {
                m = atoi(&key_buf[10]);
                break;
            }

            // If not, check if it is a material name
            if (strncasecmp(key_buf, "material_", 9) == 0) {
                k = atoi(&key_buf[9]);
                break;
            }
        }

        // Write parameter name to proper location
        if (j >= 0) {
            // Convert to proper type and write to location
            switch (par_types[j]) {
                case DOUBLE:
                    *((double *) par_ptrs[j]) = atof(value_buf);
                    break;
                case INT:
                    *((int *) par_ptrs[j]) = atoi(value_buf);
                    break;
                case STRING:
                    strncpy(par_ptrs[j], value_buf, MAX_PAR_LEN-1);
                    break;
            }
        }
        else if (m >= 0) {
            // Write model name to the proper location
            strncpy(params->model_filenames[m], value_buf, MAX_PAR_LEN-1);
        }
        else if (k >= 0) {
            // Write material name to the proper location
            strncpy(params->material_names[k], value_buf, MAX_PAR_LEN-1);
        }
        else {
            // If none applies, this value is not allowed
            fprintf(stderr, "ERROR: Parameter key %s is invalid!!!\n", key_buf);
            return(EXIT_FAILURE);
        }

        // Check if this was the 'n_models' key
        if (strncasecmp(key_buf, "n_models", MAX_PAR_LEN) == 0) {
            // Allocate memory for model_names
            params->model_filenames = (char **)malloc(sizeof(char *)*params->n_models);

            // Loop over all model names and assign memory
            for (i=0; i<params->n_models; i++) {
                params->model_filenames[i] = (char *)malloc(MAX_PAR_LEN);
            }
        }

        // Check if this was the 'n_materials' key
        if (strncasecmp(key_buf, "n_materials", MAX_PAR_LEN) == 0) {
            // Allocate memory for material_names
            params->material_names = (char **)malloc(sizeof(char *)*params->n_materials);

            // Loop over all material names and assign memory
            for (i=0; i<params->n_materials; i++) {
                params->material_names[i] = (char *)malloc(MAX_PAR_LEN);
            }
        }
    }

    // Check which material indices were given for outer_lower_material/outer_upper_material
    for (i=0; i<params->n_materials; i++) {
        if (strncasecmp(outer_lower_material_name, params->material_names[i], MAX_PAR_LEN) == 0) {
            params->outer_lower_material = i;
        }
        if (strncasecmp(outer_upper_material_name, params->material_names[i], MAX_PAR_LEN) == 0) {
            params->outer_upper_material = i;
        }
    }

    // Close file
    fclose(file);

    // Return success
    return(EXIT_SUCCESS);
}

// This function destroys the created string array in a run_params struct
void run_params_destroy(struct run_params *params){
    free(params->material_names);
    free(params->model_filenames);
}

// Undefine global variables
#undef MAX_PARS
#undef MAX_PAR_LEN
