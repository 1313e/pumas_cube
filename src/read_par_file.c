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
    // Path to inner Rubik's Cube model data
    char inner_model_filename[MAX_PAR_LEN];
    // Path to outer Rubik's Cube model data
    char outer_model_filename[MAX_PAR_LEN];
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

    // INNER_MODEL_FILENAME
    strncpy(par_keys[n_par], "inner_model_filename", MAX_PAR_LEN);
    par_ptrs[n_par] = params->inner_model_filename;
    par_types[n_par++] = STRING;

    // OUTER_MODEL_FILENAME
    strncpy(par_keys[n_par], "outer_model_filename", MAX_PAR_LEN);
    par_ptrs[n_par] = params->outer_model_filename;
    par_types[n_par++] = STRING;

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
        int i, j = -1, k = -1;
        for (i=0; i<n_par; i++) {
            // Check if this line's key matches any of the normal keys
            if (strncasecmp(key_buf, par_keys[i], MAX_PAR_LEN) == 0) {
                // If the value of par_keys matches the key, then this is the correct index
                j = i;
                par_keys[i][0] = 0;
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
        else if (k >= 0) {
            // Write material name to the proper location
            strncpy(params->material_names[k], value_buf, MAX_PAR_LEN-1);
        }
        else {
            // If neither applies, this value is not allowed
            fprintf(stderr, "ERROR: Parameter key %s is invalid!!!\n", key_buf);
            return(EXIT_FAILURE);
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

    // Close file
    fclose(file);

    // Return success
    return(EXIT_SUCCESS);
}

// This function destroys the created string array in a run_params struct
void run_params_destroy(struct run_params *params){
    // Free the array itself
    free(params->material_names);
}

// Undefine global variables
#undef MAX_PARS
#undef MAX_PAR_LEN
