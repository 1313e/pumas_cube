# -*- coding: utf-8 -*-

# %% IMPORTS
# Built-in imports
from os import path

# Package imports
import numpy as np
import pytest

# PUMAS Cube imports
from pumas_cube import (
    export_to_txt, make_hist, make_scatter, read_cube_HDF5, run_multi_cube)


# %% GLOBALS
# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain paths to materials directory and input files
input_par = path.join(dirpath, 'data', 'input.par')
materials_dir = path.join(dirpath, 'data', 'materials')
MDF_file = path.join(materials_dir, 'materials.xml')
model_file = path.join(dirpath, 'data', 'test_model.txt')


# Create fixture for providing a temporary directory with input parameters
@pytest.fixture(scope='module')
def working_dir(tmpdir_factory):
    tmpdir = tmpdir_factory.mktemp('test_multi_cube')
    input_par_file = copy_input_par(tmpdir)
    return({'output_dir': tmpdir,
            'input_par': input_par_file})


# %% HELPER FUNCTIONS
# This function creates a test Rubik's cube model with 3 layers
def create_Rubiks_cube_model(filename):     # pragma: no cover
    # Define the number of unique X, Y and Z coordinates
    nx = 15
    ny = 15
    nz = 15

    # Define the boundaries between the layers
    f_layers = np.array([1/3, 2/3])

    # Determine the X, Y and Z coordinate ranges
    x_range = np.arange(0, nx*5, 5)
    y_range = np.arange(0, ny*5, 5)
    z_range = np.arange(0, -nz*5, -5)

    # Determine the actual grid
    Z, Y, X = np.meshgrid(z_range, y_range, x_range, indexing='ij')
    x_arr = X.ravel()
    y_arr = Y.ravel()
    z_arr = Z.ravel()

    # Determine the indices of the layers
    idx = np.array(np.rint(f_layers*nz)*nx*ny, dtype=int)

    # Create a rock_id array
    rock_id_arr = np.empty(nx*ny*nz, dtype=int)
    rock_id_arr[:idx[0]] = 0
    rock_id_arr[idx[0]:idx[1]] = 1
    rock_id_arr[idx[1]:] = 2

    # Create file
    with open(filename, 'w') as file:
        # Loop over all coordinates and rock_ids
        for x, y, z, rock_id in zip(x_arr, y_arr, z_arr, rock_id_arr):
            # Write line to file using a default density of 0.00
            file.write("{:6.2f} {:6.2f} {:6.2f} 0.00 {}\n".format(
                x, y, z, rock_id))


# This function copies the input parameters file to a specific location
def copy_input_par(dest_dir):
    # Read in the input parameters file
    data = np.genfromtxt(input_par, dtype=str, comments="%")
    data_dct = dict(data)

    # Obtain absolute path to provided dest_dir
    dest_dir = path.abspath(dest_dir)

    # Modify paths
    data_dct['output_dir'] = dest_dir
    data_dct['cube_model_0'] = model_file
    data_dct['MDF_filename'] = MDF_file
    data_dct['DEDX_dir'] = materials_dir

    # Create a new data_lst list
    data_lst = [' '.join([key, value]) for key, value in data_dct.items()]

    # Convert to a single string
    data_str = '\n'.join(data_lst)

    # Save to file
    dest_file = path.join(dest_dir, 'input.par')
    with open(dest_file, 'w') as file:
        file.write(data_str)

    # Return path to created file
    return(dest_file)


# %% PYTEST CLASSES AND FUNCTIONS
# Pytest class for run_multi_cube()-function
@pytest.mark.run_on_pass
class Test_run_multi_cube(object):
    # Test default way of calling this function
    def test_default(self, working_dir):
        # Run default model
        run_multi_cube(working_dir['input_par'],
                       N=100,
                       az_rng=(0, 45),
                       el_rng=(88, 90),
                       logE_rng=(2, 2.5))

    # Test using a new range of elevations that overlaps
    def test_el_rng_overlap(self, working_dir):
        # Run default model
        run_multi_cube(working_dir['input_par'],
                       N=100,
                       az_rng=(0, 45),
                       el_rng=(85, 90),
                       logE_rng=(2, 2.5))

    # Test using a new range of elevations that does not overlap
    def test_el_rng_no_overlap(self, working_dir):
        # Run default model
        run_multi_cube(working_dir['input_par'],
                       N=100,
                       az_rng=(0, 45),
                       el_rng=(81, 83),
                       logE_rng=(2, 2.5))

    # Test if providing an invalid input_par value raises an error
    def test_invalid_input_par(self):
        with pytest.raises(OSError):
            run_multi_cube('this_is_not_a_valid_file.txt',
                           N=100,
                           az_rng=(0, 45),
                           el_rng=(85, 90),
                           logE_rng=(2, 2.5))


# Pytest class for read_cube_HDF5()-function
class Test_read_cube_HDF5(object):
    # Test default way of calling this function
    def test_default(self, working_dir):
        # Obtain data
        data = read_cube_HDF5(output_dir=working_dir['output_dir'],
                              az_rng=None,
                              elevation=85,
                              logE_rng=None)

        # Check if certain expected data can be found in this dict
        assert (0, 2.0, 2.1) in data
        assert (0, 2.4, 2.5) in data
        assert (40, 2.0, 2.1) in data
        assert (40, 2.4, 2.5) in data
        assert (60, 2.0, 2.1) not in data
        assert (0, 2.0, 2.5) not in data

    # Test if specific azimuth and logE ranges can be used
    def test_ranges(self, working_dir):
        # Obtain data
        data = read_cube_HDF5(output_dir=working_dir['output_dir'],
                              az_rng=(0, 20),
                              elevation=85,
                              logE_rng=(2, 2.3))

        # Check if certain expected data can be found in this dict
        assert (0, 2.0, 2.1) in data
        assert (0, 2.4, 2.5) not in data
        assert (15, 2.0, 2.1) in data
        assert (15, 2.4, 2.5) not in data
        assert (40, 2.0, 2.1) not in data
        assert (0, 2.0, 2.5) not in data

    # Test if specific azimuth value can be used
    def test_az_value(self, working_dir):
        # Obtain data
        data = read_cube_HDF5(output_dir=working_dir['output_dir'],
                              az_rng=5,
                              elevation=85,
                              logE_rng=(2, 2.5))

        # Check if certain expected data can be found in this dict
        assert (0, 2.0, 2.1) not in data
        assert (0, 2.4, 2.5) not in data
        assert (5, 2.0, 2.1) in data
        assert (5, 2.4, 2.5) in data
        assert (40, 2.0, 2.1) not in data
        assert (0, 2.0, 2.5) not in data

    # Test if specific datasets can be requested
    def test_request_dsets(self, working_dir):
        # Obtain data
        data = read_cube_HDF5('position_xf', 'position_yf', 'position_zf',
                              output_dir=working_dir['output_dir'],
                              az_rng=5,
                              elevation=85,
                              logE_rng=(2, 2.5))

        # Check if certain expected data can be found in this dict
        assert 'position_xf' in data[(5, 2.0, 2.1)]
        assert 'energy_f' not in data[(5, 2.0, 2.1)]

    # Test if no datasets can be requested
    def test_no_dsets(self, working_dir):
        # Obtain data
        data = read_cube_HDF5(None,
                              output_dir=working_dir['output_dir'],
                              az_rng=5,
                              elevation=85,
                              logE_rng=(2, 2.5))

        # Check if certain expected data can be found in this dict
        assert list(data[(5, 2.0, 2.1)].keys()) == ['attrs']

    # Test if giving an invalid directory raises an error
    def test_invalid_dir(self):
        with pytest.raises(OSError):
            read_cube_HDF5(output_dir='invalid_dir',
                           az_rng=(0, 20),
                           elevation=85,
                           logE_rng=(2, 2.3))

    # Test if giving an invalid elevation raises an error
    def test_invalid_elevation(self, working_dir):
        with pytest.raises(ValueError):
            read_cube_HDF5(output_dir=working_dir['output_dir'],
                           az_rng=(0, 20),
                           elevation=80,
                           logE_rng=(2, 2.3))

    # Test if requested an invalid dataset raises an error
    def test_invalid_dset(self, working_dir):
        with pytest.raises(KeyError):
            read_cube_HDF5('invalid_dset_name',
                           output_dir=working_dir['output_dir'],
                           az_rng=(0, 20),
                           elevation=85,
                           logE_rng=(2, 2.3))


# Pytest class for make_hist()-function
class Test_make_hist(object):
    # Test default
    def test_default(self, working_dir):
        # Create plot
        make_hist('energy_f',
                  output_dir=working_dir['output_dir'],
                  az_rng=None,
                  el_rng=85,
                  logE_rng=None)

    # Test if specific azimuth, elevation and logE ranges can be used
    def test_ranges(self, working_dir):
        # Create plot
        make_hist('energy_f',
                  output_dir=working_dir['output_dir'],
                  az_rng=(0, 20),
                  el_rng=(80, 90),
                  logE_rng=(2, 2.3))

    # Test if specific azimuth value can be used
    def test_az_value(self, working_dir):
        # Create plot
        make_hist('energy_f',
                  output_dir=working_dir['output_dir'],
                  az_rng=5,
                  el_rng=85,
                  logE_rng=(2, 2.5))

    # Test if the figure can be saved
    def test_savefig(self, working_dir):
        # Create plot
        filename = path.join(working_dir['output_dir'], 'test.png')
        make_hist('energy_f',
                  output_dir=working_dir['output_dir'],
                  az_rng=None,
                  el_rng=85,
                  logE_rng=None,
                  savefig=filename)


# Pytest class for make_scatter()-function
class Test_make_scatter(object):
    # Test default
    def test_default(self, working_dir):
        # Create plot
        make_scatter(output_dir=working_dir['output_dir'],
                     az_rng=None,
                     el_rng=85,
                     logE_rng=None)

    # Test if specific azimuth, elevation and logE ranges can be used
    def test_ranges(self, working_dir):
        # Create plot
        make_scatter(output_dir=working_dir['output_dir'],
                     az_rng=(0, 20),
                     el_rng=(80, 90),
                     logE_rng=(2, 2.3))

    # Test if specific azimuth value can be used
    def test_az_value(self, working_dir):
        # Create plot
        make_scatter(output_dir=working_dir['output_dir'],
                     az_rng=5,
                     el_rng=85,
                     logE_rng=(2, 2.5))

    # Test if the figure can be saved
    def test_savefig(self, working_dir):
        # Create plot
        filename = path.join(working_dir['output_dir'], 'test.png')
        make_scatter(output_dir=working_dir['output_dir'],
                     az_rng=None,
                     el_rng=85,
                     logE_rng=None,
                     savefig=filename)


# Pytest class for export_to_txt()-function
class Test_export_to_txt(object):
    # Test default
    def test_default(self, working_dir):
        # Export to file
        filename = path.join(working_dir['output_dir'], 'test.txt')
        export_to_txt(filename,
                      output_dir=working_dir['output_dir'],
                      az_rng=None,
                      el_rng=(80, 90),
                      logE_rng=None)

    # Test if specific azimuth, elevation and logE ranges can be used
    def test_ranges(self, working_dir):
        # Export to file
        filename = path.join(working_dir['output_dir'], 'test.txt')
        export_to_txt(filename,
                      output_dir=working_dir['output_dir'],
                      az_rng=(0, 20),
                      el_rng=(80, 90),
                      logE_rng=(2, 2.3))

    # Test if specific azimuth value can be used
    def test_az_value(self, working_dir):
        # Export to file
        filename = path.join(working_dir['output_dir'], 'test.txt')
        export_to_txt(filename,
                      output_dir=working_dir['output_dir'],
                      az_rng=5,
                      el_rng=85,
                      logE_rng=(2, 2.5))
