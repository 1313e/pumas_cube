# -*- coding: utf-8 -*-

"""
SOMETHING
=========


"""


# %% IMPORTS
# Built-in imports
from itertools import chain, repeat
from os import path

# Package imports
import cmasher as cmr
import e13tools as e13
import h5py
import matplotlib.pyplot as plt
from mpi4pyd import MPI
import numpy as np

# Cython imports
import cgeometry_double_cube as ccube


# %% GLOBALS
# MPI variables
comm = MPI.HYBRID_COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
is_controller = bool(not rank)
is_worker = bool(rank)

# Write template for obtaining the name of the HDF5-file
HDF5_file = "data/double_N{0:g}_El{1:02g}.hdf5"

# Set logE spacing
logE_spc = 0.1

# Set which datasets should be exported to txt
dsets_export = [
    'position_xf', 'position_yf', 'position_zf',
    'energy_i', 'energy_f']
N_dsets = len(dsets_export)+7


# %% FUNCTION DEFINITIONS
# Function that reads in an HDF5-file created with 'double_geometry_cube.c'
def read_cube_HDF5(N, az_rng, elevation, logE_rng, *args):
    """
    Reads in the HDF5-file created with the provided variables and returns a
    dict containing all attributes of this file plus all datasets specified in
    `args` of all groups that satisfy the given `azimuth` and `elevation`.

    Parameters
    ----------
    N : int
        The number of muons that was simulated in the requested simulation.
    az_rng : int, 2-tuple of int or None
        All azimuth angles that must be read in from the file.
        If *None*, all azimuth angles available are read in.
    elevation : int
        The elevation angle that was simulated in the requested simulation.
    logE_rng : 2-tuple of float or None
        The range of energies that must be read in from the file.
        If *None*, all energies in the range `[-3, 4]` are read in.
    args : positional arguments
        The names of all datasets that must be read from the HDF5-file.
        If `args` is empty, all datasets are read.

    Returns
    -------
    data : dict
        Dict containing all attributes of the HDF5-file and every group/dataset
        that was requested in `args` or all datasets if `args` is empty.

    """

    # Obtain the name of the HDF5 file associated with this elevation
    filename = HDF5_file.format(N, elevation)

    # Obtain absolute path to file
    filename = path.abspath(filename)

    # Check if path exists
    if not path.exists(filename):
        # If not, raise error
        raise OSError("Provided 'filename' does not exist!")

    # Check values of azimuth and logE_rng
    if isinstance(az_rng, int):
        azimuth = [az_rng]
    else:
        if az_rng is None:
            az_rng = (0, 360)
        azimuth = np.arange(*az_rng)
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Convert provided azimuth and logE to lists
    azimuth = list(azimuth)
    logE = np.around(np.linspace(logE_rng[0], logE_rng[1],
                                 int((logE_rng[1]-logE_rng[0])/logE_spc+1)), 1)
    logE = list(zip(logE[:-1], logE[1:]))

    # Create empty dict
    data = {'attrs': {}}

    # Open HDF5-file
    with h5py.File(filename, mode='r') as file:
        # Read in all the attributes
        for attr in file.attrs:
            if attr.endswith('model_data'):
                data['attrs'][attr] = file.attrs[attr].decode('utf-8')
            else:
                data['attrs'][attr] = file.attrs[attr]

        # Loop over all azimuth requested
        for az in azimuth:
            # Try to obtain group
            try:
                group = file[f"Az{az:03g}"]
            except KeyError:
                continue

            # Obtain attributes of this group
            attrs = {}
            for attr in group.attrs:
                attrs[attr] = group.attrs[attr]

            # Loop over all logE requested
            for logE_rng in logE:
                # Try to obtain subgroup
                try:
                    subgroup =\
                        group[f"logE{logE_rng[0]:+04.1f}_{logE_rng[1]:+04.1f}"]
                except KeyError:
                    continue

                # Check if args is not empty
                if not args:
                    # If it is, all datasets are required
                    args = list(subgroup.keys())

                # Create a new dict for every combination requested
                azloge_dct = {'attrs': dict(attrs)}

                # Read in all the group attributes
                for attr in subgroup.attrs:
                    azloge_dct['attrs'][attr] = subgroup.attrs[attr]

                # Read in all requested datasets
                for name in args:
                    # Try to obtain this dataset
                    dset = subgroup.get(name)

                    # Check if dset exists
                    if dset is None:
                        # If not, raise error
                        raise KeyError(
                            f"Requested dataset '{name}' does not exist! "
                            f"Valid dataset names are {list(group.keys())}")
                    else:
                        # Else, store this dataset in the dict
                        azloge_dct[name] = dset[()]

                # Add azloge_dct to data
                data[(az, *logE_rng)] = azloge_dct

    # Return data
    return(data)


# Function that runs the double cube model
def run_double_cube(N=10000, az_rng=(0, 360), el_rng=(40, 90),
                    logE_rng=(-3, 4)):
    """
    Run the double Rubik's cube model in MPI.

    This function automatically takes care of angles that have already been
    simulated, and can be restarted from any point.

    Optional
    --------
    N : int. Default: 10000
        The number of muons to simulate per simulation.
    az_rng : 2-tuple of int. Default: (0, 360)
        The lower and upper limit of all detection azimuth angles to simulate.
    el_rng : 2-tuple of int. Default: (40, 90)
        The lower and upper limit of all detection elevation angles to
        simulate.
    logE_rng : 2-tuple of float. Default: (-3, 4)
        The lower and upper limit of all detection energies to simulate in
        logarithmic form.

    """

    # Determine Az+logE
    az = np.arange(*az_rng)
    logE = np.around(np.linspace(logE_rng[0], logE_rng[1],
                                 int((logE_rng[1]-logE_rng[0])/logE_spc+1)), 1)
    logE = list(zip(logE[:-1], logE[1:]))
    _, Az = np.meshgrid(np.zeros([len(logE)]), az)
    logE = list(chain(*repeat(logE, len(az))))
    az = Az.ravel()
    AzLogE = set(zip(az, logE))

    # Set required elevations
    if isinstance(el_rng, int):
        el_all = [el_rng]
    else:
        el_all = np.arange(*el_rng)

    # Spread all requested elevations to the number of ranks
    if is_controller:
        # Determine the index gap of all elevations when spread of N cores
        idx = [el_all.shape[0]/size for _ in range(size-1)]

        # Determine the actual indices
        idx2 = np.array(np.cumsum(idx), dtype=int)

        # Split them up over N cores
        el_split = np.split(el_all, idx2)

        # Scatter to other ranks
        el = comm.scatter(el_split)
    else:
        el = comm.scatter([])

    # Initialize structs
    ccube.init_structs()

    # MPI barrier
    comm.Barrier()

    # Loop over all elevations
    for e in el:
        # Obtain the name of the HDF5-file
        filename = HDF5_file.format(N, e)

        # Obtain all groups that are already in this file if it exists
        AzLogE_known = set()
        if path.exists(filename):
            with h5py.File(filename, 'r') as file:
                for x in file.keys():
                    az = int(x[2:5])
                    AzLogE_known.update(
                        {(az, (float(y[4:8]), float(y[9:])))
                         for y in file[x].keys()})

        # Remove AzLogE_known from AzLogE
        az_logE = AzLogE.difference(AzLogE_known)

        # Loop over all assigned AzLogE and execute
        for azloge in az_logE:
            ccube.run_double_cube(N, azloge[0], e, *azloge[1], 0)

    # Exit
    ccube.destroy_structs()


# This function creates a figure showing the end positions of all muons
def make_figure(N, az_rng, elevation, logE_rng):
    """
    Creates a figure showing the final positions of all muons simulated with
    the given arguments.

    Parameters
    ----------
    N : int
        The number of muons that was simulated in the requested simulation.
    az_rng : int, 2-tuple of int or None
        All azimuth angles that must be read in from the file.
        If *None*, all azimuth angles available are read in.
    elevation : int
        The elevation angle that was simulated in the requested simulation.
    logE_rng : 2-tuple of float or None
        The range of energies that must be read in from the file.
        If *None*, all energies in the range `[-3, 4]` are read in.

    """

    # Check values of azimuth and logE_rng
    if isinstance(az_rng, int):
        azimuth = [az_rng]
    else:
        if az_rng is None:
            az_rng = (0, 360)
        azimuth = np.arange(*az_rng)
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Load in the data from the HDF5-file
    data = read_cube_HDF5(N, az_rng, elevation, logE_rng,
                          'position_xf', 'position_yf', 'position_zf')
    attrs = data['attrs']

    # Obtain all keys in data
    keys = list(data.keys())
    keys.remove('attrs')

    # Determine lowest and highest value in the Z-direction
    vmin = min([min(data[key]['position_zf']) for key in keys])
    vmax = max([max(data[key]['position_zf']) for key in keys])

    # Create 3D scatter plot of exit locations of all muons
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for key in keys:
        data_i = data[key]
        ax.scatter(data_i['position_xf'], data_i['position_yf'],
                   data_i['position_zf'], vmin=vmin, vmax=vmax,
                   s=0.01, cmap=cmr.ember, c=data_i['position_zf'])
    ax.scatter(*[[x] for x in attrs['detector_pos']], marker='x', s=100,
               label="Detector")
    ax.text(*attrs['detector_pos'], "%s" % (tuple(attrs['detector_pos']),),
            fontsize=9)

    # Title
    fig.suptitle(
        r"$N_{\mathrm{par}} = %s, Az = [%s, %s]\degree, "
        r"El = %s\degree, E_{\mathrm{det}} = [%s, %s]\,\mathrm{GeV}, "
        r"E_{\mathrm{max}} = %s\,\mathrm{GeV}$"
        % (e13.f2tex(N*(len(keys)), sdigits=2),
           int(azimuth[0]), int(azimuth[-1]+1), int(elevation),
           e13.f2tex(10**logE_rng[0], sdigits=2),
           e13.f2tex(10**logE_rng[1], sdigits=2),
           e13.f2tex(attrs['energy_threshold'], sdigits=2),
           ))

    # Labels
    ax.set_xlabel("X-coordinate (MGA54)[m]")
    ax.set_ylabel("Y-coordinate (GDA94)[m]")
    ax.set_zlabel("Z-coordinate (AHD)[m]")

    # Legend
    ax.legend(loc='best')


# This function reads in data and exports it to txt
def export_to_txt(filename, N, az_rng, el_rng, logE_rng):
    """
    Exports the data associated with the given arguments in a text file
    `filename`.

    Parameters
    ----------
    filename : str
        The path of the file the data must be stored in.
    N : int
        The number of muons that was simulated in the requested simulation.
    az_rng : int, 2-tuple of int or None
        All azimuth angles that must be read in from the file.
        If *None*, all azimuth angles available are read in.
    el_rng : int, 2-tuple of int
        The range of elevation angles that must be read in from the files.
    logE_rng : 2-tuple of float or None
        The range of energies that must be read in from the file.
        If *None*, all energies in the range `[-3, 4]` are read in.

    """

    # Obtain the absolute path to the provided filename
    filename = path.abspath(filename)

    # Set required elevations
    if isinstance(el_rng, int):
        el_all = [el_rng]
    else:
        el_all = np.arange(*el_rng)

    # Create data_dct
    data_dct = {}

    # Create empty variable for detector position
    det_pos = []

    # Obtain data for every elevation requested
    for elevation in el_all:
        # Read in this elevation
        data_dct[elevation] = read_cube_HDF5(N, az_rng, elevation, logE_rng,
                                             *dsets_export)

        # Save detector position if not done before
        if not det_pos:
            det_pos = data_dct[elevation]['attrs']['detector_pos']

        # Remove all attributes from this dict
        data_dct[elevation].pop('attrs')

    # Determine lengths
    N_el = len(data_dct)
    N_azE = len(data_dct[elevation])

    # Create empty array to store flattened data in
    data_flat = np.empty([N*N_el*N_azE, N_dsets])

    # Set detector position
    data_flat[:, 4] = det_pos[0]
    data_flat[:, 5] = det_pos[1]
    data_flat[:, 6] = det_pos[2]

    # Loop over all dicts and flatten their data into data_flat
    for i, (el, el_dct) in enumerate(data_dct.items()):
        for j, ((az, logE_min, logE_max), dct) in enumerate(el_dct.items()):
            # Determine indices
            idx = N*(i*N_azE+j)

            # Add az, el, logE_min and logE_max to the proper columns
            data_flat[idx:idx+N, 0] = az
            data_flat[idx:idx+N, 1] = el
            data_flat[idx:idx+N, 2] = logE_min
            data_flat[idx:idx+N, 3] = logE_max

            # Add remaining data to the columns
            for k, dset_name in enumerate(dsets_export):
                data_flat[idx:idx+N, k+7] = dct[dset_name]

    # Create header
    header = ("azimuth elevation logE_min logE_max "
              "position_xi position_yi position_zi {}").format(
              ' '.join(dsets_export))

    # Create format
    fmt = "%i %i %+4.1f %+4.1f {}".format(' '.join(["%#10.10g"]*(N_dsets-4)))

    # Save data
    np.savetxt(filename, data_flat, fmt, header=header)


# %% MAIN FUNCTION
if(__name__ == '__main__'):
    # All processes run the cube
    run_double_cube(N=100, el_rng=(40, 41), logE_rng=(3, 4))

    # Make a figure
#    make_figure(100, None, 80, (0, 1))

    pass
