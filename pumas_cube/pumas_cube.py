# -*- coding: utf-8 -*-

"""
PUMAS Cube
==========

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
import pumas_cube.cgeometry_multi_cube as ccube

# All declaration
__all__ = ['export_to_txt', 'make_hist', 'make_scatter', 'read_cube_HDF5',
           'run_multi_cube']


# %% GLOBALS
# MPI variables
comm = MPI.HYBRID_COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
is_controller = bool(not rank)
is_worker = bool(rank)

# Write template for obtaining the name of the HDF5-file
HDF5_file = "{0}/multi_El{1:02g}.hdf5"

# Set logE spacing
logE_spc = 0.1

# Set which datasets should be exported to txt
dsets_export = [
    'position_xf', 'position_yf', 'position_zf',
    'energy_i', 'energy_f']
N_dsets = len(dsets_export)+4

# Define dict of units for all datasets
unit_dct = {
    'det_position': r'm',
    'elevation': r'\degree',
    'energy_threshold': r'GeV',
    'n_particles': '',
    'solid_angle': r'\degree',
    'azimuth': r'\degree',
    'avg_flux': r'GeV^{-1}m^{-2}s^{-2}sr^{-1}',
    'avg_flux_err': r'GeV^{-1}m^{-2}s^{-2}sr^{-1}',
    'energy_max': r'GeV',
    'energy_min': r'GeV',
    'logE_max': r'GeV dex',
    'logE_min': r'GeV dex',
    'random_seed': '',
    'azimuth_i': r'\degree',
    'charge': '',
    'decayed': '',
    'density_f': r'g/cm^3',
    'direction_xf': '',
    'direction_xi': '',
    'direction_yf': '',
    'direction_yi': '',
    'direction_zf': '',
    'direction_zi': '',
    'distance_f': r'm',
    'elevation_i': r'\degree',
    'energy_f': r'GeV',
    'energy_i': r'GeV',
    'event': '',
    'flux_f': r'GeV^{-1}m^{-2}s^{-2}sr^{-1}',
    'grammage_f': r'kg/m^2',
    'position_xf': r'm',
    'position_yf': r'm',
    'position_zf': r'm',
    'rock_id_f': '',
    'time_f': r'm/c',
    'weight_f': '',
    'weight_i': ''}


# %% FUNCTION DEFINITIONS
# Function that reads in an HDF5-file created with 'multi_geometry_cube.c'
def read_cube_HDF5(output_dir, az_rng, elevation, logE_rng, *args):
    """
    Reads in the HDF5-file created with the provided variables and returns a
    dict containing all attributes of this file plus all datasets specified in
    `args` of all groups that satisfy the given `azimuth` and `elevation`.

    Parameters
    ----------
    output_dir : str
        The path towards the directory where the output HDF5-files are stored.
    az_rng : int, 2-tuple of int or None
        All azimuth angles that must be read in from the file.
        If *None*, all azimuth angles available are read in.
    elevation : int
        The elevation angle that was simulated in the requested simulation.
    logE_rng : 2-tuple of float or None
        The range of energies that must be read in from the file.
        If *None*, all energies in the range `[-3, 4]` are read in.

    Optional
    --------
    args : positional arguments
        The names of all datasets that must be read from the HDF5-file.
        If `args` is empty, all datasets are read.

    Returns
    -------
    data : dict
        Dict containing all attributes of the HDF5-file and every group/dataset
        that was requested in `args` or all datasets if `args` is empty.
        Groups themselves are stored as dicts as well.

    """

    # Obtain the name of the HDF5 file associated with this elevation
    filename = HDF5_file.format(output_dir, elevation)

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
            if attr.startswith('cube_model_'):
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


# Function that runs the multi cube model
def run_multi_cube(input_par, N=10000, az_rng=(0, 360), el_rng=(40, 90),
                    logE_rng=(-3, 4)):
    """
    Run the multi Rubik's cube model in MPI.

    This function automatically takes care of angles that have already been
    simulated, and can be restarted from any point.

    Parameters
    ----------
    input_par : str
        The path towards the input parameters file to use for all runs.

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

    # Read parameter file to obtain the output_dir
    data = np.genfromtxt(input_par, dtype=str, comments="%")
    data_dct = dict(data)
    output_dir = data_dct['output_dir']

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
    ccube.init_structs(input_par)

    # MPI barrier
    comm.Barrier()

    # Loop over all elevations
    for e in el:
        # Obtain the name of the HDF5-file
        filename = HDF5_file.format(output_dir, e)

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
            ccube.run_multi_cube(N, azloge[0], e, *azloge[1], 0)

    # Exit
    ccube.destroy_structs()


# This function creates a figure showing the end positions of all muons
def make_hist(dset, output_dir, az_rng, el_rng, logE_rng, savefig=None,
              nbins=100, figsize=(6.4, 4.8)):
    """
    Creates a histogram plot showing the `dset` values of all muons simulated
    with the given arguments.

    Parameters
    ----------
    dset : str
        The dataset for which a histogram must be made.
    output_dir : str
        The path towards the directory where the output HDF5-files are stored.
    az_rng : int, 2-tuple of int or None
        All azimuth angles that must be read in from the file.
        If *None*, all azimuth angles available are read in.
    el_rng : int, 2-tuple of int
        The range of elevation angles that must be read in from the files.
        Note that using a range of elevations can become disorganized very
        quickly.
    logE_rng : 2-tuple of float or None
        The range of energies that must be read in from the file.
        If *None*, all energies in the range `[-3, 4]` are read in.

    Optional
    --------
    savefig : str or None. Default: None
        If not *None*, the path where the scatter plot must be saved to.
        Else, the plot will simply be shown.
    nbins : int. Default: 100
        The number of bins to use in the histogram.
    figsize : tuple of float. Default: (6.4, 4.8)
        The size of the figure in inches.

    """

    # Set required elevations
    if isinstance(el_rng, int):
        el_all = [el_rng]
    else:
        el_all = np.arange(*el_rng)

    # Check values of azimuth and logE_rng
    if isinstance(az_rng, int):
        azimuth = [az_rng]
    else:
        if az_rng is None:
            az_rng = (0, 360)
        azimuth = np.arange(*az_rng)
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Create data_dct
    attrs = None
    data_dct = {}

    # Obtain data for every elevation requested
    for elevation in el_all:
        # Read in this elevation
        dct = read_cube_HDF5(
            output_dir, az_rng, elevation, logE_rng, dset)

        # Obtain attrs if not obtained already
        if attrs is None:
            attrs = dct['attrs']

        # Remove all attributes from this dict
        dct.pop('attrs')

        # Add to data_dct
        data_dct[elevation] = dct

    # Determine lengths
    N = int(attrs['n_particles'])
    N_el = len(data_dct)
    N_azE = len(data_dct[elevation])
    N_total = N*N_el*N_azE

    # Create empty array to store flattened data in
    data_flat = np.empty(N_total)

    # Loop over all dicts and flatten their data into data_flat
    for i, el_dct in enumerate(data_dct.values()):
        for j, dct in enumerate(el_dct.values()):
            # Determine indices
            idx = N*(i*N_azE+j)

            # Add this dataset to the data_flat
            data_flat[idx:idx+N] = dct[dset]

    # Create histogram plot of all muons
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()

    # Create histogram
    ax.hist(data_flat, bins=nbins, log=True)

    # Title
    fig.suptitle(
        r"$N_{\mathrm{par}} = %s, Az = [%s, %s]\degree, El = [%s, %s]\degree, "
        r"E_{\mathrm{det}} = [10^{%s}, 10^{%s}]\,\mathrm{GeV}$"
        % (e13.f2tex(N_total, sdigits=2),
           int(azimuth[0]), int(azimuth[-1]+1),
           int(el_all[0]), int(el_all[-1]+1),
           logE_rng[0], logE_rng[1],
           ))

    # Labels
    ax.set_xlabel(r"%s [$\mathrm{%s}$]" % (dset, unit_dct[dset]))
    ax.set_ylabel("Count")

    # If savefig is not None, save the figure
    if savefig is not None:
        plt.savefig(savefig, dpi=100)
        plt.close(fig)

    # Else, simply show it
    else:
        plt.show()


# This function creates a figure showing the end positions of all muons
def make_scatter(output_dir, az_rng, el_rng, logE_rng, savefig=None,
                 figsize=(6.4, 4.8)):
    """
    Creates a scatter plot showing the final positions of all muons simulated
    with the given arguments.

    Parameters
    ----------
    output_dir : str
        The path towards the directory where the output HDF5-files are stored.
    az_rng : int, 2-tuple of int or None
        All azimuth angles that must be read in from the file.
        If *None*, all azimuth angles available are read in.
    el_rng : int, 2-tuple of int
        The range of elevation angles that must be read in from the files.
        Note that using a range of elevations can become disorganized very
        quickly.
    logE_rng : 2-tuple of float or None
        The range of energies that must be read in from the file.
        If *None*, all energies in the range `[-3, 4]` are read in.

    Optional
    --------
    savefig : str or None. Default: None
        If not *None*, the path where the scatter plot must be saved to.
        Else, the plot will simply be shown.
    figsize : tuple of float. Default: (6.4, 4.8)
        The size of the figure in inches.

    """

    # Set required elevations
    if isinstance(el_rng, int):
        el_all = [el_rng]
    else:
        el_all = np.arange(*el_rng)

    # Check values of azimuth and logE_rng
    if isinstance(az_rng, int):
        azimuth = [az_rng]
    else:
        if az_rng is None:
            az_rng = (0, 360)
        azimuth = np.arange(*az_rng)
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Create data_dct
    attrs = None
    vmin = np.infty
    vmax = -np.infty
    data_dct = {}

    # Obtain data for every elevation requested
    for elevation in el_all:
        # Read in this elevation
        dct = read_cube_HDF5(
            output_dir, az_rng, elevation, logE_rng,
            'position_xf', 'position_yf', 'position_zf')

        # Obtain attrs if not obtained already
        if attrs is None:
            attrs = dct['attrs']

        # Remove all attributes from this dict
        dct.pop('attrs')

        # Determine vmin and vmax
        vmin = min(vmin, *[min(dct[key]['position_zf']) for key in dct.keys()])
        vmax = max(vmax, *[min(dct[key]['position_zf']) for key in dct.keys()])

        # Add to data_dct
        data_dct[elevation] = dct

    # Determine lengths
    N = attrs['n_particles']
    N_el = len(data_dct)
    N_azE = len(data_dct[elevation])
    N_total = int(N*N_el*N_azE)

    # Create 3D scatter plot of exit locations of all muons
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection='3d')

    # Loop over all dicts in data_dct
    for data in data_dct.values():
        # Obtain all keys
        keys = list(data.keys())

        # Loop over all keys
        for key in keys:
            # Create scatter plot for this elevation
            data_i = data[key]
            ax.scatter(data_i['position_xf'], data_i['position_yf'],
                       data_i['position_zf'], vmin=vmin, vmax=vmax,
                       s=0.01, cmap=cmr.ember, c=data_i['position_zf'])

    # Detector position
    ax.scatter(*[[x] for x in attrs['det_position']], marker='x', s=100,
               label="Detector")
    ax.text(*attrs['det_position'], "%s" % (tuple(attrs['det_position']),),
            fontsize=9)

    # Title
    fig.suptitle(
        r"$N_{\mathrm{par}} = %s, Az = [%s, %s]\degree, El = [%s, %s]\degree, "
        r"E_{\mathrm{det}} = [10^{%s}, 10^{%s}]\,\mathrm{GeV}$"
        % (e13.f2tex(N_total, sdigits=2),
           int(azimuth[0]), int(azimuth[-1]+1),
           int(el_all[0]), int(el_all[-1]+1),
           logE_rng[0], logE_rng[1],
           ))

    # Labels
    ax.set_xlabel("X-pos (MGA54)[m]")
    ax.set_ylabel("Y-pos (GDA94)[m]")
    ax.set_zlabel("Z-pos (AHD)[m]")

    # Legend
    ax.legend(loc='best')

    # If savefig is not None, save the figure
    if savefig is not None:
        plt.savefig(savefig, dpi=100)
        plt.close(fig)

    # Else, simply show it
    else:
        plt.show()


# This function reads in data and exports it to txt
def export_to_txt(filename, output_dir, az_rng, el_rng, logE_rng):
    """
    Exports the data associated with the given arguments in a text file
    `filename`.

    The actual exported data is defined in :obj:`~dsets_export`.

    Parameters
    ----------
    filename : str
        The path of the file the data must be stored in.
    output_dir : str
        The path towards the directory where the output HDF5-files are stored.
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
    N = None
    data_dct = {}

    # Obtain data for every elevation requested
    for elevation in el_all:
        # Read in this elevation
        data_dct[elevation] = read_cube_HDF5(
            output_dir, az_rng, elevation, logE_rng, *dsets_export)

        # Obtain N if not obtained already
        if N is None:
            N = int(data_dct[elevation]['attrs']['n_particles'])

        # Remove all attributes from this dict
        data_dct[elevation].pop('attrs')

    # Determine lengths
    N_el = len(data_dct)
    N_azE = len(data_dct[elevation])
    N_total = N*N_el*N_azE

    # Create empty array to store flattened data in
    data_flat = np.empty([N_total, N_dsets])

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
                data_flat[idx:idx+N, k+4] = dct[dset_name]

    # Create header
    header = ("azimuth elevation logE_min logE_max {}").format(
              ' '.join(dsets_export))

    # Create format
    fmt = "%i %i %+4.1f %+4.1f {}".format(' '.join(["%#10.10g"]*(N_dsets-4)))

    # Save data
    np.savetxt(filename, data_flat, fmt, header=header)
