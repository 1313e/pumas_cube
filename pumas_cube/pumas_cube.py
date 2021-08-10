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
from pumas_cube import __version__ as __version__
import pumas_cube.cgeometry_multi_cube as ccube

# All declaration
__all__ = ['calc_flux', 'export_to_txt', 'make_flux_plot', 'make_hist',
           'make_scatter', 'read_cube_HDF5', 'run_multi_cube']


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
attr_unit_dct = {
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
    'logE_max': r'GeV\,dex',
    'logE_min': r'GeV\,dex',
    'random_seed': ''}

dset_unit_dct = {
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


# %% HELPER FUNCTIONS
# This function calculates the fraction of the muon flux for a given charge
def _charge_fraction(charge):
    charge_ratio = 1.2766
    if(charge < 0):
        return(1/(1+charge_ratio))
    elif(charge > 0):
        return(charge_ratio/(1+charge_ratio))
    else:           # pragma: no cover
        return(1)


# This function provides the Gaisser's flux model
def _flux_gaisser(cos_theta, kinetic_energy, charge):
    Emu = kinetic_energy+0.10566
    ec = 1.1*Emu*cos_theta
    rpi = 1+ec/115
    rK = 1+ec/850
    return(1.4e03*pow(Emu, -2.7)*(1/rpi+0.054/rK)*_charge_fraction(charge))


# This function calculates Volkova's parametrization of the cosine of theta
def _cos_theta_star(cos_theta):
    p = [0.102573, -0.068287, 0.958633, 0.0407253, 0.817285]
    cs2 = (cos_theta*cos_theta+p[0]*p[0]+p[1]*pow(cos_theta, p[2])+p[3]*pow(
        cos_theta, p[4]))/(1+p[0]*p[0]+p[1]+p[3])
    return(np.sqrt(cs2) if cs2 > 0 else 0)


# This function provides the Guan et al. parametrization of the sea level flux
def _flux_gccly(cos_theta, kinetic_energy, charge):
    Emu = kinetic_energy+0.10566
    cs = _cos_theta_star(cos_theta)
    return(pow(1+3.64/(Emu*pow(cs, 1.29)), -2.7)*_flux_gaisser(
        cs, kinetic_energy, charge))


# %% FUNCTION DEFINITIONS
# This function calculates the flux for a given portion of a simulation
def calc_flux(output_dir, az_rng, el_rng, logE_rng):
    """
    Calculates the average flux with associated error for the given portion of
    the simulation in `output_dir`.

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

    Returns
    -------
    avg_flux, avg_flux_err : float
        The calculated average flux and its associated error in units of
        :math:`GeV^{-1}m^{-2}s^{-2}sr^{-1}`.

    """

    # Check if MPI worker and return if so
    if is_worker:       # pragma: no cover
        return

    # Check values of az_rng, el_rng and logE_rng
    if isinstance(az_rng, (int, np.integer)):
        az_rng = (az_rng, az_rng+1)
    else:
        if az_rng is None:
            az_rng = (0, 360)
    if isinstance(el_rng, (int, np.integer)):
        el_rng = (el_rng, el_rng+1)
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Initialize some variables
    w = 0
    w2 = 0
    n_muons = 0
    deg = np.pi/180

    # Calculate the solid angle and reweight factor of this simulation
    solid_angle = deg*np.fabs(np.sin(el_rng[1]*deg)-np.sin(el_rng[0]*deg)) *\
        (az_rng[1]-az_rng[0])
    ksel = 2*solid_angle*(logE_rng[1]-logE_rng[0])*np.log(10)

    # Loop over all requested elevations
    for el in range(*el_rng):
        # Read in data
        data = read_cube_HDF5('direction_zf', 'energy_f', 'charge',
                              'weight_i', 'energy_i', 'weight_f', 'event',
                              output_dir=output_dir,
                              az_rng=az_rng,
                              elevation=el,
                              logE_rng=logE_rng)

        # Remove attrs from the keys
        keys = list(data.keys())
        keys.remove('attrs')

        # Loop over all individual simulations
        for key in keys:
            # Loop over all data sets in this simulation
            data_i = data[key]
            for zf, energy, charge, weight, e0, w0, e in zip(
                    data_i['direction_zf'],
                    data_i['energy_f'],
                    data_i['charge'],
                    data_i['weight_f'],
                    data_i['energy_i'],
                    data_i['weight_i'],
                    data_i['event']):
                # Calculate the flux
                n_muons += 1
                if not(e == 1):
                    weight *= ksel*e0/w0  # Reweight to the selection window
                    fi = _flux_gccly(-zf, energy, charge)
                    wi = weight*fi
                    w += wi
                    w2 += wi*wi

    # Reweight flux and error to become the average values
    if w:
        w /= n_muons
        sigma = np.sqrt(((w2/n_muons)-w*w)/n_muons)/solid_angle
        w /= solid_angle
    else:
        sigma = 0

    # Return average flux and error
    return(w, sigma)


# Function that reads in an HDF5-file created with 'multi_geometry_cube.c'
def read_cube_HDF5(*args, output_dir, az_rng, elevation, logE_rng):
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
    args : positional arguments or None
        The names of all datasets that must be read from the HDF5-file.
        If `args` is empty, all datasets are read.
        If `args` is *None*, no datasets are read.
        See :obj:`~dset_unit_dct` for all dataset names.

    Returns
    -------
    data : dict
        Dict containing all attributes of the HDF5-file and every group/dataset
        that was requested in `args` or all datasets if `args` is empty.
        Groups themselves are stored as dicts as well.

    """

    # Check if MPI worker and return if so
    if is_worker:           # pragma: no cover
        return

    # Check if path exists
    if not path.exists(output_dir):
        # If not, raise error
        raise OSError("Provided input argument 'output_dir' does not exist!")

    # Obtain the name of the HDF5 file associated with this elevation
    filename = HDF5_file.format(output_dir, elevation)

    # Obtain absolute path to file
    filename = path.abspath(filename)

    # Check if elevation exists
    if not path.exists(filename):
        # If not, raise error
        raise ValueError("Provided input argument 'elevation' is invalid (%s)!"
                         % (elevation))

    # Check values of azimuth and logE_rng
    if isinstance(az_rng, (int, np.integer)):
        azimuth = [az_rng]
    else:
        if az_rng is None:
            az_rng = (0, 360)
        azimuth = np.arange(*az_rng)
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Convert provided logE to list
    n_logE = np.rint((logE_rng[1]-logE_rng[0])/logE_spc+1).astype(int)
    logE = np.around(np.linspace(logE_rng[0], logE_rng[1], n_logE), 1)
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

                # Check if args is None or not empty
                if(args == (None,)):
                    args = []
                elif not args:
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
def run_multi_cube(input_par, *, N=10000, az_rng=(0, 360), el_rng=(40, 90),
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

    # Check if given input_par exists
    if not path.exists(input_par):
        # If not, raise error
        raise OSError("Provided input argument 'input_par' does not exist!")

    # Read parameter file to obtain the output_dir
    data = np.genfromtxt(input_par, dtype=str, comments="%")
    data_dct = dict(data)
    output_dir = data_dct['output_dir']

    # Determine Az+logE
    az = np.arange(*az_rng)
    n_logE = np.rint((logE_rng[1]-logE_rng[0])/logE_spc+1).astype(int)
    logE = np.around(np.linspace(logE_rng[0], logE_rng[1], n_logE), 1)
    logE = list(zip(logE[:-1], logE[1:]))
    _, Az = np.meshgrid(np.zeros([len(logE)]), az)
    logE = list(chain(*repeat(logE, len(az))))
    az = Az.ravel()
    AzLogE = set(zip(az, logE))

    # Set required elevations
    el_all = np.arange(*el_rng)

    # Spread all requested elevations to the number of ranks
    if is_controller:
        # Determine the index gap of all elevations when spread of N cores
        idx = [el_all.shape[0]/size for _ in range(size-1)]

        # Determine the actual indices
        idx2 = np.cumsum(idx).astype(int)

        # Split them up over N cores
        el_split = np.split(el_all, idx2)

        # Scatter to other ranks
        el = comm.scatter(el_split)
    else:               # pragma: no cover
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

        # Open the HDF5-file to add the version of PUMAS_Cube as an attribute
        with h5py.File(filename, 'r+') as file:
            file.attrs['pumas_cube_version'] = __version__

    # Exit
    ccube.destroy_structs()


# This function creates a figure showing the end positions of all muons
def make_hist(dset, *, output_dir, az_rng, el_rng, logE_rng, savefig=None,
              nbins=100, figsize=(6.4, 4.8)):
    """
    Creates a histogram plot showing the `dset` values of all muons simulated
    with the given arguments.

    Parameters
    ----------
    dset : str
        The dataset for which a histogram must be made.
        See :obj:`~dset_unit_dct` for all valid dataset names.
    output_dir : str
        The path towards the directory where the output HDF5-files are stored.
    az_rng : int, 2-tuple of int or None
        All azimuth angles that must be read in from the file.
        If *None*, all azimuth angles available are read in.
        These values can be negative to count counter-clockwise.
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
        If not *None*, the filepath where the scatter plot must be saved to.
        Else, the plot will simply be shown.
    nbins : int. Default: 100
        The number of bins to use in the histogram.
    figsize : tuple of float. Default: (6.4, 4.8)
        The size of the figure in inches.

    """

    # Check if MPI worker and return if so
    if is_worker:           # pragma: no cover
        return

    # Set required elevations
    if isinstance(el_rng, (int, np.integer)):
        elevation = [el_rng]
    else:
        elevation = np.arange(*el_rng)

    # Check values of azimuth and logE_rng
    if isinstance(az_rng, (int, np.integer)):
        azimuth = [az_rng % 360]
    else:
        if az_rng is None:
            az_rng = (0, 360)
        azimuth = np.arange(*az_rng) % 360
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Create data_dct
    attrs = None
    data_dct = {}

    # Obtain data for every elevation requested
    for el in elevation:
        # Try to read this elevation
        try:
            dct = read_cube_HDF5(
                dset,
                output_dir=output_dir,
                az_rng=az_rng,
                elevation=el,
                logE_rng=logE_rng)
        except ValueError:
            # If this elevation does not exist, continue
            continue

        # Obtain attrs if not obtained already
        if attrs is None:
            attrs = dct['attrs']

        # Remove all attributes from this dict
        dct.pop('attrs')

        # Add to data_dct
        data_dct[el] = dct

    # Determine lengths
    N = int(attrs['n_particles'])
    N_el = len(data_dct)
    N_azE = len(data_dct[el])
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
    plt.title(
        r"$N_{\mathrm{par}} = %s, Az = [%s, %s]\degree, El = [%s, %s]\degree, "
        r"E_{\mathrm{det}} \in [10^{%s}, 10^{%s}]\,\mathrm{GeV}$"
        % (e13.f2tex(N_total, sdigits=2),
           int(azimuth[0]), int(azimuth[-1]+1),
           int(elevation[0]), int(elevation[-1]+1),
           logE_rng[0], logE_rng[1],
           ))

    # Labels
    ax.set_xlabel(r"%s [$\mathrm{%s}$]" % (dset, dset_unit_dct[dset]))
    ax.set_ylabel("Count")
    plt.tight_layout()

    # If savefig is not None, save the figure
    if savefig is not None:
        plt.savefig(savefig, dpi=100)
        plt.close(fig)

    # Else, simply show it
    else:
        plt.show()


# This function creates a figure showing the end positions of all muons
def make_scatter(*, output_dir, az_rng, el_rng, logE_rng, savefig=None,
                 figsize=(6.4, 4.8), cmap='cmr.ember'):
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
        These values can be negative to count counter-clockwise.
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
        If not *None*, the filepath where the scatter plot must be saved to.
        Else, the plot will simply be shown.
    figsize : tuple of float. Default: (6.4, 4.8)
        The size of the figure in inches.
    cmap : str or :obj:`~matplotlib.colors.Colormap` object
        The registered name of the colormap in :mod:`matplotlib.cm` or its
        corresponding :obj:`~matplotlib.colors.Colormap` object.

    """

    # Check if MPI worker and return if so
    if is_worker:           # pragma: no cover
        return

    # Obtain colormap
    cmap = plt.get_cmap(cmap)

    # Set required elevations
    if isinstance(el_rng, (int, np.integer)):
        elevation = [el_rng]
    else:
        elevation = np.arange(*el_rng)

    # Check values of azimuth and logE_rng
    if isinstance(az_rng, (int, np.integer)):
        azimuth = [az_rng % 360]
    else:
        if az_rng is None:
            az_rng = (0, 360)
        azimuth = np.arange(*az_rng) % 360
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Create data_dct
    attrs = None
    vmin = np.infty
    vmax = -np.infty
    data_dct = {}

    # Obtain data for every elevation requested
    for el in elevation:
        # Try to read in this elevation
        try:
            dct = read_cube_HDF5(
                'position_xf', 'position_yf', 'position_zf',
                output_dir=output_dir,
                az_rng=az_rng,
                elevation=el,
                logE_rng=logE_rng)
        except ValueError:
            # If this elevation does not exist, continue
            continue

        # Obtain attrs if not obtained already
        if attrs is None:
            attrs = dct['attrs']

        # Remove all attributes from this dict
        dct.pop('attrs')

        # Determine vmin and vmax
        vmin = min(vmin, *[min(dct[key]['position_zf']) for key in dct.keys()])
        vmax = max(vmax, *[min(dct[key]['position_zf']) for key in dct.keys()])

        # Add to data_dct
        data_dct[el] = dct

    # Determine lengths
    N = attrs['n_particles']
    N_el = len(data_dct)
    N_azE = len(data_dct[el])
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
                       s=0.01, cmap=cmap, c=data_i['position_zf'])

    # Detector position
    ax.scatter(*[[x] for x in attrs['det_position']], marker='x', s=100,
               label="Detector")
    ax.text(*attrs['det_position'], "%s" % (tuple(attrs['det_position']),),
            fontsize=9)

    # Title
    plt.title(
        r"$N_{\mathrm{par}} = %s, Az = [%s, %s]\degree, El = [%s, %s]\degree, "
        r"E_{\mathrm{det}} \in [10^{%s}, 10^{%s}]\,\mathrm{GeV}$"
        % (e13.f2tex(N_total, sdigits=2),
           int(azimuth[0]), int(azimuth[-1]+1),
           int(elevation[0]), int(elevation[-1]+1),
           logE_rng[0], logE_rng[1],
           ))

    # Labels
    ax.set_xlabel("X-pos (MGA54)[m]")
    ax.set_ylabel("Y-pos (GDA94)[m]")
    ax.set_zlabel("Z-pos (AHD)[m]")

    # Legend
    ax.legend(loc='best')
    plt.tight_layout()

    # If savefig is not None, save the figure
    if savefig is not None:
        plt.savefig(savefig, dpi=100)
        plt.close(fig)

    # Else, simply show it
    else:
        plt.show()


# This function creates a figure showing the average flux from above the sim
def make_flux_plot(*, output_dir, az_rng, el_rng, logE_rng, savefig=None,
                   figsize=(6.4, 4.8), cmap='cmr.ocean'):
    """
    Creates a plot showing a topdown of the simulation in `output_dir` with the
    average flux of every square degree for the chosen angle and energy ranges.

    Parameters
    ----------
    output_dir : str
        The path towards the directory where the output HDF5-files are stored.
    az_rng : int, 2-tuple of int or None
        All azimuth angles that must be read in from the file.
        If *None*, all azimuth angles available are read in.
        These values can be negative to count counter-clockwise.
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
        If not *None*, the filepath where the scatter plot must be saved to.
        Else, the plot will simply be shown.
    figsize : tuple of float. Default: (6.4, 4.8)
        The size of the figure in inches.
    cmap : str or :obj:`~matplotlib.colors.Colormap` object
        The registered name of the colormap in :mod:`matplotlib.cm` or its
        corresponding :obj:`~matplotlib.colors.Colormap` object.

    """

    # Check if MPI worker and return if so
    if is_worker:       # pragma: no cover
        return

    # Obtain colormap
    cmap = plt.get_cmap(cmap)

    # Check values of az_rng, el_rng and logE_rng
    if isinstance(az_rng, (int, np.integer)):
        az_rng = np.array((az_rng, az_rng+1))
        azimuth = az_rng % 360
    else:
        if az_rng is None:
            az_rng = (0, 360)
        azimuth = np.arange(*az_rng) % 360
    if isinstance(el_rng, (int, np.integer)):
        el_rng = (el_rng, el_rng+1)
        elevation = el_rng
    else:
        elevation = np.arange(*el_rng)
    if logE_rng is None:
        logE_rng = (-3, 4)

    # Create empty array for flux data
    flux_data = np.empty([len(azimuth), len(elevation)])

    # Loop over all azimuth and elevation angles, and calculate the flux
    for i, az in enumerate(azimuth):
        for j, el in enumerate(elevation):
            flux_data[i, j] = calc_flux(
                output_dir, az, el, logE_rng)[0]

    # Convert angles to format to use in a polar plot
    az = np.arange(az_rng[0], az_rng[-1]+1)*(2*np.pi/360)
    el = np.arange(el_rng[0], el_rng[-1]+1)

    # Determine the tick labels for the elevation angles
    el_axis = np.abs(el-90)
    el_axis_spc = np.ceil(abs(el_axis[-1]-el_axis[0])/5).astype(int)
    el_axis = el_axis[::el_axis_spc]

    # Create an angle meshgrid
    Az, El = np.meshgrid(az, el_axis)

    # Create new figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection='polar')

    # Set properties of this polar plot to mimic azimuth and elevation
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_thetamin(az_rng[0])
    ax.set_thetamax(az_rng[-1])
    ax.set_rgrids(el_axis, [r"$%s\degree$" % (abs(90-x)) for x in el_axis])

    # Create colormesh of flux values
    pc = ax.pcolormesh(Az, El, flux_data.T, cmap=cmap, shading='flat')
    cbar = fig.colorbar(pc)
    cbar.set_label(r"Avg. flux [$\mathrm{GeV^{-1} m^{-2} s^{-2} sr^{-1}}$]")

    # Set figure title
    plt.title(
        r"Average flux for $E_{\mathrm{det}} \in [10^{%s}, 10^{%s}]\,"
        r"\mathrm{GeV}$" % (logE_rng[0], logE_rng[1]))
    plt.tight_layout()

    # If savefig is not None, save the figure
    if savefig is not None:
        plt.savefig(savefig, dpi=100)
        plt.close(fig)

    # Else, simply show it
    else:
        plt.show()


# This function reads in data and exports it to txt
def export_to_txt(filename, *, output_dir, az_rng, el_rng, logE_rng):
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

    # Check if MPI worker and return if so
    if is_worker:       # pragma: no cover
        return

    # Obtain the absolute path to the provided filename
    filename = path.abspath(filename)

    # Set required elevations
    if isinstance(el_rng, (int, np.integer)):
        elevation = [el_rng]
    else:
        elevation = np.arange(*el_rng)

    # Create data_dct
    N = None
    data_dct = {}

    # Obtain data for every elevation requested
    for el in elevation:
        # Try to read in this elevation
        try:
            data_dct[el] = read_cube_HDF5(
                *dsets_export,
                output_dir=output_dir,
                az_rng=az_rng,
                elevation=el,
                logE_rng=logE_rng)
        except ValueError:
            # If this elevation does not exist, continue
            continue

        # Obtain N if not obtained already
        if N is None:
            N = int(data_dct[el]['attrs']['n_particles'])

        # Remove all attributes from this dict
        data_dct[el].pop('attrs')

    # Determine lengths
    N_el = len(data_dct)
    N_azE = len(data_dct[el])
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
