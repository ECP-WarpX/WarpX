"""Methods for storing parameters from a run.

This should all be made to Pickle easily, so that we can load the parameters at
a later time.
"""
import collections
import copy
import datetime
import logging
import os
import platform
import subprocess

import dill
import mewarpx
import numpy as np

import mewarpx.utils_store.util as mwxutil

from mewarpx.mwxrun import mwxrun
from functools import reduce

logger = logging.getLogger(__name__)

class SimInfo(object):

    """Store basic simulation parameters used throughout analysis.

    This must contain:
        SimInfo.nxyz (int nx, int ny, int nz)
        SimInfo.pos_lims (floats xmin, xmax, ymin, ymax, zmin, zmax)
        SimInfo.geom (str 'XZ', 'RZ' or 'XYZ')
        SimInfo.dt (float)
        SimInfo.periodic (bool)

    Arguments:
        nxyz (list): tuple of ints, holding the grid dimensions
        pos_lims (list) : tuple of floats, holding the position limits along each
            dimension (xmax, ymin, ymax, zmin, zmax)
        geom (str): geometry string 'XZ', 'RZ', or 'XYZ' for 2D grid,
            cylindrical, or 3D grid
        dt (float): simulation timestep interval
        periodic (bool): use periodic boundary conditions

    Note:
        This base class has been created to provide additional plotting
        functionality, and allow post-run creation of this object if needed.
    """

    def __init__(self, nxyz, pos_lims, geom, dt, periodic=True):
        self.nxyz = nxyz
        self.pos_lims = pos_lims
        self.geom = geom
        self.dt = dt
        self.periodic = periodic

    def get_vec(self, axis):
        if axis == 'r':
            return self.get_rvec()
        axis_dict = {0: 0, 1: 1, 2: 2, 'x': 0, 'y': 1, 'z': 2}
        axis = axis_dict[axis]
        npts = self.nxyz[axis]
        xmin = self.pos_lims[2*axis]
        xmax = self.pos_lims[2*axis + 1]
        # There is one more point on the grid than cell number
        return np.linspace(xmin, xmax, npts + 1)

    def get_rvec(self):
        raise NotImplementedError
        # npts = self.nxyz[0]
        # xmin = 0.
        # xmax = self.pos_lims[1]
        # There is one more point on the grid than cell number
        # return np.linspace(xmin, xmax, npts + 1)



class WarpXSimInfo(SimInfo):

    """Create a SimInfo object when running a warpx simulation."""

    def __init__(self):
        self.nxyz = [mwxrun.nx, mwxrun.ny, mwxrun.nz]
        self.pos_lims = [mwxrun.xmin, mwxrun.xmax,
                         mwxrun.ymin, mwxrun.ymax,
                         mwxrun.zmin, mwxrun.zmax]
        self.geom = mwxrun.geom_str
        self.dt = mwxrun.get_dt()
        # self.periodic = (
            #   warp.w3d.boundxy is warp.periodic
            # The rwall in RZ simulations can be set even if boundxy is
            # periodic. It defaults to 1e36 so this should be fine.
        #    and warp.top.prwall > 1e10
        # )


class RunInfo(SimInfo):

    """Store parameters for a WARP run.

    Attributes:
        diagnames_dict (dict): Dictionary of paths to diagnostic outputs.
        component_labels (collections.OrderedDict): Ordered dictionary of
        (name, pretty_label) for each device component.
    """

    diagnames_dict = {'part_diag_dir': 'diags/xzsolver/hdf5',
                      'field_base_path': 'diags/fields/'}
    component_labels = collections.OrderedDict([
        ("cathode", "Cathode"), ("anode", "Anode"),
        ("accgrid", "Acceleration Grid"), ("absgrid", "Absorption Grid"),
        ("focusgrid", "Focus Grid"), ("dielectric", "Dielectrics"),
        ("aux", "Aux Electrode"), ("cs_ionization", "Cs ionization"),
        ("ar_ionization", "Ar ionization"),
        ("ionization", "Ionization"),
    ])

    def __init__(self, injector_dict, surface_dict,
                 electrode_params, user_var_dict=None, local_vars=None,
                 **kwargs):
        """**Initialization**: Set basic run parameters

        Call just before creating final diagnostics. This uses built-in WARP
        variables, in addition to the passed quantities, to save relevant
        parameters for diagnostic calculations and postprocessing.

        Arguments:
            injector_dict (dict): Dictionary where keys are strings for groups
                of surfaces or interaction physics (see component_labels keys
                for main possibilities) and values are injectors or lists of
                injectors.  Used by various diagnostics, so should be complete.
                ``cathode`` is a required key.
            surface_dict (dict): Dictionary where keys are groups of surfaces
                which scrape particles and values are surfaces or lists of
                surfaces.  Used by various diagnostics, so should be complete.
                ``cathode`` is a required key.
            user_var_dict (dict): A dictionary of all parameters set by user.
                This will depend on the given system being investigated.
            local_vars (dict): Instead of user_var_dict, pass locals() from the
                executing script. Then all UPPER_CASE variables will be saved to
                user_var_dict internally.
            electrode_params (dict): A dictionary of electrode parameters:

                * `CATHODE_A` (Amp/m^2/K^2)
                * `CATHODE_TEMP` (K)
                * `ANODE_TEMP` (K) (if back-emission used)
                * `ANODE_A` (Amp/m^2/K^2) (if back-emission used)
                * Additional elements if desired
            run_param_dict (dict): Optional, a dictionary of run parameters such
                as diag_timesteps and total_timesteps.
            run_file (str): Path to the main simulation file, in order to store
                it in run_param_dict for easy access later.
            area (float): Optional, cross-sectional area in m^2 to use for flux
                calculations from this simulation. If not supplied, will be
                calculated automatically from the geometry.
        """
        # Parameter dictionaries and injectors
        if local_vars is not None:
            if user_var_dict is not None:
                raise ValueError(
                    "Specify local_vars or user_var_dict, not both.")
            self.user_var_dict = self._process_local_vars(local_vars)
        else:
            self.user_var_dict = user_var_dict
        self.run_param_dict = kwargs.pop("run_param_dict", None)
        self.run_file = kwargs.pop("run_file", None)
        self.electrode_params = electrode_params

        self.init_species(mwxrun.simulation.species)

        self.init_injectors_surfaces(injector_dict, surface_dict)

        # Store cross-sectional area for flux calculations
        self.area = kwargs.pop('area', None)
        self.geom = mwxrun.geom_str
        self.pos_lims = [mwxrun.xmin, mwxrun.xmax,
                         mwxrun.ymin, mwxrun.ymax,
                         mwxrun.zmin, mwxrun.zmax]
        self.nxyz = [mwxrun.nx, mwxrun.ny, mwxrun.nz]
        self.dt = mwxrun.get_dt()
        self.periodic = (
            #TODO: Should be implemented when needed
            #warp.w3d.boundxy is warp.periodic
            # The rwall in RZ simulations can be set even if boundxy is
            # periodic. It defaults to 1e36 so this should be fine.
            #and warp.top.prwall > 1e10
        )

        if self.area is None:
            self.area = mwxrun.get_domain_area()

        # All kwargs should have been popped, so raise an error if that isn't
        # the case.
        if len(kwargs) > 0:
            raise ValueError(
                f'Unrecognized kwarg(s) {list(kwargs.keys())}'
            )
        # Shapes/injectors just adds all the lists of shapes together into one
        # list
        self.shapes = reduce(
            (lambda x, y: x + y), list(self.surface_dict.values()))
        self.injectors = reduce(
            (lambda x, y: x + y), list(self.injector_dict.values()))
        # _reserved_names ensures results dictionary will not have duplicated
        # entries that overwrite each other.
        self._reserved_names = (["cathode_all", "anode_all",
                                 "aside", "aside_nodiel"] +
                                [x + "_all" for x in self.component_list])
        self._check_voltages_wfs(check_names=True)

        # Split out saving other stuff into functions
        self._save_sim_info()

    def init_species(self, species_list):
        """Handle creation and checking of species-related data."""
        self.species_list = species_list

    def init_injectors_surfaces(self, injector_dict, surface_dict):
        """Handle creation and checking of injectors & surfaces."""
        # Handle injector and surface dicts
        if 'cathode' not in injector_dict:
            raise ValueError("'cathode' must be specified in injector_dict")
        if 'cathode' not in surface_dict:
            raise ValueError("'cathode' must be specified in surface_dict")

        self.component_list = []
        # Ordering: We start with all known components in default order, then go
        # to the ordering of injector_dict and surface_dict for any remaining
        # labels.
        for cname in list(self.component_labels.keys()):
            if cname in injector_dict or cname in surface_dict:
                self.component_list.append(cname)

        for cdict in [injector_dict, surface_dict]:
            for cname in list(cdict.keys()):
                if cname not in self.component_list:
                    print(("Adding non-standard component {} to runinfo "
                           "data.").format(cname))
                    self.component_list.append(cname)

        # Now create ordered dictionaries internally for ease-of-use.
        self.injector_dict = collections.OrderedDict()
        self.surface_dict = collections.OrderedDict()

        for cname in self.component_list:
            if cname in injector_dict:
                self.injector_dict[cname] = mwxutil.return_iterable(
                    injector_dict[cname])
            if cname in surface_dict:
                self.surface_dict[cname] = mwxutil.return_iterable(
                    surface_dict[cname])

    @staticmethod
    def _process_local_vars(local_vars):
        """Return a new dictionary with uppercase entries.

        Arguments:
            local_vars (dict): A dictionary obtained by calling locals()

        Returns:
            user_var_dict (dict): A dictionary with uppercase entries, and
                copies of their values.
        """
        # locals() may or may not be a copy of the namespace. We want stuff
        # passed to RunInfo to be static, so we'll copy it all.
        user_var_dict = {}
        for key in list(local_vars.keys()):
            # IPython uses variables starting with _ for internal logic, eg '_'
            # or '_3', so we check for and avoid those too.
            if (key == key.upper()) and (key != 'MPI') and (key[0] != '_'):
                user_var_dict[key] = copy.deepcopy(local_vars[key])
        return user_var_dict

    def _check_voltages_wfs(self, check_names=False):
        """Check that for each electrode, all voltages/WFs/names are present.

        Allows no WF present for dielectric and focus grids only. V and WF must
        be equal for all cathode shapes.

        Arguments:
            check_names (bool): If True, make sure names are not in
                ``_reserved_names``, and update ``_reserved_names`` accordingly.
        """
        for component in self.surface_dict:
            for shape in self.surface_dict[component]:
                if shape.name == '':
                    raise ValueError(f'No name assigned for shape {shape}')
                if check_names:
                    if shape.name in self._reserved_names:
                        raise ValueError(
                            f'Repeated or reserved name {shape.name} '
                            f'used for shape {shape}'
                        )
                    self._reserved_names.append(shape.name)
                if shape.WF == 0.0:
                    raise ValueError(f'Work function for shape {shape.name}'
                                     f'; assign WF to this shape.'
                    )
                if mwxrun.me == 0:
                    if callable(shape.V):
                        voltage = shape.V(mwxrun.get_t())
                    else:
                        voltage = shape.V
                    print(f'Saved shape {shape.name} of type {component} '
                          f'with WF {shape.WF} eV and voltage {voltage:.2f} V '
                          f'(at t = {mwxrun.get_t()*1e9:.2f} ns)'
                    )

    def _save_sim_info(self):
        """Sets up the run_param_dict with helpful information.

        Note:
            In general, no guarantee is given that all elements are present, so
            use .get() functions with defaults if processing these elements of
            run_param_dict.
        """
        if self.run_param_dict is None:
            self.run_param_dict = {}
        # Main run info
        self.run_param_dict['datetime'] = datetime.datetime.utcnow()
        self.run_param_dict['nproc'] = mwxrun.n_procs
        # TODO: Implement when needed
        # self.run_param_dict['decomp'] = warp.warpoptions.options.decomp
        if self.run_file is not None:
            with open(self.run_file, 'r') as rfile:
                self.run_param_dict['run_file'] = rfile.read()
        # Other info
        # TODO: warpx doesn't have version info yet
        #self._save_version_info()
        self._save_comp_info()

    def _save_version_info(self):
        """Save info about file versions."""
        self.run_param_dict['mewarpx_version'] = mewarpx.__version__
        # Version info is a tuple rather than a string
        self.run_param_dict['mewarpx_version_info'] = mewarpx.__version_info__
        cwd = os.getcwd()
        try:
           os.chdir(mwxutil.mewarpx_dir)
           self.run_param_dict['mewarpx_gitver'] = subprocess.check_output(
               ['git', 'describe'])
           self.run_param_dict['mewarpx_commit'] = subprocess.check_output(
               ['git', 'log', '-n', '1', '--pretty=%h'])
        except subprocess.CalledProcessError as e:
            logger.warning(
                f'Failed to retrieve git information with error {e}'
            )
        finally:
           os.chdir(cwd)

    def _save_comp_info(self):
        """Save info about the computer."""
        self.run_param_dict['uname'] = platform.uname()
        if platform.system() == 'Linux':
            try:
                self.run_param_dict['cpuinfo'] = subprocess.check_output(
                    ['cat', '/proc/cpuinfo'])
            except subprocess.CalledProcessError as e:
                logger.warning(
                    f'Failed to retrieve processor information with error {e}'
                )
            try:
                self.run_param_dict['lscpu'] = subprocess.check_output(
                    ['lscpu'])
            except subprocess.CalledProcessError as e:
                logger.warning(
                    f'Failed to retrieve lscpu information with error {e}'
                )
            try:
                # https://stackoverflow.com/questions/20010199/how-to-determine-if-a-process-runs-inside-lxc-docker
                cgroupinfo = subprocess.check_output(
                    ['cat', '/proc/1/cgroup'])
                self.run_param_dict['cgroupinfo'] = cgroupinfo
                self.run_param_dict['ecs'] = b'ecs' in cgroupinfo
                self.run_param_dict['docker'] = (b'docker' in cgroupinfo
                                                 or self.run_param_dict['ecs'])
            except subprocess.CalledProcessError as e:
                logger.warning(
                    f'Failed to retrieve docker information with error {e}'
                )

    def save(self, diagdir='diags', filename='runinfo.dpkl'):
        """Save run info as a pickle file into the diagnostics directory. The
        inclusion of injectors forces the use of dill rather than pickle, both
        owing to customization of Norcross ionization but also apparently
        scipy.interp1d.

        Arguments:
            diagdir (str): Folder to save to, default 'diags'
            filename (str): File to save to, default 'runinfo.dpkl'
        """
        if mwxrun.me == 0:
            filepath = os.path.join(diagdir, filename)
            if not os.path.isdir(diagdir):
                mwxutil.mkdir_p(diagdir)
            with open(filepath, 'w+b') as pfile:
                dill.dump(self, pfile)
