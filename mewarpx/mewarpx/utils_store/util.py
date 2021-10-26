"""
Utility functions for mewarpx.
"""
import errno
import inspect
import os
import warnings
import collections
import logging

import numpy as np

from pywarpx import geometry

import mewarpx
from mewarpx.utils_store import mwxconstants as constants

logger = logging.getLogger(__name__)

# http://stackoverflow.com/questions/50499/in-python-how-do-i-get-the-path-and-name-of-the-file-t
mewarpx_dir = os.path.join(os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe()))), "..")


def check_version(script_version, script_physics_version):
    """Check that a run script is compatible with this mewarpx.

    If this mewarpx distribution is *lower* than that of the run script, throw
    an error. If this mewarpx distribution API or physics version is higher
    than that of the run script, throw a warning. Else do nothing.

    Arguments:
        script_version (tuple): A tuple of ints representing the mewarpx
            version of the run script: (API_version, feature_version,
            patch_version).
        script_physics_version (int): Integer representing the physics version of the run script.
    """
    mewarpx_version = mewarpx.__version_info__
    mewarpx_physics_version = mewarpx.__physics_version__

    # Tuple comparison works well in Python and does what we want!
    if mewarpx_version < script_version:
        raise ValueError(
            f"This version of mewarpx {mewarpx_version} is older than the "
            f"version {script_version} this script was designed for."
        )

    # I'm not sure of any instance where mewarpx physics version would be <
    # script physics version but software version would not be, but still
    # safest to do the check.
    if mewarpx_physics_version < script_physics_version:
        raise ValueError(
            f"This physics version of mewarpx {mewarpx_physics_version} is "
            f"older than the version {script_physics_version} this script was "
            "written for."
        )

    # Warnings only printed if API or physics versions are out of date.
    if mewarpx_version[0] > script_version[0]:
        logger.warning(
            f"This version of mewarpx {mewarpx_version} is a newer API "
            f"version than the version {script_version} this script was "
            "designed for. Incompatibilities may be present."
        )

    if mewarpx_physics_version > script_physics_version:
        logger.warning(
            f"This physics version of mewarpx {mewarpx_physics_version} is "
            f"newer than the version {script_physics_version} this script was "
            "written for. Results may be different now."
        )


def init_libwarpx(ndim, rz):
    """_libwarpx requires the geometry be set before importing.

    This complicates a lot of our code if we need to delay importing it - so
    instead we import it here.

    Very Bad Things will happen if ndim and rz here are different than is
    used in the rest of the simulation!

    Arguments:
        ndim (int): Number of dimensions. Ignored for RZ.
        rz (bool): True for RZ simulations, else False.
    """
    geometry.coord_sys = 1 if rz else 0
    geometry.prob_lo = [0]*ndim
    import pywarpx._libwarpx
    # This just quiets linters like pyflakes by using the otherwise-unused
    # variable
    assert pywarpx._libwarpx


def compute_step(simulation, interval=None):
    diag_periods = []
    for diag in simulation.diagnostics:
        diag_periods.append(diag.period)

    step_interval = None
    if interval:
        step_interval = interval
    elif (interval is None) and (len(diag_periods) > 0):
        step_interval = max(diag_periods)
    else:
        step_interval = simulation.max_steps

    for period in diag_periods:
        if step_interval % period != 0:
            warnings.warn(f'Diagnostic interval {step_interval} not divisible '
                          f'by the minimun diagnostic period {period}! Extra '
                          f'diagnostic data may be outputted!')
    return step_interval


def get_velocities(num_samples, T, m, emission_type='thermionic',
                   transverse_fac=1.0, rseed=None):
    """Generate array of random [vx, vy, vz] for cathode-emitted electrons.

    Arguments:
        num_samples (int): Number of particles to generate velocities for
        T (float): Temperature for the electrons (usually material temp) (K)
        m (float): Mass of elementary particle (kg)
        emission_type (str): Use "thermionic" for a thermionic emitter oriented
            along +zhat, and use "random" for a purely thermal distribution
            with no preferred direction. "half_maxwellian" is used at present
            for surface ionization, again along +zhat. Defaults to
            "thermionic".
        transverse_fac (float): Scale the particle's x and y average energies
            by this factor, scales z average energy to conserve total average
            particle energy in the distribution. Default 1., Min 0., Max 2.
        rseed (positive int): If specified, seed the random number generator.
            Used for testing. The random number generator is set back at the
            end of the function.

    Returns:
        velocities (np.ndarray): array of shape (num_samples, 3) with (vx, vy,
        vz) for each electron.
    """
    if (emission_type != 'thermionic') and not np.isclose(transverse_fac, 1.0):
        return ValueError('transverse_fac is a support argument only for '
                          'thermionic emissiion models!')

    if rseed is not None:
        nprstate = np.random.get_state()
        np.random.seed(rseed)
    sigma = np.sqrt(constants.kb_J * T / m)

    if transverse_fac < 0.:
        warnings.warn('WARNING: transverse_fac is out of bounds\n'
                      'Constraining to minimum value of 0.')
        beta = 0.
    elif transverse_fac > 2.:
        warnings.warn('WARNING: transverse_fac is out of bounds\n'
                      'Constraining to maximum value of 2.')
        beta = np.sqrt(2.)
    else:
        beta = np.sqrt(transverse_fac)

    alpha = np.sqrt(2. - beta**2.)

    # vx and vy follow Maxwellian distributions
    vx = sigma * np.random.randn(num_samples) * beta
    vy = sigma * np.random.randn(num_samples) * beta

    if emission_type == 'random':
        vz = sigma * np.random.randn(num_samples) * beta
    elif emission_type == 'thermionic':
        # vz is truncated with k_B*T/m << Phi assumed. See
        # "Emission Distribution from a Thermionic Cathode" document
        # on Bloomfire.
        P = np.random.rand(num_samples)
        vz = np.sqrt(-2 * sigma**2 * np.log(1 - P)) * alpha
    elif emission_type == 'half_maxwellian':
        vz = np.abs(sigma * np.random.randn(num_samples) * beta)
    else:
        raise ValueError(f'Unsupported emission type "{emission_type}".')

    if rseed is not None:
        np.random.set_state(nprstate)
    return vx, vy, vz


def get_positions(num_samples, xmin, xmax, ymin=0, ymax=0, z=0,
                  rseed=None):
    """Provide random samples of [x, y, z] for electrons in simulation.
    In x and y, positions are uniformly distributed. In z, positions are
    placed at the emitter.

    Arguments:
        num_samples (int): Number of particles to generate positions for
        xmin (float): Min position in x (meters)
        xmax (float): Max position in x (meters)
        ymin (float): Min position in y (meters)
        ymax (float): Max position in y (meters)
        z (float): Position of the emitter on the z-axis (meters)
        rseed (positive int): If specified, seed the random number generator.
            Used for testing. The random number generator is set back at the
            end of the function.

    Returns:
        positions (np.ndarray): Array of shape (num_samples, 3) with positions.
    """
    if rseed is not None:
        nprstate = np.random.get_state()
        np.random.seed(rseed)

    # Random x and y positions
    x = xmin + (xmax - xmin) * np.random.rand(num_samples)
    y = ymin + (ymax - ymin) * np.random.rand(num_samples)
    z = np.ones_like(x) * z

    if rseed is not None:
        np.random.set_state(nprstate)

    return x, y, z


def return_iterable(x, depth=1):
    """Return x if x is iterable, None if x is None, [x] otherwise.

    Useful for arguments taking either a list or single value. Strings are a
    special case counted as 'not iterable'.

    Arguments:
        depth (int): This many levels must be iterable. So if you need an
            iterable of an iterable, this is 2.
    """
    if x is None:
        return None
    elif depth > 1:
        # First make sure it's iterable to one less than the required depth.
        x = return_iterable(x, depth=depth-1)
        # Now check that it's iterable to the required depth. If not, we just
        # need to nest it in one more list.
        x_flattened = x
        while depth > 1:
            if all([(isinstance(y, collections.abc.Iterable)
                     and not isinstance(y, str))
                    for y in x_flattened]):
                x_flattened = [z for y in x_flattened for z in y]
                depth -= 1
            else:
                return [x]
        return x

    elif isinstance(x, str):
        return [x]
    elif isinstance(x, collections.abc.Iterable):
        return x
    else:
        return [x]


# https://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
def mkdir_p(path):
    """Make directory and parent directories if they don't exist.

    Do not throw error if all directories already exist.
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def ideal_gas_density(p, T):
    """Calculate neutral gas density (in 1/m^3) from the ideal gas law using
    pressure in Torr.

    Arguments:
        p (float): Gas pressure (Torr)
        T (float): Mean gas temperature (K)

    Returns:
        N (float): Number density of gas atoms/molecules (1/m^3)
    """
    return (p * constants.torr_cgs) / (constants.kb_cgs * T) * 1e6


def J_RD(T, WF, A):
    """Returns the Richardson-Dushmann thermionic emission given a temperature
    and effective work function. Constant coefficient of emission (A) is
    assumed.

    Arguments:
        T (float): temperature of the cathode in K
        WF (float): work function of the cathode in eV
        A (float): coefficient of emission in Amp/m^2/K^2

    Returns:
        J (float): current density in Amp/m^2
    """
    return A*T**2*np.exp(-1.*WF/(constants.kb_eV*T))


def recursive_update(d, u):
    """Recursively update dictionary d with keys from u.
    If u[key] is not a dictionary, this works the same as dict.update(). If
    u[key] is a dictionary, then update the keys of that dictionary within
    d[key], rather than replacing the whole dictionary.
    """
    # https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            d[k] = recursive_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def plasma_Debye_length(T, n):
    """Returns the thermal Debye length of a plasma.

    Arguments:
        T (float): plasma temperature in K
        n (float): plasma density in m^-3

    Returns:
        lambda (float): Debye length in m
    """
    return np.sqrt(
        constants.epsilon_0 * constants.kb_J * T / (n * constants.e**2)
    )
