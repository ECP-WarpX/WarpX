Create a custom laser using lasy package
==================================

LASY (LAser SYmple manipulator) is a Python package that can be used to define complex laser pulse profiles
for which either the analytical form or the experimental measurement is known.
Most of the contents of this user guide can be found in `lasy docs <https://lasydoc.readthedocs.io/en/latest/>`_ and one can see examples at ``Examples/Tests/laser_injection_from_file``.

Install lasy package and dependencies
-------------------

To install the code you will need to first clone the repository to your local machine. Change into the new directory and then run the install command as given below.

.. code-block:: python

    git clone https://github.com/LASY-org/lasy.git
    cd lasy
    python3 -m pip install -v .

Run your first example
-------------------

We will try a simple example by generating a Gaussian pulse at focus.

First let's load the required functions from the library:

.. code-block:: python

    from lasy.laser import Laser
    from lasy.profiles import GaussianProfile

Next, we define the physical parameters of the laser pulse and create the laser profile object:

.. code-block:: python

    wavelength     = 1.e-6    # Laser wavelength in meters
    polarization   = (1,0)    # Linearly polarized in the x direction
    energy         = 1.0      # Energy of the laser pulse in joules
    spot_size      = 12e-6    # Waist of the laser pulse in meters
    pulse_duration = 10e-15   # Pulse duration of the laser in seconds
    t_peak         = 0.0      # Location of the peak of the laser pulse in time

    laser_profile = GaussianProfile(wavelength, polarization, energy, spot_size, pulse_duration, t_peak)

Now we create a full laser object containing the above physical parameters together with the computational settings:

.. code-block:: python

    dimensions     = "xyt"                              # Use 3D geometry
    lo             = (-25e-6, -25e-6, -20e-15)          # Lower bounds of the simulation box
    hi             = (+25e-6, +25e-6, +20e-15)          # Upper bounds of the simulation box
    num_points     = (100, 100, 100)                    # Number of points (nx,ny,nt) in each dimension

    laser = Laser(dimensions, lo, hi, num_points, laser_profile)
    laser.normalize(energy, kind="energy")              # Optionnal: normalize energy of the laser pulse contained in grid

Finally, we create a lasy file containing the laser profile, using openPMD standard:

.. code-block:: python

    file_prefix    = 'gaussianlaser3d' # The file name will start with this prefix
    file_format    = 'h5'          # Format to be used for the output file

    laser.write_to_file(file_prefix, file_format)

or simply:

.. code-block:: python

    laser.write_to_file("gaussianlaser3d")  # Use h5 format by default


That's it! The laser pulse profile has been created and the complexe envelope is now stored in a lasy file called ``gaussianlaser3d.h5`` and can be used with either a 3D,2D,RZ or 1D WarpX executable!
Now let's take a look at the laser parameters that need to be specified in the inputs parameters:

.. code-block:: python

    #######################################################################################################################################
    ######################################################## INPUTS FILE ##################################################################
    #######################################################################################################################################
    lasers.names        = lasy_laser
    lasy_laser.profile      = from_file                     # Specify that we want to read the lasy file instead of a handwritten laser
    lasy_laser.lasy_file_name = "gaussianlaser3d.h5"        # Name of the lasy file
    lasy_laser.time_chunk_size = 50                         # Read the lasy file in time chunks (accepts a value between 2 and nt)
    lasy_laser.position     = 0. 0. 0.                      # This point is on the laser plane
    lasy_laser.direction    = 0. 0. 1.                      # The plane normal direction
    lasy_laser.polarization = 0. 1. 0.                      # The main polarization vector
    lasy_laser.e_max        = 1.e14                         # Maximum amplitude of the laser field (in V/m)
    lasy_laser.wavelength = 1.0e-6                          # The wavelength of the laser (in meters)
    lasy_laser.delay = 0.0                                  # Delay (>0) or anticipate (<0) the laser by a specific amount of time

    # Other inputs parameters

or inside a PICMI script:

.. code-block:: python




Customize your laser profile by using combined longitunal and transverse profiles
-------------------

lasy allows you to define a custom laser pulse profile using different longitunal and transverse profiles.

.. code-block:: python

    from lasy.laser import Laser
    from lasy.profiles import (
        CombinedLongitudinalTransverseProfile,
        GaussianProfile,
    )
    from lasy.profiles.longitudinal import GaussianLongitudinalProfile
    from lasy.profiles.transverse import LaguerreGaussianTransverseProfile

    # Units
    um = 1.e-6
    fs = 1.e-15

    # Parameters of the Laguerre Gaussian beam
    wavelength = 1.*um
    w0 = 12.*um
    tt = 10.*fs
    t_c = 20.*fs
    laser_energy = 1.0
    pol = (1, 0)

    # Create a Laguerre Gaussian laser in RZ geometry
    profile = CombinedLongitudinalTransverseProfile(
    wavelength,pol,laser_energy,
    GaussianLongitudinalProfile(wavelength, tt, t_peak=0),
    LaguerreGaussianTransverseProfile(w0, p=0, m=1),
    )
    dim = "rt"
    lo = (0e-6, -20e-15)
    hi = (+25e-6, +20e-15)
    npoints = (100,100)
    laser = Laser(dim, lo, hi, npoints, profile, n_azimuthal_modes=2)
    laser.normalize(laser_energy, kind="energy")
    laser.write_to_file("laguerrelaserRZ")

Customize your laser profile by using NumPy arrays
-------------------

Profile defined from a NumPy array directly (only supported for 3D arrays):

.. code-block:: python

    from lasy.laser import Laser
    from lasy.profiles import FromArrayProfile

    # Units
    um = 1.e-6

    # Parameters of the Laguerre Gaussian beam
    wavelength = 1.*um
    pol = (1, 0)

    # Array of the electric field of the laser pulse.
    E_field = custom_numpy_array
    # Python dictionary containing the axes vectors. Keys are ‘x’, ‘y’, ‘t’. Values are the 1D arrays of each axis. array.shape = (axes[‘x’].size, axes[‘y’].size, axes[‘t’].size)
    axes
    laser_profile = FromArrayProfile(wavelength, pol, E_field, axes, axes_order=['x', 'y', 't'])

    dimensions     = "xyt"                              # Use 3D geometry
    lo             = (-25e-6, -25e-6, -20e-15)          # Lower bounds of the simulation box
    hi             = (+25e-6, +25e-6, +20e-15)          # Upper bounds of the simulation box
    num_points     = (100, 100, 100)                    # Number of points (nx,ny,nt) in each dimension

    laser = Laser(dimensions, lo, hi, num_points, laser_profile)
    laser.write_to_file("numpylaser3D")

Customize your laser profile by using an openPMD file
-------------------
Profile defined from an openPMD file:

.. code-block:: python

    from lasy.laser import Laser
    from lasy.profiles import FromOpenPMDProfile

    # Parameters of the Laguerre Gaussian beam
    pol = (1, 0)

    # Load the external openPMD file containing the laser pulse.
    path = './diags/' # Path to openPMDTimeSeries
    iteration = 450
    field = 'E'

    FromOpenPMDProfile(path, iteration, pol, field)

    dimensions     = "xyt"                              # Use 3D geometry
    lo             = (-25e-6, -25e-6, -20e-15)          # Lower bounds of the simulation box
    hi             = (+25e-6, +25e-6, +20e-15)          # Upper bounds of the simulation box
    num_points     = (100, 100, 100)                    # Number of points (nx,ny,nt) in each dimension

    laser = Laser(dimensions, lo, hi, num_points, laser_profile)
    laser.write_to_file("openPMDlaser3D")
