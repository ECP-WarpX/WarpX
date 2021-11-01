Introduction
============

MEWarpX is a set of tools for running and postprocessing
`WarpX <https://warpx.readthedocs.io/en/latest/index.html>`_ simulations for
thermionic applications. This includes functionality that couples to WarpX's
python interface to inject *thermionically emitted* electrons at every simulation
step (see :class:`~mewarpx.emission.ThermionicInjector`) from an
:class:`~mewarpx.emission.Emitter`. Emitters can include both domain
boundaries and embedded conductors. The velocity of thermionically emitted
electrons are sampled from the distribution function derived in the supplemental
material of `Groenewald et. al. (2021)
<https://journals.aps.org/pre/abstract/10.1103/PhysRevE.103.023207>`_.
Furthermore, both the emitted and absorbed charge at every simulation step can
be recorded for all conductors in the simulation. This allows tracking of
electrical currents in the simulation as a function of time, through the
:mod:`~mewarpx.diags_store.flux_diagnostic` module.

Using a fully kinetic approach (such as a particle-in-cell code) to simulate the
underlying particle dynamics in a thermionic converter is important since the
particle populations have been observed to be strongly non-Maxwellian (see
page 23 on `Lietz et. al. (2021)
<http://plasma-pici-doe.umich.edu/files/LTP_2021_booklet_v10a.pdf>`_ results).

The package is openly available from the `MEWarpX Github repo
<https://github.com/ModernElectron/WarpX>`_.

For more details about thermionics and their applications as well as other
details about Modern Electron please visit the `company website
<https://modernelectron.com/>`_.

Usage
-----

To use Modern Electron WarpX Tools in a project::

    import mewarpx

Examples scripts of simulations that use the package can be found in the
`Examples <https://github.com/ModernElectron/WarpX/tree/memaster/mewarpx/examples>`_
directory. More usage information should come in the future.

Installation
------------
Prerequisites include: numpy, scipy, WarpX.

See :ref:`development-page` for information on how to do an efficient compile
and install of everything. That is the configuration we expect to consistently
use.

If you only want to install the python package component:

  * As a user, you can run ``make devel`` or ``make devel-all``. These use the
    ``pip`` package manager, but they also set it so that a symlink to the
    current directory is inserted into your python user environment; changes
    made to the source will automatically be used in new ``import mewarpx``
    calls.
  * In a virtual environment, run ``make devel-system`` or ``make
    devel-all-system`` for the same effect.

``make devel-all`` and ``make devel-all-system`` will also install all optional
dependencies. This, notably, will let you run ``make test``, or equivalently
``pytest``, to check that the whole distribution is running properly (which
requires WarpX compilation, however, again see :ref:`development-page` for
details).

For building all docs in latex format, ``latexmk``, and a number of latex
packages, must be installed as well. Installing a **full** ``texlive``
distribution will satisfy package requirements.

Compatibility
-------------

Targeting python 3.6+

Authors
-------

`mewarpx` was written by `Modern Electron <peter.scherpelz@modernelectron.com>`_.
