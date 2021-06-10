mewarpx README
==============

A set of tools for running and postprocessing WarpX simulations for thermionic
applications.

Usage
-----

To use Modern Electron WarpX Tools in a project::

    import mewarpx

Then see :ref:`package-overview` or :ref:`usage-page` for more details.

Installation
------------
Prerequisites include: numpy, scipy, WarpX.

At this stage I recommend running ``make devel`` or ``make devel-all``. These
use the ``pip`` package manager, but they also set it so that a symlink to the
current directory is inserted into your python user environment; changes made
to the source will automatically be used in new ``import mewarpx`` calls.

``make devel-all`` will also install all optional dependencies. This, notably,
will let you run ``make test``, or equivalently ``pytest``, to check that the
whole distribution is running properly. See :ref:`testing-page` for more
details and options in testing.

For building all docs in latex format, ``latexmk``, and a number of latex
packages, must be installed as well. Installing a **full** ``texlive``
distribution will satisfy package requirements.


Compatibility
-------------

Targeting python 3.6+

Authors
-------

`mewarpx` was written by `Modern Electron <peter.scherpelz@modernelectron.com>`_.
