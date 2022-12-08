:orphan:

WarpX
-----

WarpX is an advanced **electromagnetic Particle-In-Cell** code.

It supports many features including:

    - Perfectly-Matched Layers (PML)
    - Boosted-frame simulations
    - Mesh refinement

For details on the algorithms that WarpX implements, see the :ref:`theory section <theory>`.

WarpX is a *highly-parallel and highly-optimized code*, which can run on GPUs and multi-core CPUs, and includes load balancing capabilities.
In addition, WarpX is also a *multi-platform code* and runs on Linux, macOS and Windows.

.. _contact:

Contact us
^^^^^^^^^^

If you are starting using WarpX, or if you have a user question, please pop in our `Gitter chat room <https://gitter.im/ECP-WarpX/community>`__ and get in touch with the community.

The `WarpX GitHub repo <https://github.com/ECP-WarpX/WarpX>`__ is the main communication platform.
Have a look at the action icons on the top right of the web page: feel free to watch the repo if you want to receive updates, or to star the repo to support the project.
For bug reports or to request new features, you can also open a new `issue <https://github.com/ECP-WarpX/WarpX/issues>`__.

We also have a `discussion page <https://github.com/ECP-WarpX/WarpX/discussions>`__ on which you can find already answered questions, add new questions, get help with installation procedures, discuss ideas or share comments.

.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   section#installation,
   section#usage,
   section#theory,
   section#data-analysis,
   section#development,
   section#maintenance,
   section#epilogue {
       display:none;
   }
   </style>

.. toctree::
   :hidden:

   coc
   acknowledge_us
   highlights

Installation
------------
.. toctree::
   :caption: INSTALLATION
   :maxdepth: 1
   :hidden:

   install/users
   install/cmake
   install/hpc
..   install/changelog
..   install/upgrade

Usage
-----
.. toctree::
   :caption: USAGE
   :maxdepth: 1
   :hidden:

   usage/how_to_run
   usage/parameters
   usage/python
   usage/examples
   usage/pwfa
   usage/workflows
   usage/faq

Data Analysis
-------------
.. toctree::
   :caption: DATA ANALYSIS
   :maxdepth: 1
   :hidden:

   dataanalysis/formats
   dataanalysis/yt
   dataanalysis/openpmdviewer
   dataanalysis/openpmdapi
   dataanalysis/paraview
   dataanalysis/visit
   dataanalysis/visualpic
   dataanalysis/picviewer
   dataanalysis/backtransformed_diags
   dataanalysis/reduced_diags

Theory
------
.. toctree::
   :caption: THEORY
   :maxdepth: 1
   :hidden:

   theory/intro
   theory/picsar_theory
   theory/amr
   theory/PML
   theory/boosted_frame
   theory/input_output
   theory/collisions

Development
-----------
.. toctree::
   :caption: DEVELOPMENT
   :maxdepth: 1
   :hidden:

   developers/contributing
   developers/workflows
   developers/developers
   developers/doxygen
   developers/gnumake
   developers/faq
.. good to have in the future:
..   developers/repostructure

Maintenance
-----------
.. toctree::
   :caption: MAINTENANCE
   :maxdepth: 1
   :hidden:

   maintenance/release
   maintenance/performance_tests

Epilogue
--------
.. toctree::
   :caption: EPILOGUE
   :maxdepth: 1
   :hidden:

   glossary
   acknowledgements
