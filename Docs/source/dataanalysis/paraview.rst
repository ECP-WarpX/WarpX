.. _dataanalysis-visualization-paraview:

3D Visualization: ParaView
==========================

*TODO*

WarpX results can be visualized by ParaView, an open source visualization and analysis software.
ParaView can be downloaded and installed from httpshttps://www.paraview.org.

Tutorials
---------

* https://www.paraview.org/Wiki/The_ParaView_Tutorial
* https://www.youtube.com/results?search_query=paraview+introduction
* https://www.youtube.com/results?search_query=paraview+tutorial

openPMD
-------

WarpX' openPMD files can be visualized with ParaView 5.9+.
ParaView supports ADIOS1, ADIOS2 and HDF5 files, as it implements (like WarpX) against `openPMD-api <https://github.com/openPMD/openPMD-api>`__.

For openPMD output, WarpX automatically creates an ``.pmd`` file per diagnostics, which can be opened with ParaView.

Plotfiles (AMReX)
-----------------

*TODO* (cross-ref AMReX capabilities)
