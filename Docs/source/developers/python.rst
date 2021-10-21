.. _development-python:

Python interface
================

The Python interface provides low and high level access to much of the data in WarpX.
With the low level access, a user has direct access to the underlying memory contained
in the MultiFabs and in the particle arrays.
The high level provides a more user friendly interface.

High level interface
--------------------

There are two python modules that provide convenient access to the fields and the particles.

Fields
~~~~~~

The ``fields`` module provides wrapper around most of the MultiFabs that are defined in the WarpX class.
For a list of all of the available wrappers, see the file ``Python/pywarpx/fields.py``.
For each MultiFab, there is a function that will return a wrapper around the data.
For instance, the function ``ExWrapper`` returns a wrapper around the ``x`` component of the MultiFab vector ``Efield_aux``.

.. code-block:: python

   from pywarpx import fields
   Ex = fields.ExWrapper()

By default, this wraps the MultiFab for level 0. The ``level`` argument can be specified for other levels.
By default, the wrapper only includes the valid cells. To include the ghost cells, set the argument ``include_ghosts=True``.

The wrapper provides access to the data via global indexing.
Using standard array indexing (with exceptions) with square brackets, the data can be accessed using indices that are relative to the full domain (across the MultiFab and across processors).
With multiple processors, the result is broadcast to all processors.
This example will return the ``Bz`` field at all points along ``x`` at the specified ``y`` and ``z`` indices.

.. code-block:: python

   from pywarpx import fields
   Bz = fields.BzWrapper()
   Bz_along_x = Bz[:,5,6]

The same global indexing can be done to set values. This example will set the values over a range in ``y`` and ``z`` at the
specified ``x``. The data will be scattered appropriately to the underlying FABs.

.. code-block:: python

   from pywarpx import fields
   Jy = fields.JyFPWrapper()
   Jy[5,6:20,8:30] = 7.

The code does error checking to ensure that the specified indices are within the bounds of the global domain.
Note that negative indices are handled differently than with numpy arrays because of the possibility of having ghost cells.
With ghost cells, the lower ghost cells are accessed using negative indices (since ``0`` is the index of the lower bound of the
valid cells). Without ghost cells, a negative index will always raise an out of bounds error since there are no ghost cells.

Under the covers, the wrapper object has a list of numpy arrays that have pointers to the underlying data, one array for each FAB.
When data is being fetched, it loops over that list to gather the data.
The result is then gathered among all processors.
Note that the result is not writeable, in the sense that changing it won’t change the underlying data since it is a copy.
When the data is set, using the global indexing, a similar process is done where the processors loop over their FABs and set the data at the appropriate indices.

The wrappers are always up to date since whenever an access is done (either a get or a set), the list of numpy arrays for the FABs is regenerated.
In this case, efficiency is sacrificed for consistency.

If it is needed, the list of numpy arrays associated with the FABs can be obtained using the wrapper method ``_getfields``.
Additionally, there are the methods ``_getlovects`` and ``_gethivects`` that get the list of the bounds of each of the arrays.

Particles
~~~~~~~~~

This is still in development.
