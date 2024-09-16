.. _accelerator_lattice:

Accelerator lattice
===================

The files in this directory handle the accelerator lattice. These are fields of various types and configurations.
The lattice is laid out along the z-axis.

The AcceleratorLattice has the instances of the accelerator element types and handles the input of the data.

The LatticeElementFinder manages the application of the fields to the particles. It maintains index lookup tables
that allow rapidly determining which elements the particles are in.

The classes for each element type are in the subdirectory LatticeElements.

Host and device classes
-----------------------

The LatticeElementFinder and each of the element types have two classes, one
that lives on the host and one that can be trivially copied to the device.
This dual structure is needed because of the complex data structures
describing both the accelerator elements and the index lookup tables. The
host level classes manage the data structures, reading in and setting up the
data. The host classes copy the data to the device and maintain the pointers
to that data on the device. The device level classes grab pointers to the
appropriate data (on the device) needed when fetching the data for the particles.

External fields
---------------

The lattice fields are applied to the particles from the GetExternalEBField
class. If a lattice is defined, the GetExternalEBField class gets the lattice
element finder device level instance associated with the grid being operated
on. The fields are applied from that instance, which calls the "get_field"
method for each lattice element type that is defined for each particle.

Adding new element types
------------------------

A number of places need to be touched when adding a new element types. The
best method is to look for every place where the "quad" element is referenced
and duplicate the code for the new element type. Changes will only be needed
within the AcceleratorLattice directory.
