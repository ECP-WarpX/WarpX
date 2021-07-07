.. _developers-repo-structure:

WarpX Structure
===============

Repo Organization
-----------------

All the WarpX source code is located in ``Source/``.
All sub-directories have a pretty straightforward name.
The PIC loop is part of the WarpX class, in function ``WarpX::EvolveEM`` implemented in ``Source/WarpXEvolveEM.cpp``.
The core of the PIC loop (i.e., without diagnostics etc.) is in ``WarpX::OneStep_nosub`` (when subcycling is OFF) or ``WarpX::OneStep_sub1`` (when subcycling is ON, with method 1).

Code organization
-----------------

The main WarpX class is WarpX, implemented in ``Source/WarpX.cpp``.

Build System
------------

WarpX uses the :ref:`CMake build system generator <building-cmake>`.
Each sub-folder contains a file ``CMakeLists.txt`` with the names of the source files (``.cpp``) that are added to the build.
Do not list header files (``.H``) here.

For experienced developers, we also support AMReX' :ref:`GNUmake build script collection <developers-gnumake>`.
The file ``Make.package`` in each sub-folder has the same purpose as the ``CMakeLists.txt`` file, please add new ``.cpp`` files to both dirs.

C++ Includes
------------

All WarpX header files need to be specified relative to the ``Source/`` directory.

- e.g. ``#include "Utils/WarpXConst.H"``
- files in the same directory as the including header-file can be included with ``#include "FileName.H"``

By default, in a ``MyName.cpp`` source file we do not include headers already included in ``MyName.H``. Besides this exception, if a function or a class
is used in a source file, the header file containing its declaration must be included, unless the inclusion of a facade header is more appropriate. This is
sometimes the case for AMReX headers. For instance ``AMReX_GpuLaunch.H`` is a façade header for ``AMReX_GpuLaunchFunctsC.H`` and ``AMReX_GpuLaunchFunctsG.H``, which
contain respectively the CPU and the GPU implemetation of some methods, and which should not be included directly.
Whenever possible, forward declarations headers are included instead of the actual headers, in order to save compilation time (see dedicated section below). In WarpX forward
declaration headers have the suffix ``*_fwd.H``, while in AMReX they have the suffix ``*Fwd.H``.
The include order (see `PR #874 <https://github.com/ECP-WarpX/WarpX/pull/874#issuecomment-607038803>`__ and `PR #1947 <https://github.com/ECP-WarpX/WarpX/pull/1947>`__) and `proper quotation marks <https://gcc.gnu.org/onlinedocs/cpp/Include-Syntax.html>`__ are:

In a ``MyName.cpp`` file:

1. ``#include "MyName.H"`` (its header) then
2. (further) WarpX header files ``#include "..."`` then
3. WarpX forward declaration header files ``#include "..._fwd.H"``
4. AMReX header files ``#include <...>`` then
5. AMReX forward declaration header files ``#include <...Fwd.H>`` then
6. PICSAR header files ``#include <...>`` then
7. other third party includes ``#include <...>`` then
8. standard library includes, e.g. ``#include <vector>``

In a ``MyName.H`` file:

1. ``#include "MyName_fwd.H"`` (the corresponding forward declaration header, if it exists) then
2. WarpX header files ``#include "..."`` then
3. WarpX forward declaration header files ``#include "..._fwd.H"``
4. AMReX header files ``#include <...>`` then
5. AMReX forward declaration header files ``#include <...Fwd.H>`` then
6. PICSAR header files ``#include <...>`` then
7. other third party includes ``#include <...>`` then
8. standard library includes, e.g. ``#include <vector>``

Each of these groups of header files should ideally be sorted alphabetically, and a blank line should be placed between the groups.

For details why this is needed, please see `PR #874 <https://github.com/ECP-WarpX/WarpX/pull/874#issuecomment-607038803>`_, `PR #1947 <https://github.com/ECP-WarpX/WarpX/pull/1947>`_, the `LLVM guidelines <https://llvm.org/docs/CodingStandards.html#include-style>`_, and `include-what-you-use <https://github.com/include-what-you-use/include-what-you-use/blob/master/docs/WhyIWYU.md>`_.

Forward Declaration Headers
---------------------------
Forward declarations can be used when a header file needs only to know that a given class exists, without any further detail (e.g., when only a pointer to an instance of
that class is used). Forward declaration headers are a convenient way to organize forward declarations. If a forward declaration is needed for a given class ``MyClass``, declared in ``MyClass.H``,
the forward declaration should appear in a header file named ``MyClass_fwd.H``, placed in the same folder containing ``MyClass.H``.
Below we provide a simple example:

``MyClass_fwd.H``:

.. code-block:: cpp

  class MyClass;

``MyClass.H``:

.. code-block:: cpp

   #include "MyClass_fwd.H"
   #include "someHeader.H"
   class MyClass{/* stuff */};

``MyClass.cpp``:

.. code-block:: cpp

   #include "MyClass.H"
   class MyClass{/* stuff */};

Usage: in ``SimpleUsage.H``

.. code-block:: cpp

   #include "MyClass_fwd.H"
   #include <memory>

