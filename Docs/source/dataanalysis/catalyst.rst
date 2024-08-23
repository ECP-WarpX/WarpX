.. _visualization-catalyst:

In situ Visualization with Catalyst 2
=====================================
Catalyst 2 (further referred to as just Catalyst) is a lightweight in-situ visualization and analysis framework API developed for simulations and other scientific data producers. It has a lightweight implementation
(or **stub**) and an SDK to develop custom implementations of Catalyst. ParaView comes with its own implementation (known as **ParaView Catalyst**) for leveraging ParaView's
visualization and analysis capabilities, which is what this document will focus on.


Enabling Catalyst
-----------------
In order to use Catalyst with WarpX, you must `build Catalyst 2 <https://catalyst-in-situ.readthedocs.io/en/latest/build_and_install.html>`_ and `build <https://github.com/Kitware/ParaView/blob/master/Documentation/dev/build.md>`__ or `install <https://www.paraview.org/download/>`__ ParaView 5.9+. Afterward, AMReX must be built with ``AMReX_CONDUIT=TRUE``,
``AMReX_CATALYST=TRUE``, ``Conduit_DIR=/path/to/conduit``, and ``Catalyst_DIR=/path/to/catalyst`` (``/path/to/catalyst`` should be the directory containing ``catalyst-config.cmake``, not the path to the implementation).

Once AMReX is appropriately built, WarpX can be built with the following options:

.. code-block:: cmake

    WarpX_amrex_internal=FALSE
    AMReX_DIR="/path/to/amrex/build"

If they cannot be found, ``Conduit_DIR`` and ``Catalyst_DIR`` will have to be set again. Ensure that AMReX is built with all required options, some common ones being:

.. code-block:: cmake

    AMReX_MPI=TRUE
    AMReX_MPI_THREAD_MULTIPLE=TRUE
    AMReX_LINEAR_SOLVERS=TRUE
    AMReX_PARTICLES=TRUE
    AMReX_PARTICLES_PRECISION=DOUBLE
    AMReX_PIC=TRUE
    AMReX_TINY_PROFILE=TRUE


Inputs File Configuration
-------------------------
Once WarpX has been compiled with Catalyst support, it will need to be enabled and configured at runtime.
This is done using our usual inputs file (read with ``amrex::ParmParse``).
The supported parameters are part of the :ref:`FullDiagnostics <running-cpp-parameters-diagnostics>` with ``<diag_name>.format`` parameter set to ``catalyst``.

In addition to configuring the diagnostics, the following parameters must be included:
    * ``catalyst.script_paths``: The locations of the pipeline scripts, separated by either a colon or semicolon (e.g. ``/path/to/script1.py;/path/to/script2.py``).
    * ``catalyst.implementation`` (default ``paraview``): The name of the implementation being used (case sensitive).
    * ``catalyst.implementation_search_paths``: The locations to search for the given implementation. The specific file being searched for will be ``catalyst_{implementation}.so``.

Because the scripts and implementations are global, Catalyst does not benefit from nor differentiate between multiple diagnostics.


Visualization/Analysis Pipeline Configuration
---------------------------------------------
Catalyst uses the files specified in ``catalyst.script_paths`` to run all analysis.

The following script, :code:`simple_catalyst_pipeline.py`, automatically detects the type of data for both the mesh and particles, then creates an extractor for them. In most
cases, these will be saved as ``.VTPC`` files which can be read with the ``XML Partitioned Dataset Collection Reader``.

.. code-block:: python

    from paraview.simple import *
    from paraview import catalyst

    # Helper function
    def create_extractor(data_node, filename="Dataset"):
        VTK_TYPES = ["vtkImageData", "vtkRectilinearGrid", "vtkStructuredGrid", "vtkPolyData", "vtkUnstructuredGrid", "vtkUniformGridAMR", "vtkMultiBlockDataSet", "vtkPartitionedDataSet", "vtkPartitionedDataSetCollection", "vtkHyperTreeGrid"]
        FILE_ASSOCIATIONS = ["VTI", "VTR", "VTS", "VTP", "VTU", "VTH", "VTM", "VTPD", "VTPC", "HTG"]
        clientside_data = data_node.GetClientSideObject().GetOutputDataObject(0) # Gets the dataobject from the default output port

        # Loop is required because .IsA() detects valid classes that inherit from the VTK_TYPES
        for i, vtk_type in enumerate(VTK_TYPES):
            if (clientside_data.IsA(vtk_type)):
                filetype = FILE_ASSOCIATIONS[i]
                extractor = CreateExtractor(filetype, data_node, registrationName=f"_{filetype}")
                extractor.Writer.FileName = filename + "_{timestep:}" + f".{filetype}"
                return extractor

        raise RuntimeError(f"Unsupported data type: {clientside_data.GetClassName()}")

    # Camera settings
    paraview.simple._DisableFirstRenderCameraReset()  # Prevents the camera from being shown

    # Options
    options = catalyst.Options()

    options.CatalystLiveTrigger = "TimeStep"  # "Python", "TimeStep", "TimeValue"
    options.EnableCatalystLive = 0  # 0 (disabled), 1 (enabled)
    if (options.EnableCatalystLive == 1):
        options.CatalystLiveURL = "localhost:22222"  # localhost:22222 is default

    options.ExtractsOutputDirectory = "datasets"  # Base for where all files are saved
    options.GenerateCinemaSpecification = 0 # 0 (disabled), 1 (enabled), generates additional descriptor files for cinema exports
    options.GlobalTrigger = "TimeStep"  # "Python", "TimeStep", "TimeValue"

    meshSource = PVTrivialProducer(registrationName="mesh")  # "mesh" is the node where the mesh data is stored
    create_extractor(meshSource, filename="meshdata")
    particleSource = PVTrivialProducer(registrationName="particles")  # "particles" is the node where particle data is stored
    create_extractor(particleSource, filename="particledata")

    # Called on catalyst initialize (after Cxx side initialize)
    def catalyst_initialize():
        return

    # Called on catalyst execute (after Cxx side update)
    def catalyst_execute(info):
        print(f"Time: {info.time}, Timestep: {info.timestep}, Cycle: {info.cycle}")
        return

    # Callback if global trigger is set to "Python"
    def is_activated(controller):
        return True

    # Called on catalyst finalize (after Cxx side finalize)
    def catalyst_finalize():
        return

    if __name__ == '__main__':
        paraview.simple.SaveExtractsUsingCatalystOptions(options)


For the case of ParaView Catalyst, pipelines are run with ParaView's included ``pvbatch`` executable and use the ``paraview`` library to modify the data. While pipeline scripts
could be written manually, this is not advised for anything beyond the script above. It is much more practical to use ParaView's built in ``Save Catalyst State`` button.

The process for creating a pipeline is as follows:
    1. Run at least one step of simulation and save the data in a ParaView compatible format, then open it in ParaView.
    2. Set up the desired scene, including filters, camera and views, and extractors.
    3. Press ``Save Catalyst State``, or the multicolored flask icon in the top left corner, and save it to a desired location.
    4. Open the script and replace the used producer with ``PVTrivialProducer``, setting the ``registrationName`` to either ``mesh`` or ``particles`` based on what data is used.

As an example for step four, here are a few lines from a script directly exported from ParaView:

.. code-block:: python

    # create a new 'XML Image Data Reader'
    meshdatavti = XMLImageDataReader(registrationName='meshdata.vti', FileName=['/path/to/meshdata.vti'])
    meshdatavti.CellArrayStatus = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez']
    meshdatavti.TimeArray = 'None'

    # Calculator sample filter
    calculator1 = Calculator(registrationName='Calculator1', Input=meshdatavti)
    calculator1.AttributeType = 'Cell Data'
    calculator1.ResultArrayName = 'BField'
    calculator1.Function = 'sqrt(Bx^2 + By^2 + Bz^2)'

In order to use it with the mesh data coming from the simulation, the above code would be changed to:

.. code-block:: python

    # create the producer
    meshdata = PVTrivialProducer(registrationName='mesh')
    meshdata.CellArrayStatus = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez']
    meshdata.TimeArray = 'None'

    # Calculator sample filter
    calculator1 = Calculator(registrationName='Calculator1', Input=meshdata)
    calculator1.AttributeType = 'Cell Data'
    calculator1.ResultArrayName = 'BField'
    calculator1.Function = 'sqrt(Bx^2 + By^2 + Bz^2)'

Steps one is advised so that proper scaling and framing can be done, however in certain cases it may not be possible. If this is the case, a dummy object can be used instead
(such as a wavelet or geometric shape scaled appropriately) and the rest of the steps can be followed as usual.

Replay
------

Catalyst 2 supports replay capabilities, which can be read about `here <https://catalyst-in-situ.readthedocs.io/en/latest/catalyst_replay.html>`_.

.. note::

   * TODO: Add more extensive documentation on replay
