.. _visualization-ascent:

In situ Visualization with Ascent
=================================
Ascent is a system designed to meet the in-situ visualization and analysis needs of simulation code teams running multi-physics calculations on many-core HPC architectures.
It provides rendering runtimes that can leverage many-core CPUs and GPUs to render images of simulation meshes.

Compiling with GNU Make
-----------------------
After building and installing Ascent according to the instructions at `Building Ascent <https://ascent.readthedocs.io/en/latest/BuildingAscent.html>`_, you can enable support for it in WarpX by changing the line

.. code-block:: bash

   USE_ASCENT_INSITU=FALSE

in ``GNUmakefile`` to

.. code-block:: bash

   USE_ASCENT_INSITU=TRUE

Furthermore, you must ensure that either the :code:`ASCENT_DIR` shell environment variable contains the directory where Ascent is installed or you must specify this location when invoking make, i.e.,

.. code-block:: bash

   make -j 8 USE_ASCENT_INSITU=TRUE ASCENT_DIR=/path/to/ascent/install


Inputs File Configuration
-------------------------

Once WarpX has been compiled with Ascent support, it will need to be enabled and configured at runtime.
This is done using our usual inputs file (read with ``amrex::ParmParse``).
The supported parameters are part of the :ref:`FullDiagnostics <running-cpp-parameters-diagnostics>` with ``<diag_name>.format`` parameter set to ``ascent``.


Visualization/Analysis Pipeline Configuration
---------------------------------------------
Ascent uses the file :code:`ascent_actions.yaml` to configure analysis and visualization pipelines.
Ascent looks for the :code:`ascent_actions.yaml` file in the current working directory.

For example, the following :code:`ascent_actions.yaml` file extracts an isosurface of the field ``Ex`` for ``15`` levels and saves the resulting images to :code:`levels_<nnnn>.png`.
`Ascent Actions <https://ascent.readthedocs.io/en/latest/Actions/index.html>`_ provides an overview over all available analysis and visualization actions.

.. code-block:: json

    -
      action: "add_pipelines"
      pipelines:
        p1:
          f1:
            type: "contour"
            params:
               field: "Ex"
               levels: 15
    -
      action: "add_scenes"
      scenes:
        scene1:
          image_prefix: "levels_%04d"
          plots:
            plot1:
              type: "pseudocolor"
              pipeline: "p1"
              field: "Ex"

Here is another :code:`ascent_actions.yaml` example that renders isosurfaces and particles:

.. code-block:: json

    -
      action: "add_pipelines"
      pipelines:
        p1:
          f1:
            type: "contour"
            params:
               field: "Bx"
               levels: 3
    -
      action: "add_scenes"
      scenes:
        scene1:
          plots:
            plot1:
              type: "pseudocolor"
              pipeline: "p1"
              field: "Bx"
            plot2:
              type: "pseudocolor"
              field: "particle_electrons_Bx"
              points:
                radius: 0.0000005
          renders:
            r1:
              camera:
                azimuth: 100
                elevation: 10
              image_prefix: "out_render_3d_%06d"


Finally, here is a more complex :code:`ascent_actions.yaml` example that creates the same images as the prior example, but adds a trigger that creates a Cinema Database at cycle ``300``:

.. code-block:: json

    -
      action: "add_triggers"
      triggers:
        t1:
          params:
            condition: "cycle() == 300"
            actions_file: "trigger.yaml"
    -
      action: "add_pipelines"
      pipelines:
        p1:
          f1:
            type: "contour"
            params:
               field: "jy"
               iso_values: [ 1000000000000.0, -1000000000000.0]
    -
      action: "add_scenes"
      scenes:
        scene1:
          plots:
            plot1:
              type: "pseudocolor"
              pipeline: "p1"
              field: "jy"
            plot2:
              type: "pseudocolor"
              field: "particle_electrons_w"
              points:
                radius: 0.0000002
          renders:
            r1:
              camera:
                azimuth: 100
                elevation: 10
              image_prefix: "out_render_jy_part_w_3d_%06d"


When the trigger condition is meet, ``cycle() == 300``, the actions in :code:`trigger.yaml` are also executed:

.. code-block:: json

    -
      action: "add_pipelines"
      pipelines:
        p1:
          f1:
            type: "contour"
            params:
               field: "jy"
               iso_values: [ 1000000000000.0, -1000000000000.0]
    -
      action: "add_scenes"
      scenes:
        scene1:
          plots:
            plot1:
              type: "pseudocolor"
              pipeline: "p1"
              field: "jy"
            plot2:
              type: "pseudocolor"
              field: "particle_electrons_w"
              points:
                radius: 0.0000001
          renders:
            r1:
              type: "cinema"
              phi: 10
              theta: 10
              db_name: "cinema_out"

You can view the Cinema Database result by opening :code:`cinema_databases/cinema_out/index.html`.

Replay
-------------------------

`Replay <https://ascent.readthedocs.io/en/latest/Utilities.html#getting-data-for-replay>`_ is a utility that allows the user to ‘replay’ Conduit Blueprint HDF5 files saved by Ascent or exported by VisIt (starting in version 3.0 beta) back into Ascent. Replay enables the user or developer to pick specific time steps and load them into Ascent 

We'll guide you through the procedure of how to use Replay through a WarpX LWFA example

1. Get Blueprint Files
-------------------------
To use replay, you first need Conduit Blueprint HDF5 files which could be extracted by replay. The following ascent actions file can be used to extract Conduit Blueprint HDF5 files.

.. code-block:: json

    -
      action: action: "add_extracts"
      extracts:
        e1:
          type: "relay"
          params:
            path: "my_output_file_name"
            protocol: "blueprint/mesh/hdf5"

inputs_3d is LWFA input deck. The data size is 256X256X1024 and the simulation will run for 5000 cycles. WarpX will call ascent every 200 cyces. submit the following job submit script will extract 25 Buleprint HDF5:

.. code-block:: bash

   #!/bin/bash

   #BSUB -P aph114
   #BSUB -W 00:29
   #BSUB -nnodes 2
   #BSUB -alloc_flags smt4
   #BSUB -J WarpX
   #BSUB -o WarpXo.%J
   #BSUB -e WarpXe.%J

   module load gcc
   module load cuda

   export OMP_NUM_THREADS=1
   jsrun -r 6 -a 1 -g 1 -c 7 -l GPU-CPU -d packed -b rs --smpiargs="-gpu" ./warpx inputs_3d

2. example Actions
-------------------------

Images of Ex contour with particles could be obtained by replaying on the follow action:

.. code-block:: json

    -
      action: "add_pipelines"
      pipelines:
        clipped_volume:
          f0:
            type: "contour"
            params:
              field: "Ex"
              levels: 16
          f1:
            type: "clip"
            params:
              topology: topo # name of the amr mesh
              multi_plane:
                point1:
                  x: 0.0
                  y: 0.0
                  z: 0.0
                normal1:
                  x: 0.0
                  y: -1.0
                  z: 0.0
                point2:
                  x: 0.0
                  y: 0.0
                  z: 0.0
                normal2:
                  x: -0.7
                  y: -0.7
                  z: 0.0
        sampled_particles:
          f1:
            type: histsampling
            params:
              field: particle_electrons_uz
              bins: 64
              sample_rate: 0.90
          f2:
            type: "clip"
            params:
              topology: particle_electrons # particle data
              multi_plane:
                point1:
                  x: 0.0
                  y: 0.0
                  z: 0.0
                normal1:
                  x: 0.0
                  y: -1.0
                  z: 0.0
                point2:
                  x: 0.0
                  y: 0.0
                  z: 0.0
                normal2:
                  x: -0.7
                  y: -0.7
                  z: 0.0

    -
      action: "add_scenes"
      scenes:
        scene1:
          plots:
            p0:
              type: "pseudocolor"
              field: "particle_electrons_uz"
              pipeline: "sampled_particles"
            p1:
              type: "pseudocolor"
              field: "Ex"
              pipeline: "clipped_volume"
          renders:
            image1:
              bg_color: [1.0, 1.0, 1.0]
              fg_color: [0.0, 0.0, 0.0]
              image_prefix: "test%06d"
              camera:
                azimuth: 20
                elevation: 30
                zoom: 2.5
                                                                                                                                       
