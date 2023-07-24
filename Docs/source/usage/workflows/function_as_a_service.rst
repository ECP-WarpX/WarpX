.. _function_as_a_service:

Function-as-a-Service (FaaS)
============================

This section is work in process ...

In this section, we explore how start-to-end scientific workflows with WarpX can be run in a comfortable *function-as-a-service* (*FaaS*) approach from the comfort of a Jupyter notebook on your local machine.
We are using `Globus Compute <https://globus-compute.readthedocs.io/en/latest/index.html>`__ (formerly `funcX <https://funcx.org/>`__), a powerful *FaaS* platform that enables seamless and efficient execution of scientific applications on remote computing resources.
With Globus Compute, you can effortlessly launch WarpX simulations as functions, enabling high-performance computation without the need to manage the underlying infrastructure.
A particular advantage is that Globus Compute decouples the cloud-hosted management functionality from the edge-hosted execution functionality.
Its software can be deployed, by users or administrators, on arbitrary laptops, clouds, clusters, and supercomputers, in effect turning them into function serving systems.

The strategy is that for each major task type in the workflow we deploy a Globus Compute `endpoint <https://globus-compute.readthedocs.io/en/latest/endpoints.html>`__.
We could have separate endpoints for simulations runs, e.g. spawning SLURM jobs on a supercomputer, doing analysis, visualization, and archiving.

Installation
------------

As an example, we set up a simple workflow on :ref:`building-perlmutter`.
It is useful to create a dedicated Python environment and encapsulate loading the WarpX dependencies and this environment from a profile file.


.. code-block:: sh

    source perlmutter_gpu_warpx_globus_compute.profile

We proceed to install Globus Compute into our environment.

.. code-block:: sh

    python3 -m pip install globus-compute-sdk
    python3 -m pip install globus-compute-endpoint


Deploying Endpoints
-------------------

Globus Compute provides a tutorial endpoint to quickly try out first remote computations.
However, we would like a dedicated endpoint for running simulations.
Let's write a YAML configuration file.

.. code-block:: yaml

    display_name: pm_sim_ep

    engine:
      type: HighThroughputEngine
      available_accelerators: 4
      max_workers_per_node: 4
      cores_per_worker: 32
      worker_debug: false

      address:
        ifname: nmnb0
        type: address_by_interface

      provider:
        # Slurm scheduler could be slow at times,
        # increase the command timeouts
        cmd_timeout: 120
        init_blocks: 0

        launcher:
          # Explicit process affinity binding
          overrides: --cpu-bind=cores
          type: SrunLauncher

        # Scale between 0-1 blocks with 2 nodes per block
        max_blocks: 1
        min_blocks: 0
        nodes_per_block: 1
        partition: # Partition / QOS

        # string to prepend to #SBATCH blocks in the submit
        # script to the scheduler
        scheduler_options: |
          #SBATCH -A ntrain4_g
          #SBATCH -C gpu
          ##SBATCH -q debug
          #SBATCH -q regular
          #SBATCH --reservation=NERSC_hackathon_2023_day4_gpu
          #SBATCH --gpu-bind=none
          #SBATCH --gpus-per-node=4

        type: SlurmProvider
        walltime: 00:10:00

        # Command to be run before starting a worker
        # e.g., loading our environment here
        worker_init: "source ${HOME}/perlmutter_gpu_warpx_globus_compute.profile"


Running Simulations
-------------------

Analysis
--------

... to be written ...
