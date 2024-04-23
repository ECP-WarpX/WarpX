.. _examples-beam-beam_collision:

Beam-beam collision
====================

This example shows how to simulate the collision between two ultra-relativistic particle beams.
This is representative of what happens at the interaction point of a linear collider.
We consider a right-propagating electron bunch colliding against a left-propagating positron bunch.

We turn on the Quantum Synchrotron QED module for photon emission (also known as beamstrahlung in the collider community) and
the Breit-Wheeler QED module for the generation of electron-positron pairs (also known as coherent pair generation in the collider community).

This solver computes the average velocity of each species, and solves the corresponding relativistic Poisson equation (see the WarpX documentation for ``warpx.do_electrostatic = relativistic`` for more details). This solver accurately reproduces the subtle cancellation that occurs for some components of ``E + v x B``, which is crucial in simulations of relativistic particles.


This example is based on the following paper :cite:t:`ex-Yakimenko2019`.


Run
---

The PICMI input file is not available for this example yet.

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. literalinclude:: inputs
   :language: ini
   :caption: You can copy this file from ``Examples/Physics_applications/beam-beam_collision/inputs``.


Visualize
---------

The figure below shows the number of photons emitted per beam particle (left) and the number of secondary pairs generated per beam particle (right). We present the reduced diagnostics using the ``plot_reduced.py`` script and then further compare different results for the reduced diagnostics with some literature:

* (red) simplified WarpX simulation as the example stored in the directory ``/Examples/Physics_applications/beam-beam_collision``;
* (blue) large-scale WarpX simulation (high resolution and ad hoc generated tables ;
* (black) literature results from :cite:t:`ex-Yakimenko2019`.

The small-scale simulation has been performed with a resolution of ``nx = 64, ny = 64, nz = 128`` grid cells, while the large-scale one has a much higher resolution of ``nx = 512, ny = 512, nz = 1024``. Moreover, the large-scale simulation uses dedicated QED lookup tables instead of the builtin tables. To generate the tables within WarpX, the code must be compiled with the flag ``-DWarpX_QED_TABLE_GEN=ON``. For the large-scale simulation we have used the following options:

.. code-block:: ini

   qed_qs.lookup_table_mode = generate
   qed_bw.lookup_table_mode = generate
   qed_qs.tab_dndt_chi_min=1e-3
   qed_qs.tab_dndt_chi_max=2e3
   qed_qs.tab_dndt_how_many=512
   qed_qs.tab_em_chi_min=1e-3
   qed_qs.tab_em_chi_max=2e3
   qed_qs.tab_em_chi_how_many=512
   qed_qs.tab_em_frac_how_many=512
   qed_qs.tab_em_frac_min=1e-12
   qed_qs.save_table_in=my_qs_table.txt
   qed_bw.tab_dndt_chi_min=1e-2
   qed_bw.tab_dndt_chi_max=2e3
   qed_bw.tab_dndt_how_many=512
   qed_bw.tab_pair_chi_min=1e-2
   qed_bw.tab_pair_chi_max=2e3
   qed_bw.tab_pair_chi_how_many=512
   qed_bw.tab_pair_frac_how_many=512
   qed_bw.save_table_in=my_bw_table.txt

.. figure:: https://user-images.githubusercontent.com/17280419/291749626-aa61fff2-e6d2-45a3-80ee-84b2851ea0bf.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDMwMzQzNTEsIm5iZiI6MTcwMzAzNDA1MSwicGF0aCI6Ii8xNzI4MDQxOS8yOTE3NDk2MjYtYWE2MWZmZjItZTZkMi00NWEzLTgwZWUtODRiMjg1MWVhMGJmLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFJV05KWUFYNENTVkVINTNBJTJGMjAyMzEyMjAlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjMxMjIwVDAxMDA1MVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWFiYzY2MGQyYzIyZGIzYzUxOWI3MzNjZTk5ZDM1YzgyNmY4ZDYxOGRlZjAyZTIwNTAyMTc3NTgwN2Q0YjEwNGMmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.I96LQpjqmFXirPDVnBlFQIkCuenR6IuOSY0OIIQvtCo
   :alt: Beam-beam collision benchmark against :cite:t:`ex-Yakimenko2019`.
   :width: 100%

   Beam-beam collision benchmark against :cite:t:`ex-Yakimenko2019`.

.. tab-set::

   .. tab-item:: Full Diagnostics

      This example can be run as a python script to visualize the fields evolution of the collision between two ultra-relativistic particle beams:

      - **Python** script: ``python3 plot_full.py``

      The python script loads WarpX simulation stored data (diags) using OpenPMDTimeSeries and iterates over each time step ``n = 65``, after which the fields ``E,B,rho`` components in ``x`` and ``y`` directions are extracted. There after, the plots to visualize the evolution of electric field ``E``, magnetic field ``B``, and charge density ``rho`` components of the two ultra-relativistic colliding particle beams are generated.

      .. code-block:: python

         # You can copy this file from Examples/Physics_applications/beam-beam_collisions/Plot_full_.py
         # Contents of plot_full_.py
         from openpmd_viewer import OpenPMDTimeSeries
         import numpy as np
         import matplotlib.pyplot as plt
         from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

         series = OpenPMDTimeSeries('./diags/diag1')
         steps = series.iterations
         print(steps)

         ylabel = 'y [m]'
         xlabel = 'z [m]'
         #slice_axis = 'x'
         slice_axis = 'y'

         #loop through E,B,Rho
         for n in steps:

             fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(40, 10), dpi=100., sharex=True, sharey=True)

             #E field
             Ex, info = series.get_field(field='E', coord='x', iteration=n, plot=False, slice_across=slice_axis)
             Ey, info = series.get_field(field='E', coord='y', iteration=n, plot=False, slice_across=slice_axis)
             Ez, info = series.get_field(field='E', coord='z', iteration=n, plot=False, slice_across=slice_axis)

             #B field
             Bx, info = series.get_field(field='B', coord=slice_axis, iteration=n, plot=False, slice_across=slice_axis)
             By, info = series.get_field(field='B', coord='y', iteration=n, plot=False, slice_across=slice_axis)
             Bz, info = series.get_field(field='B', coord='z', iteration=n, plot=False, slice_across=slice_axis)

             # Rho
             rho_beam1, info = series.get_field(field='rho_beam1', iteration=n, plot=False, slice_across=slice_axis)
             rho_beam2, info = series.get_field(field='rho_beam2', iteration=n, plot=False, slice_across=slice_axis)
             rho_ele1, info = series.get_field(field='rho_ele1', iteration=n, plot=False, slice_across=slice_axis)
             rho_pos1, info = series.get_field(field='rho_pos1', iteration=n, plot=False, slice_across=slice_axis)
             rho_ele2, info = series.get_field(field='rho_ele2', iteration=n, plot=False, slice_across=slice_axis)
             rho_pos2, info = series.get_field(field='rho_pos2', iteration=n, plot=False, slice_across=slice_axis)

             xmin = info.z.min()
             xmax = info.z.max()
             if slice_axis == 'x':
                 ymin = info.y.min()
                 ymax = info.y.max()
             elif slice_axis == 'y':
                 ymin = info.x.min()
                 ymax = info.x.max()


               #E field plots
             im1 = ax[0,0].imshow(np.transpose(Ex), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
             ax[0,0].set_title(f'E$_x$', fontsize=20)
             divider1 = make_axes_locatable(ax[0,0])
             cax = divider1.append_axes('right', size='5%', pad=0.05)
             fig.colorbar(im1, cax=cax, orientation='vertical')


             im2 = ax[0,1].imshow(np.transpose(Ey), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
             ax[0,1].set_title(f'E$_y$', fontsize=20)
             divider2 = make_axes_locatable(ax[0,1])
             cax = divider2.append_axes('right', size='5%', pad=0.05)
             fig.colorbar(im2, cax=cax, orientation='vertical')

             im3 = ax[0,2].imshow(np.transpose(Ez), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
             ax[0,2].set_title(f'E$_z$', fontsize=20)
             divider3 = make_axes_locatable(ax[0,2])
             cax = divider3.append_axes('right', size='5%', pad=0.05)
             fig.colorbar(im3, cax=cax, orientation='vertical')

             #B field plots
             im4 = ax[1,0].imshow(np.transpose(Bx), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
             ax[1,0].set_title(f'B$_x$', fontsize=20)
             divider4 = make_axes_locatable(ax[1,0])
             cax = divider4.append_axes('right', size='5%', pad=0.05)
             fig.colorbar(im4, cax=cax, orientation='vertical')


             im5 = ax[1,1].imshow(np.transpose(By), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
             ax[1,1].set_title(f'B$_y$', fontsize=20)
             divider5 = make_axes_locatable(ax[1,1])
             cax = divider5.append_axes('right', size='5%', pad=0.05)
             fig.colorbar(im5, cax=cax, orientation='vertical')

             im6 = ax[1,2].imshow(np.transpose(Bz), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
             ax[1,2].set_title(f'B$_z$', fontsize=20)
             divider6 = make_axes_locatable(ax[1,2])
             cax = divider6.append_axes('right', size='5%', pad=0.05)
             fig.colorbar(im6, cax=cax, orientation='vertical')

             #Rho plots
             im7 = ax[0,3].imshow(np.transpose(rho_beam1+rho_beam2), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
             ax[0,3].set_title(f'rho$_{{\tbeams{{}}}}$', fontsize=20)
             divider7 = make_axes_locatable(ax[0,3])
             cax = divider7.append_axes('right', size='5%', pad=0.05)
             fig.colorbar(im7, cax=cax, orientation='vertical')


             im8 = ax[1,3].imshow(np.transpose(rho_ele1+rho_pos1+rho_ele2+rho_pos2), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
             ax[1,3].set_title(f'rho$_{{\tsecondaries{{}}}}$', fontsize=20)
             divider8 = make_axes_locatable(ax[1,3])
             cax = divider8.append_axes('right', size='5%', pad=0.05)
             fig.colorbar(im8, cax=cax, orientation='vertical')


             ax[1,0].set_ylabel(ylabel, fontsize=20)
             ax[0,0].set_ylabel(ylabel, fontsize=20)
             ax[1,1].set_xlabel(xlabel, fontsize=20)
             ax[1,2].set_xlabel(xlabel, fontsize=20)
             ax[1,3].set_xlabel(xlabel, fontsize=20)
             ax[1,0].set_xlabel(xlabel, fontsize=20)


             fig.suptitle(f'Iteration {n:0}', fontsize=20)
             plt.tight_layout()

             image_file_name = 'FULL_'+slice_axis+f'{n:03d}.png'
             plt.savefig(image_file_name, dpi=100, bbox_inches='tight')
             plt.close()

   .. tab-item:: Reduced Diagnostics

      This example can be run as a Python file in the same directory:

      - **Python** script: ``python3 plot_reduced.py``
      The python script below was used to produce the reduced diagnostics which was then further benchmarked with Yakimenko.

      .. code-block:: python

         # You can copy this file from Examples/Physics_applications/beam-beam_collisions/Plot_reduced_.py
         # Contents of plot_reduced.py
         import numpy as np
         import matplotlib.pyplot as plt
         from scipy.constants import micron, c, pi, m_e, femto, e, milli, eV, physical_constants, alpha, nano
         from matplotlib import use, cm
         import matplotlib.colors
         import pandas as pd

         r_e = physical_constants['classical electron radius'][0]
         my_dpi=300
         sigmaz = 10*nano

         fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(2000./my_dpi, 1000./my_dpi), dpi=my_dpi)

         rdir= './diags/reducedfiles/'

         df_cr = pd.read_csv(f'{rdir}'+'ColliderRelevant_beam1_beam2.txt', sep=" ", header=0)
         df_pn = pd.read_csv(f'{rdir}'+'ParticleNumber.txt', sep=" ", header=0)


         times = df_cr[[col for col in df_cr.columns if f']time' in col]].to_numpy()
         steps = df_cr[[col for col in df_cr.columns if f']step' in col]].to_numpy()

         x = df_cr[[col for col in df_cr.columns if f']dL_dt' in col]].to_numpy()
         coll_index = np.argmax(x)
         coll_time = times[coll_index]

         # number of photons per beam particle
         np1 = df_pn[[col for col in df_pn.columns if f']pho1_weight' in col]].to_numpy()
         np2 = df_pn[[col for col in df_pn.columns if f']pho2_weight' in col]].to_numpy()
         Ne = df_pn[[col for col in df_pn.columns if f']beam1_weight' in col]].to_numpy()[0]
         Np = df_pn[[col for col in df_pn.columns if f']beam2_weight' in col]].to_numpy()[0]

         ax[0].plot((times-coll_time)/(sigmaz/c), (np1+np2)/(Ne+Np), lw=2)
         ax[0].set_title(r'photon number/beam particle')

         # number of NLBW particles per beam particle
         e1 = df_pn[[col for col in df_pn.columns if f']ele1_weight' in col]].to_numpy()
         e2 = df_pn[[col for col in df_pn.columns if f']ele2_weight' in col]].to_numpy()

         ax[1].plot((times-coll_time)/(sigmaz/c), (e1+e2)/(Ne+Np), lw=2)
         ax[1].set_title(r'NLBW particles/beam particle')

         for a in ax.reshape(-1):
             a.set_xlabel(r'time [$\sigma_z/c$]')
         image_file_name ='reduced.png'
         plt.tight_layout()
         plt.savefig(image_file_name,dpi=300, bbox_inches='tight')
         plt.close("all")
