.. _examples-ohm-solver-magnetic-reconnection:

Ohm Solver: Magnetic Reconnection
=================================

Hybrid-PIC codes are often used to simulate magnetic reconnection in space plasmas.
An example of magnetic reconnection from a force-free sheet is provided, based on the simulation described in :cite:t:`ex-Le2016`.

Run
---

This example can be run as a **Python** script: ``python3 PICMI_inputs.py``.

.. literalinclude:: PICMI_inputs.py
   :language: python3
   :caption: You can copy this file from ``Examples/Tests/ohm_solver_magnetic_reconnection/PICMI_inputs.py``.


Analyze
-------

We run the following script to analyze correctness:

.. note::

   This section is TODO.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``analysis.py``

   .. literalinclude:: analysis.py
      :language: python3
      :caption: You can copy this file from ``Examples/Tests/ohm_solver_magnetic_reconnection/analysis.py``.

.. figure:: https://user-images.githubusercontent.com/40245517/229639784-b5d3b596-3550-4570-8761-8d9a67aa4b3b.gif
   :alt: Magnetic reconnection.
   :width: 70%

   Magnetic reconnection.
