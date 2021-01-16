Visualizing a distribution mapping
==================================

WarpX provides via :ref:`reduced diagnostics <reduced-diagnostics>` a means for
visualization of :ref:`distribution mappings and load balance costs <loadbalancecosts>`.
Here we demonstrate the workflow for generating this data and using it to plot
distribution mappings and load balance costs.


Generating the data
-------------------

To generate 'Load Balance Costs' reduced diagnostics output, WarpX should be run
with the following lines added to the input file (the name of the reduced diagnostics
file, `LBC`, and interval in steps to output reduced diagnostics data, `100`, may
be changed as needed):

.. code-block:: python

    warpx.reduced_diags_names = LBC
    LBC.type = LoadBalanceCosts 
    LBC.frequency = 100

The line `warpx.reduced_diags_names = LBC` sets the name of the reduced diagnostics
output file to `LBC`.  The next line `LBC.type = LoadBalanceCosts` tells WarpX
that the reduced diagnostics is a `LoadBalanceCosts` diagnostic, and instructs
WarpX to record costs and rank layouts.  The final line, `LBC.frequency = 100`,
controls the interval for output of this reduced diagnostic's data.
