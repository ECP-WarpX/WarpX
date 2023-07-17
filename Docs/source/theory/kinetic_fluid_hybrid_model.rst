.. _theory-kinetic-fluid-hybrid-model:

Kinetic-fluid Hybrid Model
==========================

Many problems in plasma physics fall in a class where both electron kinetics and electromagnetic waves do not
play a critical role in the solution. Examples of such situations include the
study of collisionless magnetic reconnection and instabilities driven by ion
temperature anisotropy, to mention only two. For these kinds of problems the
computational cost of resolving the electron dynamics can be avoided by modeling
the electrons as a neutralizing fluid rather than kinetic particles. By further
using Ohm's law to compute the electric field rather than evolving it with the
Maxwell-Faraday equation, light waves can be stepped over. The simulation resolution
can then be set by the ion time and length scales (commonly the ion cyclotron
period :math:`1/\Omega_i` and ion skin depth :math:`l_i`, respectively), which
can reduce the total simulation time drastically compared to a simulation that
has to resolve the electron Debye length and CFL-condition based on the speed
of light.

Many authors have described variations of the kinetic ion & fluid electron model,
generally referred to as particle-fluid hybrid or just hybrid-PIC models. The implementation
in WarpX follows the outline from :cite:t:`c-winske2022hybrid`.
This description follows mostly from that reference.

Model
-----

The basic justification for the hybrid model is that the system to which it is
applied is dominated by ion kinetics, with ions moving much slower than electrons
and photons. In this scenario two critical approximations can be made, namely,
neutrality (:math:`n_e=n_i`) and the Maxwell-Ampere equation can be simplified by
neglecting the displacement current term :cite:p:`c-NIELSON1976`, giving,

    .. math::

        \mu_0\vec{J} = \vec{\nabla}\times\vec{B},

where :math:`\vec{J} = \vec{J}_i - \vec{J}_e` is the total electrical current, i.e.
the sum of electron and ion currents. Since ions are treated in the regular
PIC manner, the ion current, :math:`\vec{J}_i`, is known during a simulation. Therefore,
given the magnetic field, the electron current can be calculated.

If we now further assume electrons are inertialess, the electron momentum
equation yields,

    .. math::

        \frac{d(n_em_e\vec{V}_e)}{dt} = 0 = -en_e\vec{E}-\vec{J}_e\times\vec{B}-\nabla\cdot\vec{P}_e+en_e\vec{\eta}\cdot\vec{J},

where :math:`\vec{V_e}=\vec{J}_e/(en_e)`, :math:`\vec{P}_e` is the electron pressure
tensor and :math:`\vec{\eta}` is the resistivity tensor. An expression for the electric field
(generalized Ohm's law) can be obtained from the above as:

    .. math::

        \vec{E} = -\frac{1}{en_e}\left( \vec{J}_e\times\vec{B} + \nabla\cdot\vec{P}_e \right)+\vec{\eta}\cdot\vec{J}.

Lastly, if an electron temperature is given from which the electron pressure can
be calculated, the model is fully constrained and can be evolved given initial
conditions.

Implementation details
----------------------

.. note::

    Various verification tests of the hybrid model implementation can be found in
    the :ref:`examples section <examples-hybrid-model>`.

The kinetic-fluid hybrid extension mostly uses the same routines as the standard electromagnetic
PIC algorithm with the only exception that the E-field is calculated from the
above equation rather than it being updated from the full Maxwell-Ampere equation. The
function ``WarpX::HybridPICEvolveFields()`` handles the logic to update the E&M fields
when the "hybridPIC" model is used. This function is executed after particle pushing
and deposition (charge and current density) has been completed. Therefore, based
on the usual time-staggering in the PIC algorithm, when ``HybridPICEvolveFields()`` is called
at timestep :math:`t=t_n`, the quantities :math:`\rho^n`, :math:`\rho^{n+1}`, :math:`\vec{J}_i^{n-1/2}`
and  :math:`\vec{J}_i^{n+1/2}` are all known.

Field update
^^^^^^^^^^^^

The field update is done in three steps as described below.

First half step
"""""""""""""""

Firstly the E-field at :math:`t=t_n` is calculated for which the current density needs to
be interpolated to the correct time, using :math:`\vec{J}_i^n = 1/2(\vec{J}_i^{n-1/2}+ \vec{J}_i^{n+1/2})`.
The electron pressure is simply calculated using :math:`\rho^n` and the B-field is also already
known at the correct time since it was calculated for :math:`t=t_n` at the end of the last step.
Once :math:`\vec{E}^n` is calculated, it is used to push :math:`\vec{B}^n` forward in time
(using the Maxwell-Faraday equation, i.e. the same as in the regular PIC routine with ``WarpX::EvolveB()``)
to :math:`\vec{B}^{n+1/2}`.

Second half step
""""""""""""""""

Next, the E-field is recalculated to get :math:`\vec{E}^{n+1/2}`. This is done
using the known fields :math:`\vec{B}^{n+1/2}`, :math:`\vec{J}_i^{n+1/2}` and
interpolated charge density :math:`\rho^{n+1/2}=1/2(\rho^n+\rho^{n+1})` (which is
also used to calculate the electron pressure). Similarly as before, the B-field
is then pushed forward to get :math:`\vec{B}^{n+1}` using the newly calculated
:math:`\vec{E}^{n+1/2}` field.

Extrapolation step
""""""""""""""""""

Obtaining the E-field at timestep :math:`t=t_{n+1}` is a well documented issue for
the hybrid model. Currently the approach in WarpX is to simply extrapolate
:math:`\vec{J}_i` foward in time, using

    .. math::

        \vec{J}_i^{n+1} = \frac{3}{2}\vec{J}_i^{n+1/2} - \frac{1}{2}\vec{J}_i^{n-1/2}.

With this extrapolation all fields required to calculate :math:`\vec{E}^{n+1}`
are known and the simulation can proceed.

Sub-stepping
^^^^^^^^^^^^

It is also well known that hybrid PIC routines require the B-field to be
updated with a smaller timestep than needed for the particles. The update steps
as outlined above are therefore wrapped in loops that enable the B-field to be
sub-stepped. The exact number of sub-steps used can be specified by the user
through a runtime simulation parameter (see :ref:`input parameters section <running-cpp-parameters-hybrid-model>`).

.. _theory-hybrid-model-elec-temp:

Electron pressure
^^^^^^^^^^^^^^^^^

The electron pressure is assumed to be a scalar quantity and calculated using the given
input parameters, :math:`T_{e0}`, :math:`n_0` and :math:`\gamma` using

    .. math::

        P_e = n_0T_{e0}\left( \frac{n_e}{n_0} \right)^\gamma.

The isothermal limit is given by :math:`\gamma = 1` while :math:`\gamma = 5/3`
(default) produces the adiabatic limit.

.. bibliography::
    :keyprefix: c-
