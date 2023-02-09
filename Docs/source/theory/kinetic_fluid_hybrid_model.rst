.. _theory-kinetic-fluid-hybrid-model:

Kinetic-fluid Hybrid Model
==========================

Many problems in plasma physics fall in a class where electron kinetics do not
play a critical role in the solution. Examples of such situations include the
study of collisionless magnetic reconnection and instabilities driven by ion
temperature anisotropy, to mention only two. For these kinds of problems the
computational cost of resolving the electron dynamics can be avoided by modeling
the electrons as a fluid rather than kinetic particles. The simulation resolution
can then be set by the ion time and length scales (commonly the ion cyclotron
period :math:`1/\Omega_i` and ion skin depth :math:`l_i`, respectively), which
can reduce the total simulation time drastically compared to a simulation that
has to resolve the electron Debye length and CFL-condition based on the speed
of light.

Many authors have described variations of the kinetic ion & fluid electron model,
generally refered to as PIC-fluid hybrid or just hybrid models. The implementation
in WarpX follows the description from `Winske et al. (2022) <https://arxiv.org/abs/2204.01676>`_.
The following description follows mostly from that reference.

Model
-----

The basic justification for the hybrid model is that the system to which it is
applied is dominated by ion kinetics, with ions moving much slower than electrons
and photons. In this scenario two critical approximations can be made, namely,
neutrality (:math:`n_e=n_i`) and the Maxwell-Ampere equation can be simplified by
neglecting the displacement current term (Darwin approximation, see `Nielson and Lewis (1976) <https://www.sciencedirect.com/science/article/abs/pii/B9780124608160500154>`_),
giving,

    .. math::

        \mu_0\vec{J} = \vec{\nabla}\times\vec{B}.

where :math:`\vec{J} = \vec{J}_i - \vec{J}_e` is the total electrical current i.e.
the sum of electron and ion currents. Seeing as ions are treated in the regular
PIC manner, the ion current, :math:`\vec{J}_i`, is known during a simulation. Therefore,
given the magnetic field, the electron current can be calculated.

If we now furher assume electrons are inertialess, the electron momentum
equation (Ohm's law) yields,

    .. math::

        \frac{d(n_em_e\vec{V}_e)}{dt} = 0 = -en_e\vec{E}-\vec{J}_e\times\vec{B}-\nabla\cdot\vec{P}_e+en_e\vec{\eta}\cdot\vec{J},

where :math:`\vec{V_e}=\vec{J}_e/(en_e)`, :math:`\vec{P}_e` is the electron pressure
tensor and :math:`\vec{\eta}` is the resistivity. An expression for the electric field
can be obtained from the above as:

    .. math::

        \vec{E} = -\frac{1}{en_e}\left( \vec{J}_e\times\vec{B} + \nabla\cdot\vec{P}_e \right)+en_e\vec{\eta}\cdot\vec{J}.

Lastly, if an electron temperature is given from which the electron pressure can
be calculated, the model is fully constrained.

Implementation details
----------------------
