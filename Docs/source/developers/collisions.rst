.. _developers-collisions:

Collisions
==========

Monte Carlo Collisions
----------------------

The Monte Carlo collisions (MCC) module can be found in *Source/Particles/Collisions/BackgroundMCC/*.
Several types of collisions between simulation particles and a neutral background gas are supported including elastic scattering, back scattering, charge exchange, excitation collisions and impact ionization.
An instance of the class :cpp:class:`MCCProcess` is created for each type of collision included in a simulation. This class saves information about the type of collision and the collision cross-section as a function of energy.

The so-called null collision strategy is used in order to minimize the computational burden of the MCC module. This strategy is standard in PIC-MCC and a detailed description can be found elsewhere, for example in `Birdsall (IEEE Transactions on Plasma Science, vol. 19, no. 2, pp. 65-85, 1991) <https://ieeexplore.ieee.org/document/106800>`_. In short the maximum collision probability is found over a sensible range of energies and is used to pre-select the appropriate number of macroparticles for collision consideration. Only these pre-selected particles are then individually considered for a collision based on their energy and the cross-sections of all the different collisional processes included.

The MCC implementation assumes that the background neutral particles are **thermal**, and are moving at non-relativistic velocities in the lab frame. For each simulation particle considered for a collision, a velocity vector for a neutral particle is randomly chosen. The particle velocity is then boosted to the stationary frame of the neutral through a Galilean transformation. The energy of the collision is calculated using the particle utility function, ``ParticleUtils::getCollisionEnergy()``, as

    .. math:: E_{coll} = \sqrt{(\gamma mc^2 + Mc^2)^2 - (mv)^2} - (mc^2 + Mc^2)

where :math:`v` is the speed of the particle as tracked in WarpX (so really :math:`\gamma` times the particle speed), :math:`m` and :math:`M` are the rest masses of the simulation and background species, respectively. The Lorentz factor is defined in the usual way, :math:`\gamma = \sqrt{1 + v^2/c^2}`. The collision cross-sections for all scattering processes are evaluated at the energy as calculated above.

Once a particle is selected for a specific collision process, that process determines how the particle is scattered as outlined below.

Charge exchange
^^^^^^^^^^^^^^^

This is the simplest scattering process. Under charge exchange the simulation particle's velocity is simply replaced with the sampled velocity of the neutral particle. Charge exchange normally has a cooling effect on the ions in the simulation by exchanging energetic ions for lower energy neutrals.

Elastic scattering
^^^^^^^^^^^^^^^^^^

This scattering process as well as the ones below that relate to it, are all performed in the center-of-momentum (COM) frame. To this end the particle velocity is rotated to a primed coordinate system such that its velocity is along the z' axis. This is achieved by rotating :math:`\tan(\theta) = -\frac{v_y}{v_x}` degrees onto the yz plane followed by a :math:`\tan(\phi) = \frac{v_x\cos(\theta) - v_y\sin(\theta)}{v_z}` rotation around the y-axis (see function ``ParticleUtils::getRotationAngles()``). Note that this rotation is not actually necessary to calculate since we know the result will give :math:`\vec{v}' = (0, 0, v)`, but it is necessary to store the angles in order to do the reverse rotation later. The Lorentz transformation to the COM frame is then simple to do since the particle velocity is aligned with the relative motion between the inertial frames:

    .. math:: u_{z} = \frac{u - v_{COM}}{1 - \frac{uv_{COM}}{c^2}}

where :math:`v_{COM} = mv/(\gamma m + M)` is the speed of the COM frame. Also note that since WarpX stores the particle velocity components multiplied by the Lorentz factor, it is necessary to divide by that factor before doing the Lorentz transformation. The particle velocity in this "collision" frame is then scattered by the appropriate angle based on the scattering model used. The Wentzel-Moliere model described in `Fernandez-Varea et al. (1992) <https://doi.org/10.1016/0168-583X(93)95827-R>`_ and `Lim (2007) <https://search.library.berkeley.edu/permalink/01UCS_BER/s4lks2/cdi_proquest_miscellaneous_35689087>`_ has been implemented to determine the elastic scattering angle as a function of collision energy. This model reduces to an isotropic model at low impact energy but favors small angle scattering at high energy as can be seen in the histogram below.

.. figure:: https://user-images.githubusercontent.com/40245517/169352031-9d66b9d9-c5ac-40a4-9d51-133b658c5cf6.png
   :alt: Wentzel-Moliere model
   :width: 80%

After the direction of the velocity vector has been appropriately changed, it is transformed back to the lab frame through a series of inverse transformations. Firstly, a Lorentz transform is done with relative velocity :math:`-v_{COM}` (see function ``ParticleUtils::doLorentzTransform()``), followed by the reverse rotation to realign the coordinate system with the lab frame. This is achieved via

    .. math::
        \vec{v} = \begin{bmatrix}
                \cos(\theta)\cos(\phi) & -\sin(\phi) & -\sin(\theta)\cos(\phi) \\
                \cos(\theta)\sin(\phi) & \cos(\phi) & -\sin(\theta)\sin(\phi) \\
                \sin(\theta) & 0 & \cos(\theta)
            \end{bmatrix} \begin{bmatrix}
                u_x \\
                u_y \\
                u_z
            \end{bmatrix}

Lastly, the reverse Galilean transformation is done using the starting neutral velocity to get back to the lab frame.

Back scattering
^^^^^^^^^^^^^^^

The process is the same as for elastic scattering above expect the scattering angle is fixed at :math:`\pi`.

Excitation
^^^^^^^^^^

The process is also the same as for elastic scattering except the excitation energy cost is subtracted from the particle energy. This is done by reducing the velocity before a scattering angle is chosen.

It is straight forward to determine the energy a projectile loses during an elastic collision with another body, as a function of scattering angle, through energy and momentum conservation. See for example `Lim (2007) <https://search.library.berkeley.edu/permalink/01UCS_BER/s4lks2/cdi_proquest_miscellaneous_35689087>`_ for a derivation. The result is that given a projectile with mass :math:`m`, a target with mass :math:`M`, a scattering angle :math:`\theta`, and collision energy :math:`E`, the post collision energy of the projectile is given by

    .. math:: E_{final} = E - [(E + mc^2)\sin^2\theta + Mc^2 - \cos(\theta)\sqrt{M^2c^4 - m^2c^4\sin^2\theta}] \\ \times\frac{E(E+2mc^2)}{(E+mc^2+Mc^2)^2 - E(E+2mc^2)\cos^2\theta}

The impact of incorporating relativistic effects in the MCC routine can be seen in the plots below where high energy collisions are considered with both a classical and relativistic implementation of MCC. It is observed that the classical version of MCC reproduces the classical limit of the above equation but especially for ions, this result differs substantially from the fully relativistic result.

.. figure:: https://user-images.githubusercontent.com/40245517/170900079-74e505a5-2790-44f5-ac84-5847eda954e6.png
   :alt: Classical v relativistic MCC
   :width: 96%
