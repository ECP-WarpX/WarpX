.. _theory-collisions:

Collisions
==========

Monte Carlo Collisions
----------------------

The Monte Carlo collisions (MCC) module can be found in *Source/Particles/Collisions/BackgroundMCC/*.
Several types of collisions between simulation particles and a neutral background gas are supported including elastic scattering, back scattering, charge exchange, excitation collisions and impact ionization.
An instance of the class :cpp:class:`MCCProcess` is created for each type of collision included in a simulation. This class saves information about the type of collision and the collision cross-section as a function of energy.

The so-called null collision strategy is used in order to minimize the computational burden of the MCC module. This strategy is standard in PIC-MCC and a detailed description can be found elsewhere, for example in `Birdsall (IEEE Transactions on Plasma Science, vol. 19, no. 2, pp. 65-85, 1991) <https://ieeexplore.ieee.org/document/106800>`_. In short the maximum collision probability is found over a sensible range of energies and is used to pre-select the appropriate number of macroparticles for collision consideration. Only these pre-selected particles are then individually considered for a collision based on their energy and the cross-sections of all the different collisional processes included.

The MCC implementation assumes that the background neutral particles are **thermal**, and are moving at non-relativistic velocities in the lab frame. For each simulation particle considered for a collision, a velocity vector for a neutral particle is randomly chosen. The particle velocity is then boosted to the stationary frame of the neutral through a Galilean transformation. The energy of the collision is calculated using the particle utility function, ``ParticleUtils::getCollisionEnergy()``, as

    .. math::

       \begin{aligned}
        E_{coll} &= \sqrt{(\gamma mc^2 + Mc^2)^2 - (mu)^2} - (mc^2 + Mc^2) \\
                 &= \frac{2Mmu^2}{M + m + \sqrt{M^2+m^2+2\gamma mM}}\frac{1}{\gamma + 1}
       \end{aligned}

where :math:`u` is the speed of the particle as tracked in WarpX (i.e. :math:`u = \gamma v` with :math:`v` the particle speed), while :math:`m` and :math:`M` are the rest masses of the simulation and background species, respectively. The Lorentz factor is defined in the usual way, :math:`\gamma = \sqrt{1 + u^2/c^2}`. Note that if :math:`\gamma\to1` the above expression clearly reduces to the classical equation :math:`E_{coll} = \frac{1}{2}\frac{Mm}{M+m} u^2`. The collision cross-sections for all scattering processes are evaluated at the energy as calculated above.

Once a particle is selected for a specific collision process, that process determines how the particle is scattered as outlined below.

Charge exchange
^^^^^^^^^^^^^^^

This is the simplest scattering process. Under charge exchange the simulation particle's velocity is simply replaced with the sampled velocity of the neutral particle. Charge exchange usually has a cooling effect on the ions in the simulation by exchanging energetic ions for lower energy neutrals.

Elastic scattering
^^^^^^^^^^^^^^^^^^

This scattering process as well as the ones below that relate to it, are all performed in the center-of-momentum (COM) frame. Designating the COM velocity of the particle as :math:`\vec{u}_c` and its labframe velocity as :math:`\vec{u}_l`, the transformation from lab frame to COM frame is done with a general Lorentz boost (see function ``ParticleUtils::doLorentzTransform()``):

    .. math::
            \begin{bmatrix}
                \gamma_c c \\
                u_{cx} \\
                u_{cy} \\
                u_{cz}
            \end{bmatrix}
         = \begin{bmatrix}
                \gamma & -\gamma\beta_x & -\gamma\beta_y & -\gamma\beta_z \\
                -\gamma\beta_x & 1+(\gamma-1)\frac{\beta_x^2}{\beta^2} & (\gamma-1)\frac{\beta_x\beta_y}{\beta^2} & (\gamma-1)\frac{\beta_x\beta_z}{\beta^2} \\
                -\gamma\beta_y & (\gamma-1)\frac{\beta_x\beta_y}{\beta^2} & 1 +(\gamma-1)\frac{\beta_y^2}{\beta^2} & (\gamma-1)\frac{\beta_y\beta_z}{\beta^2} \\
                -\gamma\beta_z & (\gamma-1)\frac{\beta_x\beta_z}{\beta^2} & (\gamma-1)\frac{\beta_y\beta_z}{\beta^2} & 1+(\gamma-1)\frac{\beta_z^2}{\beta^2} \\
            \end{bmatrix} \begin{bmatrix}
                \gamma_l c \\
                u_{lx} \\
                u_{ly} \\
                u_{lz}
            \end{bmatrix}

where :math:`\gamma` is the Lorentz factor of the relative speed between the lab frame and the COM frame, :math:`\beta_i = v^{COM}_i/c` is the i'th component of the relative velocity between the lab frame and the COM frame with

    .. math::

        \vec{v}^{COM} = \frac{m \vec{u_c}}{\gamma_u m + M}

The particle velocity in the COM frame is then isotropically scattered using the function ``ParticleUtils::RandomizeVelocity()``. After the direction of the velocity vector has been appropriately changed, it is transformed back to the lab frame with the reversed Lorentz transform as was done above followed by the reverse Galilean transformation using the starting neutral velocity.

Back scattering
^^^^^^^^^^^^^^^

The process is the same as for elastic scattering above expect the scattering angle is fixed at :math:`\pi`, meaning the particle velocity in the COM frame is updated to :math:`-\vec{u}_c`.

Excitation
^^^^^^^^^^

The process is also the same as for elastic scattering except the excitation energy cost is subtracted from the particle energy. This is done by reducing the velocity before a scattering angle is chosen.

Particle cooling due to elastic collisions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is straight forward to determine the energy a projectile loses during an elastic collision with another body, as a function of scattering angle, through energy and momentum conservation. See for example `Lim (2007) <https://search.library.berkeley.edu/permalink/01UCS_BER/s4lks2/cdi_proquest_miscellaneous_35689087>`_ for a derivation. The result is that given a projectile with mass :math:`m`, a target with mass :math:`M`, a scattering angle :math:`\theta`, and collision energy :math:`E`, the post collision energy of the projectile is given by

    .. math::

       \begin{aligned}
       E_{final} = E - &[(E + mc^2)\sin^2\theta + Mc^2 - \cos(\theta)\sqrt{M^2c^4 - m^2c^4\sin^2\theta}] \\
       &\times \frac{E(E+2mc^2)}{(E+mc^2+Mc^2)^2 - E(E+2mc^2)\cos^2\theta}
       \end{aligned}

The impact of incorporating relativistic effects in the MCC routine can be seen in the plots below where high energy collisions are considered with both a classical and relativistic implementation of MCC. It is observed that the classical version of MCC reproduces the classical limit of the above equation but especially for ions, this result differs substantially from the fully relativistic result.

.. figure:: https://user-images.githubusercontent.com/40245517/170900079-74e505a5-2790-44f5-ac84-5847eda954e6.png
   :alt: Classical v relativistic MCC
   :width: 96%
