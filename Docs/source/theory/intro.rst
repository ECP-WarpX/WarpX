.. _theory:

Introduction
============

.. figure:: Plasma_acceleration_sim.png
   :alt: Plasma laser-driven (top) and charged-particles-driven (bottom) acceleration (rendering from 3-D Particle-In-Cell simulations). A laser beam (red and blue disks in top picture) or a charged particle beam (red dots in bottom picture) propagating (from left to right) through an under-dense plasma (not represented) displaces electrons, creating a plasma wakefield that supports very high electric fields (pale blue and yellow). These electric fields, which can be orders of magnitude larger than with conventional techniques, can be used to accelerate a short charged particle beam (white) to high-energy over a very short distance.

   Plasma laser-driven (top) and charged-particles-driven (bottom) acceleration (rendering from 3-D Particle-In-Cell simulations). A laser beam (red and blue disks in top picture) or a charged particle beam (red dots in bottom picture) propagating (from left to right) through an under-dense plasma (not represented) displaces electrons, creating a plasma wakefield that supports very high electric fields (pale blue and yellow). These electric fields, which can be orders of magnitude larger than with conventional techniques, can be used to accelerate a short charged particle beam (white) to high-energy over a very short distance.

Computer simulations have had a profound impact on the design and understanding of past and present plasma acceleration experiments :cite:p:`Tsung2006,Geddes2008,Geddes2009,Geddes2010,Huang2009`, with
accurate modeling of wake formation, electron self-trapping and acceleration requiring fully kinetic methods (usually Particle-In-Cell) using large computational resources due to the wide range of space and time scales involved. Numerical modeling complements and guides the design and analysis of advanced accelerators, and can reduce development costs significantly. Despite the major recent experimental successes :cite:p:`Leemans2014,Blumenfeld2007,Bulanov2014,Steinke2016`, the various advanced acceleration concepts need significant progress to fulfill their potential. To this end, large-scale simulations will continue to be a key component toward reaching a detailed understanding of the complex interrelated physics phenomena at play.

For such simulations,
the most popular algorithm is the Particle-In-Cell (or PIC) technique,
which represents electromagnetic fields on a grid and particles by
a sample of macroparticles.
However, these simulations are extremely computationally intensive, due to the need to resolve the evolution of a driver (laser or particle beam) and an accelerated beam into a structure that is orders of magnitude longer and wider than the accelerated beam.
Various techniques or reduced models have been developed to allow multidimensional simulations at manageable computational costs: quasistatic approximation :cite:p:`Sprangle1990,Antonsen1992,Krall1993,Mora1997,Huang2006`,
ponderomotive guiding center (PGC) models :cite:p:`Antonsen1992,Krall1993,Huang2006,Benedetti2010,Cowan2011`, simulation in an optimal Lorentz boosted frame :cite:p:`Vay2007,Bruhwiler2009,Vay2009a,Vay2009b,Vay2010,Martins2010a,Martins2010b,Martins2010c,Vay2011a,Vay2011b,Vay2011c,Yu2016`,
expanding the fields into a truncated series of azimuthal modes :cite:p:`Godfrey1985,Lifschitz2009,Davidson2015,Lehe2016,Andriyash2016`, fluid approximation :cite:p:`Krall1993,Shadwick2009,Benedetti2010` and scaled parameters :cite:p:`CormierMichel2009`.

.. bibliography::
