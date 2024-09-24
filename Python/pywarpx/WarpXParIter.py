"""
This file is part of WarpX

Copyright 2024 WarpX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def soa_real_comps(self, num_comps):
    """
    Name the WarpX ParticleReal components in SoA.

    Parameters
    ----------
    self : warpx_pybind module
      used to query freestanding functions
    num_comps : int
      number of components to generate names for.

    Returns
    -------
    A list of length num_comps with values.
    """
    names = self.real_comp_names

    if len(names) != num_comps:
        raise RuntimeError(
            "num_comps ({num_comps}) is not equal to length of real_comp_names ({len(names)})"
        )

    return names


def soa_int_comps(self, num_comps):
    """
    Name the WarpX int components in SoA.

    Parameters
    ----------
    self : warpx_pybind module
      used to query freestanding functions
    num_comps : int
      number of components to generate names for.

    Returns
    -------
    A list of length num_comps with values.
    """
    names = self.int_comp_names

    if len(names) != num_comps:
        raise RuntimeError(
            f"num_comps ({num_comps}) is not equal to length of int_comp_names ({len(names)})"
        )

    return names


def soa(self, warpx_pybind):
    """Get the StructOfArrays on the current tile"""
    soa = super(type(self), self).soa()

    # overwrite name providers
    soa.soa_real_comps = lambda num_comps: soa_real_comps(warpx_pybind, num_comps)
    soa.soa_int_comps = lambda num_comps: soa_int_comps(warpx_pybind, num_comps)

    return soa


def register_WarpXParIter_extension(warpx_pybind):
    """WarpXParIter helper methods"""

    # WarpXParticleContainer
    warpx_pybind.WarpXParIter.soa = lambda self: soa(self, warpx_pybind)
    # TODO: no constant variant used in WarpX?

    # PinnedParticleContainer
    # TOOD: do we want a custom iterator for NamedParticleContainer in general?
