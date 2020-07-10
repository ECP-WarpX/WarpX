import os
import shutil
import numpy as np
import yt
import scipy.constants as scc
import glob

"""
This file is part of the suite of scripts to use LibEnsemble on top of WarpX
simulations. It reads the output plotfiles of a simulation and returns beam
quantity relevant for the LibEnsemble optimizations ('f') as well as other
beam quantities for convenience.
"""

yt.funcs.mylog.setLevel(50)

def slice_emittance(x, xp, g, b, w):
    xav = np.average(x, aweights = w)
    xpav = np.average(xp * g * b, aweights = w)
    xstd2 = np.average((x - xav)**2, aweights = w)
    xpstd2 = np.average((xp * g * b - xpav)**2, aweights = w)
    xp2 = (np.average((x - xav) * (xp * g * b-xpav), weights = w))**2
    em = np.sqrt(xstd2 * pstd2 - xp2)
    return em

def _beam_properties(filepath):
    """
    Reads plotfile filepath and compute and return beam parameters

    Parameters
    ----------
    filepath : path to plotfile to read
    """

    # Read beam quantities from plotfile
    ds = yt.load(filepath)
    ad = ds.all_data()
    w = ad['beam', 'particle_weight'].v
    x = ad['beam', 'particle_position_x'].v
    ux = ad['beam', 'particle_momentum_x'].v/scc.m_e/scc.c
    uy = ad['beam', 'particle_momentum_y'].v/scc.m_e/scc.c
    uz = ad['beam', 'particle_momentum_z'].v/scc.m_e/scc.c

    # Compute beam parameters
    # Defined like that, the beam charge is > 0.
    charge = np.sum(w) * scc.e
    gamma = np.sqrt(1. + ux**2 + uy**2 + uz**2)
    beta = beta = np.sqrt(1.0 - 1.0 / gamma**2)
    energy_MeV = scc.physical_constants["electron mass energy equivalent in MeV"][0] * (gamma - 1.)
    energy_avg = np.average(energy_MeV, weights = w)
    energy_std = np.average((energy_MeV - energy_avg)**2, weights = w) / energy_avg
    nslices = 20
    zslices = np.linspace(np.min(z), np.max(z), nslices+1)
    exlist = np.zeros(nslices)
    for slicei in range(nslices):
        cond = [z > zslices[slicei], z < zslices[slicei+1]]
        xslice = np.select(cond, x)
        wslice = np.select(cond, w)
        xpslice = np.select(cond, ux/uz)
        gslice = np.select(cond, gamma)
        bslice = np.select(cond, beta)
        exlist[slicei] = slice_emittance(x, xpslice, gslice, bslice, wslice)
    emittance = np.mean(exlist)
    return charge, energy_avg, energy_std, emittance


def read_sim_output(workdir):
    """
    Return optimizing quantity 'f' and other parameters for convenience.

    Parameters
    ----------
    workdir : Path to directory where the simulation ran.
    """
    # Get beam properties at the beginning of the run
    datafile = 'diags/plotfiles/plt00000/'
    filepath = os.path.join(workdir, datafile)
    charge_i, _, _, emittance_i = _beam_properties(filepath)

    # Get beam properties at the end of the run
    file_list = glob.glob('diags/plotfiles/plt?????')
    file_list.sort()
    datafile = file_list[-1]
    filepath = os.path.join(workdir, datafile)
    charge_f, energy_avg, energy_std, emittance_f = _beam_properties(filepath)

    # delete simulation results, just to have smaller data
    shutil.rmtree('diags')

    # Build a quantity to minimize (f) that encompasses
    # emittance AND charge loss 1% charge loss has the
    # same impact as doubling the initial emittance.
    # we minimize f!
    f = emittance_f + emittance_i*(1.-charge_f/charge_i)*100
    warpx_out = np.array([f, energy_std, energy_avg, charge_f, emittance_f])

    return warpx_out
