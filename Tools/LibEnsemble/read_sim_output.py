import glob
import os
import shutil

import numpy as np
import scipy.constants as scc
import yt

"""
This file is part of the suite of scripts to use LibEnsemble on top of WarpX
simulations. It reads the output plotfiles of a simulation and returns beam
quantity relevant for the LibEnsemble optimizations ('f') as well as other
beam quantities for convenience.
"""

yt.funcs.mylog.setLevel(50)

def slice_emittance(x, xp, g, b, w):
    xpgb = xp * g * b
    xav = np.average(x, weights = w)
    xpav = np.average(xpgb, weights = w)
    xstd2 = np.average((x - xav)**2, weights = w)
    xpstd2 = np.average((xpgb - xpav)**2, weights = w)
    xp2 = (np.average((x - xav) * (xpgb-xpav), weights = w))**2
    em = np.sqrt(xstd2 * xpstd2 - xp2)
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

    if ('beam' in [i[0] for i in ds.field_list]):
        ad = ds.all_data()
        w = ad['beam', 'particle_weight'].v
        if (w.shape[0] <= 200):
            print('Insufficient particles in ',filepath)
            return 0.0, 0.0, 0.0, 0.0
        else:
            x = ad['beam', 'particle_position_x'].v
            z = ad['beam', 'particle_position_y'].v
            ux = ad['beam', 'particle_momentum_x'].v/scc.m_e/scc.c
            uy = ad['beam', 'particle_momentum_y'].v/scc.m_e/scc.c
            uz = ad['beam', 'particle_momentum_z'].v/scc.m_e/scc.c

            # Compute beam parameters
            # Defined like that, the beam charge is > 0.
            charge = np.sum(w) * scc.e
            gamma = np.sqrt(1. + ux**2 + uy**2 + uz**2)
            beta = np.sqrt(1.0 - 1.0 / gamma**2)
            energy_MeV = scc.physical_constants["electron mass energy equivalent in MeV"][0] * (gamma - 1.)
            energy_avg = np.average(energy_MeV, weights = w)
            energy_std = np.sqrt(np.average((energy_MeV - energy_avg)**2, weights = w)) / energy_avg * 100
            nslices = 20
            zslices = np.linspace(np.min(z), np.max(z), nslices+1)
            exlist = np.zeros(nslices)
            for slicei in range(nslices):
                cond = [(z > zslices[slicei]) & (z < zslices[slicei+1])][0]
                wslice = w[cond]
                if (wslice.shape[0] > 10):
                    xslice = x[cond]
                    xpslice = ux[cond]/uz[cond]
                    gslice = gamma[cond]
                    bslice = beta[cond]
                    exlist[slicei] = slice_emittance(xslice, xpslice, gslice, bslice, wslice)
            emittance = np.mean(exlist)
            return charge, energy_avg, energy_std, emittance
    else:
        print('No beam particles in ',filepath)
        return 0.0, 0.0, 0.0, 0.0

def read_sim_output(workdir):
    """
    Return optimizing quantity 'f' and other parameters for convenience.

    Parameters
    ----------
    workdir : Path to directory where the simulation ran.
    """
    # Get beam properties at the beginning of the run
    file_list = glob.glob('diags/plotfiles/plt?????')
    if (len(file_list) <2):
        print(workdir,' did not have final plotfile')
        return np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
    file_list.sort()
    datafile = file_list[0]
    filepath = os.path.join(workdir, datafile)
    charge_i, _, _, emittance_i = _beam_properties(filepath)

    # Get beam properties at the end of the run
    datafile = file_list[-1]
    filepath = os.path.join(workdir, datafile)
    charge_f, energy_avg, energy_std, emittance_f = _beam_properties(filepath)

    if (charge_f > 0.0):
        # delete simulation results, just to have smaller data
        shutil.rmtree('diags')

        # Build a quantity to minimize (f) that encompasses
        # emittance AND charge loss 1% charge loss has the
        # same impact as doubling the initial emittance.
        # we minimize f!
        f = emittance_f + emittance_i*(1.-charge_f/charge_i)*100
        warpx_out = np.array([f, energy_std, energy_avg, charge_f, emittance_f])

        return warpx_out
    else:
        return np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
