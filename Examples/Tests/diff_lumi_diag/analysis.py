#!/usr/bin/env python3

import os
import sys
import re 
import numpy as np
from scipy.constants import e, eV

sys.path.insert(1, "../../../../warpx/Regression/Checksum/")
import checksumAPI

sigmaz = 300e-6
sigmax = 500e-9
sigmay = 10e-9
Q = 1.2e10*e
GeV = 1e9*eV

fn = sys.argv[1]

def get_bins(fname):
    my_bins = []
    with open(fname) as f:
        for line in f: 
            if line.startswith('#'):
                data = line.split()
                for d in data[2:]:
                    m = re.search('=(.*)\(\)', d)
                    my_bins.append(m.group(1))
    return np.asarray(my_bins, dtype=float)

fname = 'diags/reducedfiles/DifferentialLuminosity_beam1_beam2.txt'
Ecom = (get_bins(fname)-1)/GeV # GeV
dEcom = Ecom[1] - Ecom[0] # GeV
dL_dEcom = np.loadtxt(fname)[-1,2:] # s^-1 m^-2 J^-1
dL_dEcom *= 1e-2 * GeV # s^-1 cm^-2 J^-1

print(dL_dEcom[np.nonzero(dL_dEcom)])

# Get name of the test
test_name = os.path.split(os.getcwd())[1]

# Run checksum regression test
checksumAPI.evaluate_checksum(test_name, fn, rtol=1e-2)
