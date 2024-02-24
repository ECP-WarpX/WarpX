#!/usr/bin/env python3

import os
import re
import sys
from scipy.special import erf
from scipy.constants import epsilon_0,pi
from openpmd_viewer import OpenPMDTimeSeries
import openpmd_api as io
import numpy as np 

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI


sigmaz = 300e-6
sigmax = 516e-9
sigmay = 7.7e-9
Q = -3.2e-9

def w(z):
    return np.exp(-z**2) * ( 1 + erf(1.j*z) )
    
def evaluate_E(x, y, z):
    ''' Basseti-Erskine formula https://cds.cern.ch/record/122227/files/198005132.pdf '''
    den = np.sqrt(2*(sigmax**2-sigmay**2))
    arg1 = (x+1j*y)/den
    term1 = w(arg1)
    arg2 = (x*sigmay/sigmax + 1j*y*sigmax/sigmay)/den
    term2 = -np.exp(-x**2/(2*sigmax**2)-y**2/(2*sigmay**2))*w(arg2)   
    factor = Q/(2.*np.sqrt(2.)*pi*epsilon_0*sigmaz*den)*np.exp(-z**2/(2*sigmaz**2))     
    E_complex = factor*(term1 + term2)
    return E_complex.imag, E_complex.real


fn = sys.argv[1]

path=os.path.join(sys.argv[1], 'diags', 'diag1')
ts = OpenPMDTimeSeries(path)

Ex, info = ts.get_field(field='E', coord='x', iteration=0, plot=False)
Ey, info = ts.get_field(field='E', coord='y', iteration=0, plot=False)
Ez, info = ts.get_field(field='E', coord='z', iteration=0, plot=False) 

grid_x = info.x[1:-1] 
grid_y = info.y[1:-1]
grid_z = info.z[1:-1]
                        
hnx = int(0.5*len(grid_x))
hny = int(0.5*len(grid_y))

for k, z in enumerate(grid_z,start=1):
    Ex_warpx = Ex[k,hny,1:-1]
    Ey_warpx = Ey[k,1:-1,hnx]
    
    Ex_theory = evaluate_E(grid_x, 0., z)[0]
    Ey_theory = evaluate_E(0., grid_y, z)[1]

    assert(np.allclose(Ex_warpx, Ex_theory, rtol=0.032, atol=0))
    assert(np.allclose(Ey_warpx, Ey_theory, rtol=0.029, atol=0))
    

# Get name of the test
test_name = os.path.split(os.getcwd())[1]

# Run checksum regression test
checksumAPI.evaluate_checksum(test_name, fn, rtol=1e-2)



