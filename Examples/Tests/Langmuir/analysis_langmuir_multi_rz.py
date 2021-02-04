#! /usr/bin/env python

# Copyright 2019 David Grote, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This is a script that analyses the simulation results from
# the script `inputs.multi.rz.rt`. This simulates a RZ periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_r = -\partial_r \phi = \epsilon \,\frac{mc^2}{e}\frac{2\,r}{w_0^2} \exp\left(-\frac{r^2}{w_0^2}\right) \sin(k_0 z) \sin(\omega_p t)
# $$ E_z = -\partial_z \phi = - \epsilon \,\frac{mc^2}{e} k_0 \exp\left(-\frac{r^2}{w_0^2}\right) \cos(k_0 z) \sin(\omega_p t)
# Unrelated to the Langmuir waves, we also test the plotfile particle filter function in this
# analysis script.
import sys
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
yt.funcs.mylog.setLevel(50)
import numpy as np
from scipy.constants import e, m_e, epsilon_0, c
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

test_name = fn[:-9] # Could also be os.path.split(os.getcwd())[1]

# Parse test name and check if current correction (psatd.current_correction) is applied
current_correction = True if re.search('current_correction', fn) else False

# Parameters (these parameters must match the parameters in `inputs.multi.rz.rt`)
epsilon = 0.01
n = 2.e24
w0 = 5.e-6
n_osc_z = 2
rmin =   0e-6; rmax = 20.e-6; Nr = 64
zmin = -20e-6; zmax = 20.e-6; Nz = 128

# Wave vector of the wave
k0 = 2.*np.pi*n_osc_z/(zmax-zmin)
# Plasma frequency
wp = np.sqrt((n*e**2)/(m_e*epsilon_0))
kp = wp/c

def Er( z, r, epsilon, k0, w0, wp, t) :
    """
    Return the radial electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Er_array = \
        epsilon * m_e*c**2/e * 2*r/w0**2 * \
            np.exp( -r**2/w0**2 ) * np.sin( k0*z ) * np.sin( wp*t )
    return( Er_array )

def Ez( z, r, epsilon, k0, w0, wp, t) :
    """
    Return the longitudinal electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Ez_array = \
        - epsilon * m_e*c**2/e * k0 * \
            np.exp( -r**2/w0**2 ) * np.cos( k0*z ) * np.sin( wp*t )
    return( Ez_array )

# Read the file
ds = yt.load(fn)
t0 = ds.current_time.to_value()
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                        dims=ds.domain_dimensions)

# Get cell centered coordinates
dr = (rmax - rmin)/Nr
dz = (zmax - zmin)/Nz
coords = np.indices([Nr, Nz],'d')
rr = rmin + (coords[0] + 0.5)*dr
zz = zmin + (coords[1] + 0.5)*dz

# Check the validity of the fields
overall_max_error = 0
Er_sim = data['Ex'].to_ndarray()[:,:,0]
Er_th = Er(zz, rr, epsilon, k0, w0, wp, t0)
max_error = abs(Er_sim-Er_th).max()/abs(Er_th).max()
print('Er: Max error: %.2e' %(max_error))
overall_max_error = max( overall_max_error, max_error )

Ez_sim = data['Ez'].to_ndarray()[:,:,0]
Ez_th = Ez(zz, rr, epsilon, k0, w0, wp, t0)
max_error = abs(Ez_sim-Ez_th).max()/abs(Ez_th).max()
print('Ez: Max error: %.2e' %(max_error))
overall_max_error = max( overall_max_error, max_error )

# Plot the last field from the loop (Ez at iteration 40)
plt.subplot2grid( (1,2), (0,0) )
plt.imshow( Ez_sim )
plt.colorbar()
plt.title('Ez, last iteration\n(simulation)')
plt.subplot2grid( (1,2), (0,1) )
plt.imshow( Ez_th )
plt.colorbar()
plt.title('Ez, last iteration\n(theory)')
plt.tight_layout()
plt.savefig(test_name+'_analysis.png')

error_rel = overall_max_error

if current_correction:
   tolerance_rel = 0.06
else:
   tolerance_rel = 0.04

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

# Check charge conservation (relative L-infinity norm of error) with current correction
if current_correction:
    divE = data['divE'].to_ndarray()
    rho  = data['rho' ].to_ndarray() / epsilon_0
    error_rel = np.amax(np.abs(divE - rho)) / max(np.amax(divE), np.amax(rho))
    tolerance = 1.e-9
    print("Check charge conservation:")
    print("error_rel = {}".format(error_rel))
    print("tolerance = {}".format(tolerance))
    assert( error_rel < tolerance )


## In the final past of the test, we verify that the diagnostic particle filter function works as
## expected in RZ geometry. For this, we only use the last simulation timestep.

## This is a generic function to test a particle filter. We reproduce the filter in python and
## verify that the results are the same as with the filtered diagnostic.
def check_particle_filter(fn, filtered_fn, filter_expression):
    ds  = yt.load( fn )
    ds_filtered  = yt.load( filtered_fn )
    ad  = ds.all_data()
    ad_filtered  = ds_filtered.all_data()

    ## Load arrays from the unfiltered diagnostic
    ids = ad['electrons', 'particle_id'].to_ndarray()
    r = ad['electrons', 'particle_position_x'].to_ndarray()
    z = ad['electrons', 'particle_position_y'].to_ndarray()
    theta = ad['electrons', 'particle_theta'].to_ndarray()
    px  = ad['electrons', 'particle_momentum_x'].to_ndarray()
    py  = ad['electrons', 'particle_momentum_y'].to_ndarray()
    pz  = ad['electrons', 'particle_momentum_z'].to_ndarray()
    w  = ad['electrons', 'particle_weight'].to_ndarray()

    ## Load arrays from the filtered diagnostic
    ids_filtered_warpx = ad_filtered['particle_id'].to_ndarray()
    r_filtered_warpx = ad_filtered['particle_position_x'].to_ndarray()
    z_filtered_warpx = ad_filtered['particle_position_y'].to_ndarray()
    theta_filtered_warpx = ad_filtered['particle_theta'].to_ndarray()
    px_filtered_warpx  = ad_filtered['particle_momentum_x'].to_ndarray()
    py_filtered_warpx  = ad_filtered['particle_momentum_y'].to_ndarray()
    pz_filtered_warpx  = ad_filtered['particle_momentum_z'].to_ndarray()
    w_filtered_warpx  = ad_filtered['particle_weight'].to_ndarray()

    ## Reproduce the filter in python: this returns the indices of the filtered particles in the
    ## unfiltered arrays.
    ind_filtered_python, = np.where(eval(filter_expression))

    ## Sort the indices of the filtered arrays by particle id.
    sorted_ind_filtered_python = ind_filtered_python[np.argsort(ids[ind_filtered_python])]
    sorted_ind_filtered_warpx = np.argsort(ids_filtered_warpx)

    ## Check that the sorted ids are exactly the same with the warpx filter and the filter
    ## reproduced in python
    assert(np.array_equal(ids[sorted_ind_filtered_python],
                             ids_filtered_warpx[sorted_ind_filtered_warpx]))

    ## Finally, we check that the sum of the particles quantities are the same to machine precision
    tolerance_checksum = 1.e-12
    check_array_sum(r[sorted_ind_filtered_python], r_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(z[sorted_ind_filtered_python], z_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(theta[sorted_ind_filtered_python], theta_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(px[sorted_ind_filtered_python], px_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(py[sorted_ind_filtered_python], py_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(pz[sorted_ind_filtered_python], pz_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(w[sorted_ind_filtered_python], w_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)

## This function checks that the absolute sums of two arrays are the same to a required precision
def check_array_sum(array1, array2, tolerance_checksum):
    sum1 = np.sum(np.abs(array1))
    sum2 = np.sum(np.abs(array2))
    assert(abs(sum2-sum1)/sum1 < tolerance_checksum)

## This function is specifically used to test the random filter. First, we check that the number of
## dumped particles is as expected. Next, we call the generic check_particle_filter function.
def check_random_filter(fn, filtered_fn, random_fraction):
    ds  = yt.load( fn )
    ds_filtered  = yt.load( filtered_fn )
    ad  = ds.all_data()
    ad_filtered  = ds_filtered.all_data()

    ## Check that the number of particles is as expected
    numparts = ad['electrons', 'particle_id'].to_ndarray().shape[0]
    numparts_filtered = ad_filtered['particle_id'].to_ndarray().shape[0]
    expected_numparts_filtered = random_fraction*numparts
    # 5 sigma test that has an intrinsic probability to fail of 1 over ~2 millions
    std_numparts_filtered = np.sqrt(expected_numparts_filtered)
    error = abs(numparts_filtered-expected_numparts_filtered)
    print("Random filter: difference between expected and actual number of dumped particles: " + str(error))
    print("tolerance: " + str(5*std_numparts_filtered))
    assert(error<5*std_numparts_filtered)

    random_filter_expression = 'np.isin(ids, ids_filtered_warpx)'
    check_particle_filter(fn, filtered_fn, random_filter_expression)

parser_filter_fn = "diags/diag_parser_filter00080"
parser_filter_expression = "(py-pz < 0) * (r<10e-6) * (z > 0)"
check_particle_filter(fn, parser_filter_fn, parser_filter_expression)

uniform_filter_fn = "diags/diag_uniform_filter00080"
uniform_filter_expression = "ids%3 == 0"
check_particle_filter(fn, uniform_filter_fn, uniform_filter_expression)

random_filter_fn = "diags/diag_random_filter00080"
random_fraction = 0.66
check_random_filter(fn, random_filter_fn, random_fraction)

checksumAPI.evaluate_checksum(test_name, fn)
