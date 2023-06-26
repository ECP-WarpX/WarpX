#!/usr/bin/env python3

# Copyright 2023 Juliette Pech
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""
This script displays an animation of the phase space over time, highlighting the appearance of the Weibel instability,
and compares the temporal evolution of energy with the associated theoretical growth rate.

Further analysis enables the experimental growth rate to be calculated using linear regression.
"""

from openpmd_viewer import OpenPMDTimeSeries
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy
from sklearn.linear_model import LinearRegression

ts = OpenPMDTimeSeries('diags/diag1')

t = ts.t
B_x_2 = []
B_y_2 = []
B_z_2 = []

for i in t:
    # Calculation of B squared
    B_x, info_x = ts.get_field( field='B', coord='x', t=i, plot=False )
    B_y, info_y = ts.get_field( field='B', coord='y', t=i, plot=False )
    B_z, info_z = ts.get_field( field='B', coord='z', t=i, plot=False )

    B_x_2 = np.append(B_x_2, np.sum(np.square(B_x)))
    B_y_2 = np.append(B_y_2, np.sum(np.square(B_y)))
    B_z_2 = np.append(B_z_2, np.sum(np.square(B_z)))

# Determination of the interval for which B_x squared can be assimilated to a
# straight line in logarithmic scale
def ln(x):
    ln = np.vectorize(np.log)
    return ln(x)

dt = t[1]-t[0]

# Derivative
dlnB_x_2dt = np.gradient(ln(B_x_2), dt)

# Cleans dlnB_x_2dt from non integer values (inf or NaN)
dlndt_cleaned = [x for x in dlnB_x_2dt if not (math.isinf(x) or math.isnan(x))]

# List of the indices for which the derivative of ln(B-x^2) is constant
indices = []
eps = 0.03e12  # Can be modified

for i in range(len(dlndt_cleaned)-1):
    if abs(dlndt_cleaned[i+1] - dlndt_cleaned[i]) < eps:
        indices.append(i)

# List of index intervals for which the derivative of ln(B-x^2) is constant
int_der_const = []
min = indices[0]
for i in range(len(indices)-2):
    if indices[i+1] - indices[i] < 4 and indices[i+2] - indices[i+1] >= 4:
        max = indices[i+1]
        int_der_const.append((min, max))
    elif indices[i+1] - indices[i] < 4 and indices[i+2] - indices[i+1] < 4 and i+2 != len(indices)-1:
        max = indices[i+2]
    elif indices[i+1] - indices[i] < 4 and indices[i+2] - indices[i+1] < 4 and i+2 == len(indices)-1: # end of the indices list
        max = indices[i+2]
        int_der_const.append((min, max))
    else :
        min = indices[i+1]

a = int_der_const[0][0]
b = int_der_const[0][1] + 1
X = t[a:b]
Y = B_x_2[a:b]

# Display of the growth rate gamma :
# f : t --> exp( 2*gamma*(t - t_0) )
beta = 0.01
# Plasma frequency : w_p = sqrt( n_0 * e^2 / m_e * eps_0 )
w_p = math.sqrt(2e25*(1.602e-19)**2/(9.109e-31*8.854e-12))
# Growth rate
gamma = beta * w_p

# Optimization of the t_0
a0 = 2*gamma
# b0 : the starting estimate for the roots of func(b) = 0
b0 = a0*1e-13

def func(i):
    return ((ln(Y) - (a0*X + i))**2).sum()
b0_opt = scipy.optimize.fsolve(func, b0)
# Other option : b0_opt = np.mean(ln(Y) - a0*X)
t_0 = - b0_opt / a0

def f(t):
    exp = np.vectorize(math.exp)
    return exp(2*gamma*(t-t_0))

# Linear regression
X = X.reshape((-1, 1))
model = LinearRegression().fit(X, ln(Y))
r_sq = model.score(X, ln(Y))

print(f"coefficient of determination: {r_sq}")
print(f"intercept: {model.intercept_}")
print(f"slope: {model.coef_}")
print("theoretical gamma times 2: {:e}".format(2*gamma))

gamma_opt = model.coef_/2
t_0_opt = - model.intercept_ / (2*gamma_opt)

def g(t):
    exp = np.vectorize(math.exp)
    return exp(2*gamma_opt*(t-t_0_opt))

# Plot with logarithmic scale
plt.figure()
plt.semilogy(t, B_x_2, label = 'B_x squared')
plt.semilogy(t, B_y_2, label = 'B_y squared')
plt.semilogy(t, B_z_2, label = 'B_z squared')
plt.semilogy(t, f(t), label = 'e^2*gamma*(t-t_0)')
plt.semilogy(t, g(t), label = 'e^2*gamma*(t-t_0) optimized')

plt.title('Logarithmic scale')
plt.xlabel('time (s)')
plt.ylabel('Bsquared')
plt.xlim(0,4e-12)
plt.ylim(0,1e8)
plt.legend()