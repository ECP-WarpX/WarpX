import matplotlib.pyplot as plt
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
from scipy.optimize import curve_fit

ts = OpenPMDTimeSeries('./diags/diag1/')
dt = 1.27e-8
t=[]
phi=[]
it = ts.iterations

for i in it:
    phi_i = ts.get_field('phi',iteration=i,plot=False)
    # Find the minimum value among all grids for this iteration
    phi_min = np.min(phi_i[0])
    phi.append(phi_min)
    t.append(dt*i)


def func(x, v0, tau):

    return v0 * (1-np.exp(-np.array(x) / tau))



popt, pcov = curve_fit(func, t, phi)

plt.plot(t,phi, label='modelisation')
plt.plot(t, func(t, *popt), 'r-',label='fit: v0=%5.3f, tau=%5.9f' % (popt[0], popt[1]))
plt.legend()
plt.savefig('min_phi_analysis2.png')

print('parameters of the fit between the curve of min(phi) over the time and the function v0(1-exp(-t/tau)):')
print('v0=%5.3f, tau=%5.9f' % (popt[0], popt[1]))


tolerance_v0=0.1
tolerance_tau=0.1
print("tolerance for v0 = "+ str(tolerance_v0 *100) + '%')
print("tolerance for tau = "+ str(tolerance_tau*100) + '%')

moy_v0=-147.075
moy_tau=0.0000038725

diff_v0=np.abs((popt[0]-moy_v0)/moy_v0)
diff_tau=np.abs((popt[1]-moy_tau)/moy_tau)

print("pourcentage error for v0 = "+ str(diff_v0 *100) + '%')
print("pourcentage error for tau = "+ str(diff_tau*100) + '%')

assert (diff_v0 < tolerance_v0) and (diff_tau < tolerance_tau), 'Test spacecraft_charging did not pass'
