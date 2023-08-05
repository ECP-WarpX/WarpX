from matplotlib import cm, use
import matplotlib.colors
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import openpmd_api as io
import pandas as pd
from scipy.constants import c, e, hbar, m_e

E_crit = m_e**2*c**3/(e*hbar)
B_crit = m_e**2*c**2/(e*hbar)

def parse_input_file(input_file):
    """
    Parse WarpX input file.

    Parameters
    ----------
    input_file : string
        Path to input file.

    Returns
    -------
    input_dict : dictionary
        Dictionary storing WarpX input parameters
        (parameter's name stored as key, parameter's value stored as value).
    """
    input_dict = dict()
    with open(input_file) as ff:
        for line in ff:
            sline = line.split('=')
            # skip lines that are commented out, blank, or continuation of previous parameters
            skip_line = sline[0].startswith('#') or sline[0].startswith('\n') or len(sline) == 1
            if not skip_line:
                key = sline[0].strip()
                val = sline[1].split()
                # The value corresponding to a given key of input_dict is a list
                # of strings, from which we remove any leftover comments
                for i in range(len(val)):
                    if val[i].startswith('#'):
                        val = val[:i]
                        break
                input_dict[key] = val
    return input_dict

input_dict = parse_input_file('inputs_3d_multiple_particles')

uex1, uex2, uex3 = [float(w) for w in input_dict['beam_e.multiple_particles_ux']]
uey1, uey2, uey3 = [float(w) for w in input_dict['beam_e.multiple_particles_uy']]
uez1, uez2, uez3 = [float(w) for w in input_dict['beam_e.multiple_particles_uz']]

xe = np.array([float(w) for w in input_dict['beam_e.multiple_particles_pos_x']])
ye = np.array([float(w) for w in input_dict['beam_e.multiple_particles_pos_y']])

we = np.array([float(w) for w in input_dict['beam_e.multiple_particles_weight']])

upx1, upx2, upx3 = [float(w) for w in input_dict['beam_p.multiple_particles_ux']]
upy1, upy2, upy3 = [float(w) for w in input_dict['beam_p.multiple_particles_uy']]
upz1, upz2, upz3 = [float(w) for w in input_dict['beam_p.multiple_particles_uz']]

xp1, xp2, xp3 = [float(w) for w in input_dict['beam_p.multiple_particles_pos_x']]
yp1, yp2, yp3 = [float(w) for w in input_dict['beam_p.multiple_particles_pos_y']]
zp1, zp2, zp3 = [float(w) for w in input_dict['beam_p.multiple_particles_pos_z']]

Ex, Ey, Ez = [float(w) for w in input_dict['particles.E_external_particle']]
Bx, By, Bz = [float(w) for w in input_dict['particles.B_external_particle']]

pex1, pex2, pex3 = uex1*m_e*c, uex2*m_e*c, uex3*m_e*c
pey1, pey2, pey3 = uey1*m_e*c, uey2*m_e*c, uey3*m_e*c
pez1, pez2, pez3 = uez1*m_e*c, uez2*m_e*c, uez3*m_e*c

ppx1, ppx2, ppx3 = upx1*m_e*c, upx2*m_e*c, upx3*m_e*c
ppy1, ppy2, ppy3 = upy1*m_e*c, upy2*m_e*c, upy3*m_e*c
ppz1, ppz2, ppz3 = upz1*m_e*c, upz2*m_e*c, upz3*m_e*c

ge1, ge2, ge3 = np.sqrt(1.+uex1**2+uey1**2+uez1**2), np.sqrt(1.+uex2**2+uey2**2+uez2**2), np.sqrt(1.+uex3**2+uey3**2+uez3**2)
gp1, gp2, gp3 = np.sqrt(1.+upx1**2+upy1**2+upz1**2), np.sqrt(1.+upx2**2+upy2**2+upz2**2), np.sqrt(1.+upx3**2+upy3**2+upz3**2)

Ee1, Ee2, Ee3 = m_e*c**2*(ge1-1.), m_e*c**2*(ge2-1.), m_e*c**2*(ge3-1.)
Ep1, Ep2, Ep3 = m_e*c**2*(gp1-1.), m_e*c**2*(gp2-1.), m_e*c**2*(gp3-1.)


def chi(ux, uy, uz, Ex=Ex, Ey=Ey, Ez=Ez, Bx=Bx, By=By, Bz=Bz):
    #print(ux, uy, uz, Ex, Ey, Ez, Bx, By, Bz)
    gamma = np.sqrt(1.+ux**2+uy**2+uz**2)
    vx = ux / gamma * c
    vy = uy / gamma * c
    vz = uz / gamma * c
    tmp1x = Ex + vy*Bz - vz*By
    tmp1y = Ey - vx*Bz + vz*Bx
    tmp1z = Ez + vx*By - vy*Bx
    tmp2 = (Ex*vx + Ey*vy + Ez*vz)/c

    chi = gamma/E_crit*np.sqrt(tmp1x**2+tmp1y**2+tmp1z**2 - tmp2**2)
    return chi


def dL_dt():
    series = io.Series("diags/diag1/openpmd_%T.bp",io.Access.read_only)
    iterations = np.asarray(series.iterations)
    lumi = []
    for n,ts in enumerate(iterations):
        it = series.iterations[ts]
        rho1 = it.meshes["rho_beam_e"]
        dV = np.prod(rho1.grid_spacing)
        rho1 = it.meshes["rho_beam_e"][io.Mesh_Record_Component.SCALAR].load_chunk()
        rho2 = it.meshes["rho_beam_p"][io.Mesh_Record_Component.SCALAR].load_chunk()
        q1 = np.unique(it.particles["beam_e"]["charge"][io.Mesh_Record_Component.SCALAR].load_chunk())
        q2 = np.unique(it.particles["beam_p"]["charge"][io.Mesh_Record_Component.SCALAR].load_chunk())
        series.flush()
        n1 = rho1/q1
        n2 = rho2/q2
        #print(np.where((n1>0), n1, 0.))
        #print(n1)
        l = 2*np.sum(n1*n2)*dV*c
        lumi.append(l)
        print('--------------')
        print(np.sum(n1*n2), np.sum(n1*n1), np.sum(n2*n2))
    return lumi

CHI_ANALYTICAL = np.array([chi(uex1, uey1, uez1), chi(uex2, uey2, uez2), chi(uex3, uey3, uez3)])
THETAX = np.array([np.arctan2(uex1, uez1), np.arctan2(uex2, uez2), np.arctan2(uex3, uez3)])
THETAY = np.array([np.arctan2(uey1, uez1), np.arctan2(uey2, uez2), np.arctan2(uey3, uez3)])

CollDiagFname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
df = pd.read_csv(CollDiagFname, sep=" ", header=0)

# CHI MAX
fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimax_pe = np.loadtxt(fname)[:,19]
chimax_cr = df[[col for col in df.columns if 'chimax_beam_e' in col]].to_numpy()
assert np.allclose(np.max(CHI_ANALYTICAL), chimax_cr, rtol=1e-8)
assert np.allclose(chimax_pe, chimax_cr, rtol=1e-8)

# CHI MIN
fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimin_pe = np.loadtxt(fname)[:,18]
chimin_cr = df[[col for col in df.columns if 'chimin_beam_e' in col]].to_numpy()
assert np.allclose(np.min(CHI_ANALYTICAL), chimin_cr, rtol=1e-8)
assert np.allclose(chimin_pe, chimin_cr, rtol=1e-8)

# CHI AVERAGE
chiave_cr = df[[col for col in df.columns if 'chiave_beam_e' in col]].to_numpy()
assert np.allclose(np.average(CHI_ANALYTICAL, weights=we), chiave_cr, rtol=1e-8)

# X AVE STD
x_ave_cr = df[[col for col in df.columns if 'xave_beam_e' in col]].to_numpy()
x_std_cr = df[[col for col in df.columns if 'xstd_beam_e' in col]].to_numpy()
x_ave = np.average(xe, weights=we)
x_std = np.sqrt(np.average((xe-x_ave)**2, weights=we))
assert np.allclose(x_ave, x_ave_cr, rtol=1e-8)
assert np.allclose(x_std, x_std_cr, rtol=1e-8)

# Y AVE STD
y_ave_cr = df[[col for col in df.columns if 'yave_beam_e' in col]].to_numpy()
y_std_cr = df[[col for col in df.columns if 'ystd_beam_e' in col]].to_numpy()
y_ave = np.average(ye, weights=we)
y_std = np.sqrt(np.average((ye-y_ave)**2, weights=we))
assert np.allclose(y_ave, y_ave_cr, rtol=1e-8)
assert np.allclose(y_std, y_std_cr, rtol=1e-8)


# THETA X MIN AVE MAX STD
thetax_min_cr = df[[col for col in df.columns if 'thetax_min_beam_e' in col]].to_numpy()
thetax_ave_cr = df[[col for col in df.columns if 'thetax_ave_beam_e' in col]].to_numpy()
thetax_max_cr = df[[col for col in df.columns if 'thetax_max_beam_e' in col]].to_numpy()
thetax_std_cr = df[[col for col in df.columns if 'thetax_std_beam_e' in col]].to_numpy()
thetax_min = np.min(THETAX)
thetax_ave = np.average(THETAX, weights=we)
thetax_max = np.max(THETAX)
thetax_std = np.sqrt(np.average((THETAX-thetax_ave)**2, weights=we))
assert np.allclose(thetax_min, thetax_min_cr, rtol=1e-8)
assert np.allclose(thetax_ave, thetax_ave_cr, rtol=1e-8)
assert np.allclose(thetax_max, thetax_max_cr, rtol=1e-8)
assert np.allclose(thetax_std, thetax_std_cr, rtol=1e-8)

# THETA Y MIN AVE MAX STD
thetay_min_cr = df[[col for col in df.columns if 'thetay_min_beam_e' in col]].to_numpy()
thetay_ave_cr = df[[col for col in df.columns if 'thetay_ave_beam_e' in col]].to_numpy()
thetay_max_cr = df[[col for col in df.columns if 'thetay_max_beam_e' in col]].to_numpy()
thetay_std_cr = df[[col for col in df.columns if 'thetay_std_beam_e' in col]].to_numpy()
thetay_min = np.min(THETAY)
thetay_ave = np.average(THETAY, weights=we)
thetay_max = np.max(THETAY)
thetay_std = np.sqrt(np.average((THETAY-thetay_ave)**2, weights=we))
assert np.allclose(thetay_min, thetay_min_cr, rtol=1e-8)
assert np.allclose(thetay_ave, thetay_ave_cr, rtol=1e-8)
assert np.allclose(thetay_max, thetay_max_cr, rtol=1e-8)
assert np.allclose(thetay_std, thetay_std_cr, rtol=1e-8)


# dL/dt
dL_dt_cr = df[[col for col in df.columns if 'dL_dt' in col]].to_numpy()
assert np.allclose(dL_dt_cr, dL_dt(), rtol=1e-8)
