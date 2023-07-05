import re

from matplotlib import cm, use
import matplotlib.colors
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import openpmd_api as io
from scipy.constants import c, e, hbar, m_e

E_crit = m_e**2*c**3/(e*hbar)
B_crit = m_e**2*c**2/(e*hbar)

# extract numbers from a string
def find_num_in_line(line):
    items = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', line)
    fitems = [float(it) for it in items]
    if len(fitems)==1:
        return fitems[0]
    else:
        return fitems

# get input parameters from warpx_used_inputs
with open('./warpx_used_inputs', 'rt') as f:
    lines = f.readlines()
    for line in lines:
        if 'my_constants.nx' in line:
            nx = find_num_in_line(line)
        if 'my_constants.ny' in line:
            ny = find_num_in_line(line)
        if 'my_constants.nz' in line:
            nz = find_num_in_line(line)
        if 'my_constants.Lx' in line:
            Lx = find_num_in_line(line)
        if 'my_constants.Ly' in line:
            Ly = find_num_in_line(line)
        if 'my_constants.Lz' in line:
            Lz = find_num_in_line(line)
        if 'beam_e.density' in line:
            n1 = find_num_in_line(line)
        if 'beam_p.density' in line:
            n2 = find_num_in_line(line)
        if 'beam_e.ux' in line:
            u1x = find_num_in_line(line)
        if 'beam_e.uy' in line:
            u1y = find_num_in_line(line)
        if 'beam_e.uz' in line:
            u1z = find_num_in_line(line)
        if 'beam_p.ux' in line:
            u2x = find_num_in_line(line)
        if 'beam_p.uy' in line:
            u2y = find_num_in_line(line)
        if 'beam_p.uz' in line:
            u2z = find_num_in_line(line)
        if 'beam_e.num_particles_per_cell_each_dim' in line:
            aux1, aux2, aux3 = find_num_in_line(line)
            Ne = aux1*aux2*aux3
        if 'beam_p.num_particles_per_cell_each_dim' in line:
            aux1, aux2, aux3 = find_num_in_line(line)
            Np = aux1*aux2*aux3
        if 'particles.E_external_particle' in line:
            Ex, Ey, Ez = find_num_in_line(line)
        if 'particles.B_external_particle' in line:
            Bx, By, Bz = find_num_in_line(line)

dV = Lx/nx * Ly/ny * Lz/nz
xmin = -0.5*Lx
ymin = -0.5*Ly
zmin = -0.5*Lz
xmax = 0.5*Lx
ymax = 0.5*Ly
midx = 0.5*(xmax - xmin)
midy = 0.5*(ymax - ymin)

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

def luminosity():
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
        #print('llllllllllllllllllllll ', l, np.sum(n1*n2), np.sum(n1*n1),np.sum(n2*n2))
        lumi.append(l)
    return lumi

def disruption():
    series = io.Series("diags/diag1/openpmd_%T.bp",io.Access.read_only)
    iterations = np.asarray(series.iterations)
    XY1_AVE, XY1_STD, XY2_AVE, XY2_STD = [], [], [], []
    for n,ts in enumerate(iterations):
        it = series.iterations[ts]

        y1 = it.particles["beam_e"]["position"]["y"].load_chunk()
        y2 = it.particles["beam_p"]["position"]["y"].load_chunk()

        x1 = it.particles["beam_e"]["position"]["x"].load_chunk()
        x2 = it.particles["beam_p"]["position"]["x"].load_chunk()

        w1 = it.particles["beam_e"]["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk()
        w2 = it.particles["beam_p"]["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk()

        series.flush()

        xy1 = np.sqrt((x1-midx)**2+(y1-midy)**2)
        xy2 = np.sqrt((x2-midx)**2+(y2-midy)**2)


        xy1_ave = np.average(xy1, weights=w1)
        xy2_ave = np.average(xy2, weights=w2)

        xy1_std = np.sqrt(np.average((xy1-xy1_ave)**2, weights=w1))
        xy2_std = np.sqrt(np.average((xy2-xy2_ave)**2, weights=w2))

        XY1_AVE.append(xy1_ave)
        XY1_STD.append(xy1_std)

        XY2_AVE.append(xy2_ave)
        XY2_STD.append(xy2_std)

    return np.asarray([XY1_AVE, XY1_STD, XY2_AVE, XY2_STD])

print('chi of electrons ----------------------------------------')
print('theory:', chi(u1x, u1y, u1z))

fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimax = np.loadtxt(fname)[:,19]
print('chimax ParticleExtrema diag = ', chimax)

fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
chimax = np.loadtxt(fname)[:,10]
print('chimax ColliderRelevant diag = ', chimax)

fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimin = np.loadtxt(fname)[:,18]
print('chimin ParticleExtrema diag = ', chimin)

fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
chimin = np.loadtxt(fname)[:,8]
print('chimin ColliderRelevant diag = ', chimin)

fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
chiave = np.loadtxt(fname)[:,9]
print('chiave ColliderRelevant diag = ', chiave)

print('-------------------------')
print('luminosity')

print('from PIC data', luminosity())

CollDiagFname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
lum = np.loadtxt(CollDiagFname)[:,2]
print('ColliderRelevant diag = ', lum)

print('theory from input', 2.*n1*n2*Lx*Ly*Lz*c)

print('-------------------------')
print('xy1 ave')
print('from PIC data', disruption()[0])
print('ColliderRelevant diag ele = ',  np.loadtxt(CollDiagFname)[:,11])

print('-------------------------')
print('xy1 std')
print('from PIC data', disruption()[1])
print('ColliderRelevant diag ele = ',  np.loadtxt(CollDiagFname)[:,12])

print('disruption')
print('from PIC data', disruption()[1]/disruption()[0])
print('ColliderRelevant diag ele = ',  np.loadtxt(CollDiagFname)[:,12]/np.loadtxt(CollDiagFname)[:,11])
