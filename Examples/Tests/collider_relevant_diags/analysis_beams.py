
from matplotlib import cm, use
import matplotlib.colors
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from numpy import exp, sqrt
import openpmd_api as io
from scipy.constants import alpha, c
from scipy.constants import e as q_e
from scipy.constants import hbar, m_e, physical_constants, pi

r_e = physical_constants["classical electron radius"][0]


E_crit = m_e**2*c**3/(q_e*hbar)
B_crit = m_e**2*c**2/(q_e*hbar)



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

input_dict = parse_input_file('inputs_3d_beams')
sigmax, = [float(x) for x in input_dict['my_constants.sigmax']]
sigmay, = [float(x) for x in input_dict['my_constants.sigmay']]
sigmaz, = [float(x) for x in input_dict['my_constants.sigmaz']]
N, = [float(x) for x in input_dict['my_constants.beam_npart']]
gamma, = [float(x) for x in input_dict['my_constants.gammab']]
charge = N * q_e
n0 = charge / (q_e * sigmax * sigmay * sigmaz * (2.*pi)**(3./2.))
Lx, Ly, Lz = 7*sigmax, 7*sigmay, 14*sigmaz
nx, = [float(x) for x in input_dict['my_constants.nx']]
ny, = [float(x) for x in input_dict['my_constants.ny']]
nz, = [float(x) for x in input_dict['my_constants.nz']]

def dL_dt_full():
    series = io.Series("diags/diag1/openpmd_%T.bp",io.Access.read_only)
    iterations = np.asarray(series.iterations)
    lumi = []
    for n,ts in enumerate(iterations):
        it = series.iterations[ts]
        rho1 = it.meshes["rho_beam_e"]
        dV = np.prod(rho1.grid_spacing)
        rho1 = it.meshes["rho_beam_e"][io.Mesh_Record_Component.SCALAR].load_chunk()
        rho2 = it.meshes["rho_beam_p"][io.Mesh_Record_Component.SCALAR].load_chunk()
        series.flush()
        n1 = -rho1/q_e
        n2 =  rho2/q_e
        l = 2*np.sum(n1*n2)*dV*c
        lumi.append(l)
    return lumi

def num_dens(x,y,z):
    return n0 * exp(-x**2/(2*sigmax**2))*exp(-y**2/(2*sigmay**2))*exp(-z**2/(2*sigmaz**2))




'''
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
chiave = np.loadtxt(fname)[:,9]
print('chiave ColliderRelevant diag = ', chiave)

print('-------------------------')
print('luminosity')

print('from PIC data', luminosity())



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
'''

CollDiagFname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'

##################
### LUMINOSITY ###
##################
dL_dt_cr = np.loadtxt(CollDiagFname)[:,2]
times = np.loadtxt(CollDiagFname)[:,1]
L_cr = np.trapz(dL_dt_cr, times)
coll_timestep = np.argmax(dL_dt_cr)



x = np.linspace(-3*sigmax,3*sigmax,128)
y = np.linspace(-3*sigmay,3*sigmay,128)
z = np.linspace(-3*sigmaz,3*sigmaz,128)
dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]
X,Y,Z = np.meshgrid(x,y,z)
D = num_dens(X,Y,Z)
dL_dt = 2.*np.sum(D**2)*dx*dy*dz*c

# see formula (2.21) in Yokoya and Chen
L_theory = N**2 / (4*pi*sigmax*sigmay)

#assert np.isclose(L_theory, L_cr, rtol=0.1, atol=1e-12)
#assert np.isclose(dL_dt, np.max(dL_dt_cr), rtol=1e-1, atol=1e-12)
#assert np.allclose(dL_dt_full(), dL_dt_cr, rtol=1e-12, atol=1e-12)

###############
### CHI MIN ###
###############
fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimin_ele_pe = np.loadtxt(fname)[:,18]

fname='diags/reducedfiles/ParticleExtrema_beam_p.txt'
chimin_pos_pe = np.loadtxt(fname)[:,18]

chimin_ele_cr = np.loadtxt(CollDiagFname)[:,8]
chimin_pos_cr = np.loadtxt(CollDiagFname)[:,3]

assert np.allclose(chimin_pos_pe, chimin_pos_cr, rtol=1e-12, atol=1e-12)
assert np.allclose(chimin_ele_pe, chimin_ele_cr, rtol=1e-12, atol=1e-12)

###############
### CHI MAX ###
###############
# see formula (3.52) in Yokoya and Chen
chimax_theory = 2*N*r_e**2*gamma/(alpha*sigmaz*(sigmax+sigmay))

fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimax_ele_pe = np.loadtxt(fname)[:,19]

fname='diags/reducedfiles/ParticleExtrema_beam_p.txt'
chimax_pos_pe = np.loadtxt(fname)[:,19]

chimax_ele_cr = np.loadtxt(CollDiagFname)[:,10]
chimax_pos_cr = np.loadtxt(CollDiagFname)[:,5]

assert np.allclose(chimin_ele_pe, chimin_ele_cr, rtol=1e-12, atol=1e-12)
assert np.allclose(chimin_pos_pe, chimin_pos_cr, rtol=1e-12, atol=1e-12)

plt.plot(chimax_theory*np.ones_like(chimax_ele_cr))
plt.plot(chimax_ele_cr)
plt.plot(chimax_pos_cr)
plt.show()

###############
### CHI AVE ###
###############
# see formula (3.52) in Yokoya and Chen
chiave_theory = 5./12.*chimax_theory

chiave_ele_cr = np.loadtxt(CollDiagFname)[:,9]
chiave_pos_cr = np.loadtxt(CollDiagFname)[:,4]

#plt.plot(chiave_theory*np.ones_like(chiave_pos_cr))
plt.plot(chiave_pos_cr)
plt.plot(chiave_ele_cr)
plt.show()
print(chiave_pos_cr)
print(chiave_theory)

##################
### LUMINOSITY ###
##################
dL_dt_cr = np.loadtxt(CollDiagFname)[:,2]
times = np.loadtxt(CollDiagFname)[:,1]
L_cr = np.trapz(dL_dt_cr, times)
coll_timestep = np.argmax(dL_dt_cr)



x = np.linspace(-3*sigmax,3*sigmax,128)
y = np.linspace(-3*sigmay,3*sigmay,128)
z = np.linspace(-3*sigmaz,3*sigmaz,128)
dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]
X,Y,Z = np.meshgrid(x,y,z)
D = num_dens(X,Y,Z)
dL_dt = 2.*np.sum(D**2)*dx*dy*dz*c

# see formula (2.21) in Yokoya and Chen
L_theory = N**2 / (4*pi*sigmax*sigmay)

assert np.isclose(L_theory, L_cr, rtol=0.1, atol=1e-12)
assert np.isclose(dL_dt, np.max(dL_dt_cr), rtol=1e-1, atol=1e-12)
assert np.allclose(dL_dt_full(), dL_dt_cr, rtol=1e-12, atol=1e-12)


##################
### DISRUPTION ###
##################
x_std_ele = np.loadtxt(CollDiagFname)[:,11]
y_std_ele = np.loadtxt(CollDiagFname)[:,12]
x_std_pos = np.loadtxt(CollDiagFname)[:,6]
y_std_pos = np.loadtxt(CollDiagFname)[:,7]

#see formula (2.13) from Yokoya and Chen
Dx = 2*N*r_e*sigmaz/(gamma*sigmax*(sigmax+sigmay))
Dy = 2*N*r_e*sigmaz/(gamma*sigmay*(sigmax+sigmay))

boh = (x_std_ele[coll_timestep]-x_std_ele[0])/x_std_ele[0]
boh1 = (y_std_ele[coll_timestep]-y_std_ele[0])/y_std_ele[0]

print(Dx, Dy, boh)




plt.plot(Dx*np.ones_like(x_std_pos))
print(Dx)
plt.plot(np.abs(x_std_ele-x_std_ele[0])/x_std_ele[0])
plt.plot(np.abs(x_std_pos-x_std_pos[0])/x_std_pos[0])


plt.plot(Dy*np.ones_like(x_std_pos))
print(Dy)
plt.plot(np.abs(y_std_ele-y_std_ele[0])/y_std_ele[0])
plt.plot(np.abs(y_std_pos-y_std_pos[0])/y_std_pos[0])

plt.axvline(x=coll_timestep)
#plt.show()

# see formula (2.35) in Yokoya and Chen
#theta = 2*N*r_e/(gamma*(sigmax+sigmay))

#plt.plot( 2.*np.sum(D**2)*dx*dy*dz*c*np.ones_like(chiave_pos_cr))
#plt.plot(dL_dt_cr)
#plt.plot(dL_dt_full())
#plt.show()
#plt.close()


fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(12,5), dpi=300)

ax[0].plot(times, chiave_pos_cr)
ax[0].plot(times, chiave_ele_cr)
ax[0].plot(times, chiave_theory*np.ones_like(chiave_pos_cr))


ax[1].plot(times, chimax_pos_cr)
ax[1].plot(times, chimax_ele_cr)
ax[1].plot(times, chimax_theory*np.ones_like(chiave_pos_cr))


ax[2].plot(times, dL_dt_cr)
ax[2].plot(times, dL_dt*np.ones_like(dL_dt_cr))
ax[2].plot(times, dL_dt_full())

plt.show()
