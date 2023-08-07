
from matplotlib import cm, use
import matplotlib.colors
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from numpy import exp, sqrt
import openpmd_api as io
import pandas as pd
from scipy.constants import alpha, c
from scipy.constants import e as q_e
from scipy.constants import hbar, m_e, physical_constants, pi

r_e = physical_constants["classical electron radius"][0]


E_crit = m_e**2*c**3/(q_e*hbar)
B_crit = m_e**2*c**2/(q_e*hbar)

def rerr(a,b):
    return np.abs(a-b)/np.abs(np.minimum(a,b))

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
df = pd.read_csv(CollDiagFname, sep=" ", header=0)



##################
### LUMINOSITY ###
##################
dL_dt_cr = df[[col for col in df.columns if 'dL_dt' in col]].to_numpy()
times = df[[col for col in df.columns if 'time' in col]].to_numpy()
L_cr = np.trapz(dL_dt_cr, times)
coll_timestep = np.argmax(dL_dt_cr)
coll_time = times[coll_timestep]

x = np.linspace(-3*sigmax,3*sigmax,128)
y = np.linspace(-3*sigmay,3*sigmay,128)
z = np.linspace(-3*sigmaz,3*sigmaz,128)
dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]
X,Y,Z = np.meshgrid(x,y,z)
D = num_dens(X,Y,Z)
dL_dt_theory = 2.*np.sum(D**2)*dx*dy*dz*c

# see formula (2.21) in Yokoya and Chen
L_theory = N**2 / (4*pi*sigmax*sigmay)

# from guinea-pig: lumi_fine
L_gp = 7.95485e+30

assert np.isclose(L_cr, L_theory, rtol=0.02, atol=0.)
assert np.isclose(L_cr, L_gp, rtol=0.003, atol=0,)
assert np.isclose(np.max(dL_dt_cr), dL_dt_theory, rtol=0.003, atol=0,)
assert np.allclose(dL_dt_full(), dL_dt_cr, rtol=1e-12, atol=0.)

###############
### CHI MIN ###
###############
fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimin_ele_pe = np.loadtxt(fname)[:,18]

fname='diags/reducedfiles/ParticleExtrema_beam_p.txt'
chimin_pos_pe = np.loadtxt(fname)[:,18]

chimin_ele_cr = np.loadtxt(CollDiagFname)[:,14]
chimin_pos_cr = np.loadtxt(CollDiagFname)[:,3]

assert np.allclose(chimin_pos_pe, chimin_pos_cr, rtol=1e-16, atol=0.)
assert np.allclose(chimin_ele_pe, chimin_ele_cr, rtol=1e-16, atol=0.)
assert np.allclose(chimin_ele_cr, chimin_ele_cr, rtol=1e-12, atol=0.)

###############
### CHI MAX ###
###############
# see formula (3.52) in Yokoya and Chen
chimax_theory = 2*N*r_e**2*gamma/(alpha*sigmaz*(sigmax+sigmay))

# from guinea-pig: upsmax (maximal value of the beamstrahlung paramater that occured)
chimax_gp = 1.51733

fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimax_ele_pe = np.loadtxt(fname)[:,19]

fname='diags/reducedfiles/ParticleExtrema_beam_p.txt'
chimax_pos_pe = np.loadtxt(fname)[:,19]

chimax_ele_cr = df[[col for col in df.columns if 'chimax_beam_e' in col]].to_numpy()
chimax_pos_cr = df[[col for col in df.columns if 'chimax_beam_p' in col]].to_numpy()

assert np.allclose(chimax_pos_cr, chimax_pos_pe, rtol=1e-16, atol=0.)
assert np.allclose(chimax_ele_cr, chimax_ele_pe, rtol=1e-16, atol=0.)
assert np.isclose(chimax_ele_pe[coll_timestep], chimax_theory, rtol=0.2, atol=0.)
assert np.isclose(chimax_pos_pe[coll_timestep], chimax_theory, rtol=0.2, atol=0.)

###############
### CHI AVE ###
###############
# see formula (3.52) in Yokoya and Chen
chiave_theory = 5./12.*chimax_theory

chiave_ele_cr = df[[col for col in df.columns if 'chiave_beam_e' in col]].to_numpy()
chiave_pos_cr = df[[col for col in df.columns if 'chiave_beam_p' in col]].to_numpy()

assert np.isclose(chiave_ele_cr[coll_timestep], chiave_theory, rtol=0.2, atol=0.)
assert np.isclose(chiave_pos_cr[coll_timestep], chiave_theory, rtol=0.2, atol=0.)

##################
### DISRUPTION ###
##################
xave_ele = df[[col for col in df.columns if ']x_ave_beam_e' in col]].to_numpy()
yave_ele = df[[col for col in df.columns if ']y_ave_beam_e' in col]].to_numpy()
xstd_ele = df[[col for col in df.columns if ']x_std_beam_e' in col]].to_numpy()
ystd_ele = df[[col for col in df.columns if ']y_std_beam_e' in col]].to_numpy()

thetaxave_ele = df[[col for col in df.columns if 'thetax_ave_beam_e' in col]].to_numpy()
thetayave_ele = df[[col for col in df.columns if 'thetay_ave_beam_e' in col]].to_numpy()
thetaxstd_ele = df[[col for col in df.columns if 'thetax_std_beam_e' in col]].to_numpy()
thetaystd_ele = df[[col for col in df.columns if 'thetay_std_beam_e' in col]].to_numpy()

thetaxmax_ele = df[[col for col in df.columns if 'thetax_max_beam_e' in col]].to_numpy()
thetaymax_ele = df[[col for col in df.columns if 'thetay_max_beam_e' in col]].to_numpy()

ref_time = times[coll_timestep] + sqrt(3)*sigmaz/c
ref_timestep = np.argmin(np.abs(times - ref_time)**2)


theta0 = 2*N*r_e/(gamma*(sigmax+sigmay))
print('theta! = ', theta0)

#see formula (2.13) from Yokoya and Chen
Dx = 2*N*r_e*sigmaz/(gamma*sigmax*(sigmax+sigmay))
Dy = 2*N*r_e*sigmaz/(gamma*sigmay*(sigmax+sigmay))

print('Dx! = ', Dx)


thetax_ave_ele_gp = 0.201504 # microradians
thetay_ave_ele_gp = 0.0879509 # microradians
thetax_std_ele_gp = 773.377 # microradians
thetay_std_ele_gp = 773.696 # microradians


##################
#### PHOTONS #####
##################

chi = np.max(chiave_ele_cr)
n_gamma = 5./2. * alpha**2 * sigmaz / r_e / gamma * chi / sqrt(1.+chi**(2./3.))
print(chi, n_gamma)


fname='diags/reducedfiles/ParticleNumber.txt'
N_pho = np.loadtxt(fname)[:,9]
N_ele = np.loadtxt(fname)[:,8]
#fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5,5), dpi=300)
#ax.plot(times, N_pho/N_ele[0])
print(N_pho[-1]/N_ele[0])
#boh = np.abs(x_std_ele-x_std_ele[0])/x_std_ele[0]
#boh1 = np.abs(y_std_ele-y_std_ele[0])/y_std_ele[0]

#print(Dx, boh[coll_timestep])
# see formula (2.35) in Yokoya and Chen
#theta = 2*N*r_e/(gamma*(sigmax+sigmay))

#plt.plot( 2.*np.sum(D**2)*dx*dy*dz*c*np.ones_like(chiave_pos_cr))
#plt.plot(dL_dt_cr)
#plt.plot(dL_dt_full())
#plt.show()
#plt.close()


fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15,5), dpi=300)

ax[0].plot(times, chiave_pos_cr, label='pos cr')
ax[0].plot(times, chiave_ele_cr, label='ele cr')
ax[0].plot(times, chiave_theory*np.ones_like(chiave_pos_cr), label='theory')
ax[0].set_title('chi ave')

ax[3].plot(times, chimin_pos_pe, label='pos pe')
ax[3].plot(times, chimin_pos_cr, label='pos cr')
ax[3].plot(times, chimin_ele_pe, label='ele pe')
ax[3].plot(times, chimin_ele_cr, label='ele cr')
ax[3].set_title('chi min')

ax[1].plot(times, chimax_pos_pe, label='pos pe')
ax[1].plot(times, chimax_ele_pe, label='ele pe')
ax[1].plot(times, chimax_pos_cr, label='pos cr')
ax[1].plot(times, chimax_ele_cr, label='ele cr')
ax[1].plot(times, chimax_theory*np.ones_like(chiave_pos_cr), label='theory')
ax[1].set_title('chi max')

ax[2].plot(times, dL_dt_cr, label='cr')
ax[2].plot(times, dL_dt_theory*np.ones_like(dL_dt_cr), label='theory')
ax[2].plot(times, dL_dt_full(), label='full')
ax[2].set_title('dL/dt')

for a in ax.reshape(-1):
    a.legend()
    a.axvline(x=times[coll_timestep], color='black')
    a.axvline(x=times[ref_timestep], color='magenta')

fig.savefig('x.png', dpi=300.)
plt.close()



fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(15,5), dpi=300)

ax[0].plot(times, xstd_ele, label='x')
ax[0].plot(times, ystd_ele, label='y')
ax[0].set_title('pos std')

ax[1].plot(times, thetaxstd_ele, label='x')
ax[1].plot(times, thetaystd_ele, label='y')
ax[1].set_title('theta std')

ax[2].plot(times, thetaxmax_ele, label='x')
ax[2].plot(times, thetaymax_ele, label='y')
ax[2].set_title('theta max')


for a in ax.reshape(-1):
    a.legend()
    a.axvline(x=times[coll_timestep], color='black')
    a.axvline(x=times[ref_timestep], color='magenta')

fig.savefig('y.png', dpi=300.)
plt.close()




fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(15,5), dpi=300)
b = xstd_ele[:,0]
t = times[:,0]


db_dt = np.diff(b)/np.diff(t)

ax[0].plot(t[:-1], db_dt/c)
ax[0].plot(times, -Dx*np.ones_like(times))
ax[0].plot(t[:-11], db_dt[:-10]/c)


ax[1].plot(times, thetaxstd_ele, label='x')
ax[1].plot(times, thetaystd_ele, label='y')
ax[1].plot(times, theta0*np.ones_like(times))
ax[1].plot(times, thetax_std_ele_gp*np.ones_like(times)/1e6)
ax[1].plot(times, thetay_std_ele_gp*np.ones_like(times)/1e6)
ax[1].set_title('theta std')

ax[2].plot(times, thetaxave_ele, label='x')
ax[2].plot(times, thetayave_ele, label='y')
#ax[2].plot(times, theta0*np.ones_like(times))
#ax[2].plot(times, thetax_ave_ele_gp*np.ones_like(times)/1e6)
#ax[2].plot(times, thetay_ave_ele_gp*np.ones_like(times)/1e6)
ax[2].set_title('theta ave')


for a in ax.reshape(-1):
    a.legend()
    a.axvline(x=times[coll_timestep], color='black')
    a.axvline(x=times[ref_timestep], color='magenta')

fig.savefig('z.png', dpi=300.)
plt.show()



plt.close()

#ax[0].plot(times, -Dx*np.ones_like(times), label='Dx theory')
#ax[0].plot(times, -Dx*c*(times-coll_time)/sigmaz, label='dx/x theory')


#fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6,5), dpi=300)
#dxstd = (xstd_ele[:,0]/ xstd_ele[0,0])


#ax.plot(times, np.diff(dxstd))
#ax.plot(times, np.zeros_like(times))



fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10,5), dpi=300)
b = xstd_ele[:,0]
t = times[:,0]


db_dt = np.diff(b)/np.diff(t)

ax[0].plot(t[:-1], db_dt/c)
ax[0].plot(times, -Dx*np.ones_like(times))
ax[0].plot(t[121:], db_dt[120:]/c, lw=5)


m = np.average(db_dt[120:]/c)

T = - b[0] / m
tt = np.linspace(0, T, 500)

ax[1].plot(tt, b[0]+m*tt)
ax[1].plot(tt, np.zeros_like(tt))

myDx = sigmax / (T * c)

print(myDx)


plt.show()
#ax[1].plot(times, -(0.78+0.2*Dx)*np.ones_like(times))
#ax[1].plot(times, -(r_e*N/gamma/sigmaz**2)*np.ones_like(times)* xstd_ele[0])

#*c/sigmaz, label='dx/x theory')

#print(-Dx*c/sigmaz)
#print(db_dt[ref_timestep])
#ax[1].plot(times[:-1], -r_e*N/(gamma*sigmaz**2)*np.ones_like(times[:-1]), label='factor')
#ax[1].plot(times[:-1], par, label='bho')


















#ax[0].plot(times, theta0*np.ones_like(times), label='theory', color='red', zorder=11)

#ax[0].plot(times, thetaxave_ele, label='y')
#ax[0].set_title('theta ave')

#ax[0].plot(times, xstd_ele)


#ax[0].plot(times, (thetaxstd_ele-thetaxstd_ele[0])/thetaxstd_ele[0], label='thetax')


#ax[0].plot(times, (thetaystd_ele-thetaystd_ele[0])/thetaystd_ele[0], label='thetay')
#ax[0].plot(times, (ystd_ele - ystd_ele[0])/ystd_ele[0], label='y')


#ax[0].plot(times, thetaxstd_ele, label='y')
#ax[0].set_title('stds')

#ax[0].plot(times, -Dx*c*(times-coll_time)/sigmaz, label='dx/x theory')


#ax[2].plot(times, xave_ele, label='x ele')
#ax[2].plot(times, yave_ele, label='y ele')
#ax[2].set_title('ave pos')

#ax[1].plot(times, xstd_ele, label='x ele')
#ax[1].plot(times, ystd_ele, label='y ele')
#ax[1].set_title('std pos')



#ax[1].plot(times, -Dx*np.ones_like(times), label='Dx theory')
#ax[1].plot(times, -Dx*c*(times-coll_time)/sigmaz, label='dx/x theory')



#ax[2].plot(times, xstd_ele, label='std')

#ax[3].fill_between(times, yave_ele-ystd_ele, yave_ele+ystd_ele, alpha=0.4)
#ax[3].plot(times, yave_ele, label='ave')
#ax[3].plot(times, ystd_ele, label='std')
#ax[3].set_title('y')
