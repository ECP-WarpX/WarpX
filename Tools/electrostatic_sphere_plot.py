import yt
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from skimage import measure
import numpy as np

yt.funcs.mylog.setLevel(50)    #Mute yt's verbose shenanigans
plotfilepath = '../Examples/Tests/ElectrostaticSphere/diags/plotfiles'

#Read header file for data to extract grid and cell data
with open(plotfilepath + '/plt00000/Header','r') as file:
    header = file.readlines()
[xmin, ymin, zmin] = [eval(header[9].strip().split()[i]) for i in range(3)]
[xmax, ymax, zmax] = [eval(header[10].strip().split()[i]) for i in range(3)]
[nx, ny, nz] = eval(header[12].strip().split()[1])

#Create grids for x, y, z, used for plotting (_avg denotes cell-centered)
xg = np.linspace(xmin, xmax, nx+2);  dx = xg[1]-xg[0]
yg = np.linspace(ymin, ymax, ny+2);  dy = yg[1]-yg[0]
zg = np.linspace(zmin, zmax, nz+2);  dz = zg[1]-zg[0]
x_avg = (xg[1::]+xg[0:-1])/2
y_avg = (yg[1::]+yg[0:-1])/2
z_avg = (zg[1::]+zg[0:-1])/2
[X2D,Y2D] = np.meshgrid(x_avg,y_avg)
[X3D,Y3D,Z3D] = np.meshgrid(x_avg,y_avg,z_avg)

ts = yt.load(plotfilepath + '/plt*')    #Load all plotfiles
n = len(ts)    #Number of plotfiles (correspond to time slices)
ds = ts[0]
data = ds.all_data()
wt = data['particle_weight'].to_ndarray()
npart = np.sum((wt > 0) ==True)

x = np.zeros([npart,n]);     px = np.zeros([npart,n]);    Ex = np.zeros([npart,n]);
y = np.zeros([npart,n]);     py = np.zeros([npart,n]);    Ey = np.zeros([npart,n]);
z = np.zeros([npart,n]);     pz = np.zeros([npart,n]);    Ez = np.zeros([npart,n]);
wt = np.zeros([npart,n])

#Create arrays for fields, charge density, and along a diagonal slice
Exf = np.zeros([nx+1, ny+1, nz+1, n])
Eyf = np.zeros([nx+1, ny+1, nz+1, n])
Ezf = np.zeros([nx+1, ny+1, nz+1, n])
rho = np.zeros([nx+1, ny+1, nz+1, n])

#Constants and initial conditions
eps_0 = 8.8541878128e-12
m_e = 9.10938356e-31
q_e = -1.60217662e-19
pi = np.pi
r_0 = 0.1    #Initial distance of "massive" particle (does not move)
dt = 1e-6    #Timestep
t = [i*dt for i in range(n)]
q_tot = -1e-15

#Define functions for exact forms of v(r), t(r), Er(r) with r as the distance
#from the origin. Note, the "massive" particle is offset by -r_0 from the
#origin with the weightless particle opposite the origin with distance r. Thus,
#the initial distance between the particles is 2*r_0 (at r = r_0).
v_exact = lambda r: np.sqrt((q_e*q_tot)/(2*pi*m_e*eps_0)*(1/r_0-1/r))
t_exact = lambda r: np.sqrt(r_0**3*2*pi*m_e*eps_0/(q_e*q_tot)) \
    * (np.sqrt(r/r_0-1)*np.sqrt(r/r_0) \
       + np.log(np.sqrt(r/r_0-1)+np.sqrt(r/r_0)))
func = lambda tau: t_exact(tau)-t[-1]
r_end = fsolve(func,r_0)[0]
E_exact = lambda r: abs(q_tot)/(4*pi*eps_0*r**2)*(abs(r)>=r_end) \
    + abs(q_tot*r)/(4*pi*eps_0*r_end**3)*(abs(r)<r_end)

for i in range(len(ts)):    #Loop over plotfiles
    ds = ts[i]
    data = ds.all_data()

    #Read particle data for each plotfile
    x_temp = data['particle_position_x'].to_ndarray()
    y_temp = data['particle_position_y'].to_ndarray()
    z_temp = data['particle_position_z'].to_ndarray()
    px_temp = data['particle_momentum_x'].to_ndarray()
    py_temp = data['particle_momentum_y'].to_ndarray()
    pz_temp = data['particle_momentum_z'].to_ndarray()
    Ex_temp = data['particle_Ex'].to_ndarray()
    Ey_temp= data['particle_Ey'].to_ndarray()
    Ez_temp = data['particle_Ez'].to_ndarray()
    wt = data['particle_weight'].to_ndarray()

    rp = wt > 1e-3
    x[:,i] = x_temp[rp];    y[:,i] = y_temp[rp];    z[:,i] = z_temp[rp]
    px[:,i] = px_temp[rp];    py[:,i] = py_temp[rp];    pz[:,i] = pz_temp[rp]
    Ex[:,i] = Ex_temp[rp];    Ey[:,i] = Ey_temp[rp];    Ez[:,i] = Ez_temp[rp]

    #Read field data for each plotfile
    fdata = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                             dims=ds.domain_dimensions)
    Exf[:,:,:,i] = fdata['Ex'].to_ndarray()
    Eyf[:,:,:,i] = fdata['Ey'].to_ndarray()
    Ezf[:,:,:,i] = fdata['Ez'].to_ndarray()
    rho[:,:,:,i] = fdata['rho'].to_ndarray()

#Compute particle distance from origin
r0 = np.max(np.sqrt(x**2 + y**2 + z**2), axis=0)
r0_avg = (r0[1::]+r0[0:-1])/2    #Approx. position at half time steps
r0_avg = np.hstack((r0[0],r0_avg))

#Compute total particle momentum from origin (classical limit: v << c)
v0 = np.max(np.sqrt(px**2 + py**2 + pz**2), axis=0)/m_e

#Compute electric field magnitude along y=0 and z=0
Etot = np.sqrt(Exf**2 + Eyf**2 + Ezf**2)
E0x = (Etot[:,31,31,:] + Etot[:,31,32,:] + Etot[:,32,31,:] + Etot[:,32,32,:])/4

#Plotting routine for time, velocity, and electric field at particle location
fig1 = plt.figure(1)
plt.plot(r0,t_exact(r0))
plt.plot(r0,t,'o')
plt.xlabel('r [m]')
plt.ylabel('t [s]')
plt.title('Time vs sphere radius')

fig2 = plt.figure(2)
plt.plot(r0_avg,v_exact(r0_avg))
plt.plot(r0_avg,v0,'o')
plt.xlabel('r [m]')
plt.ylabel('v [m/s]')
plt.title('Velocity vs sphere radius')

fig3 = plt.figure(3)
plt.plot(x_avg,E_exact(x_avg))
plt.plot(x_avg,E0x[:,-1],'o')
plt.xlabel('r [m]')
plt.ylabel('E [V/m]')
plt.title('Electric field vs distance from origin at t=0')

#Create volume plot of field magnitude at first timestep
lvl = 2e-4
verts, faces, __, __ = measure.marching_cubes_lewiner\
    (Etot[:,:,:,30], level=lvl, spacing=(dx, dy, dz))
verts[:,0] = verts[:,0] + x_avg[0]
verts[:,1] = verts[:,1] + y_avg[0]
verts[:,2] = verts[:,2] + z_avg[0]

fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111, projection='3d')
ax4.plot_trisurf(verts[:,0], verts[:,1], faces, verts[:,2],
                cmap='Spectral', lw=1)
ax4.set_xlim((x_avg[0],x_avg[-1]))
ax4.set_ylim((y_avg[0],y_avg[-1]))
ax4.set_zlim((z_avg[0],z_avg[-1]))
ax4.set_xlabel('x [m]')
ax4.set_ylabel('y [m]')
ax4.set_zlabel('z [m]')
ax4.set_title('E-field magnitude isosurface at |E|=' + str(lvl))