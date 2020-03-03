import yt
import matplotlib.pyplot as plt
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

yt.funcs.mylog.setLevel(50)    #Mute yt's verbose shenanigans
plotfilepath = '../Examples/Tests/TwoParticle_electrostatic/diags/plotfiles'

#Read header file for data to extract grid and cell data
with open(plotfilepath + '/plt00000/Header','r') as file:
    header = file.readlines()
[xmin, ymin, zmin] = [eval(header[9].strip().split()[i]) for i in range(3)]
[xmax, ymax, zmax] = [eval(header[10].strip().split()[i]) for i in range(3)]
[nx, ny, nz] = eval(header[12].strip().split()[1])

#Create grids for x, y, z, used for plotting (_avg denotes cell-centered)
x = np.linspace(xmin, xmax, nx+2);  dx = x[1]-x[0]
y = np.linspace(ymin, ymax, ny+2);  dy = y[1]-y[0]
z = np.linspace(zmin, zmax, nz+2);  dz = z[1]-z[0]
x_avg = (x[1::]+x[0:-1])/2
y_avg = (y[1::]+y[0:-1])/2
z_avg = (z[1::]+z[0:-1])/2
[X2D,Y2D] = np.meshgrid(x_avg,y_avg)
[X3D,Y3D,Z3D] = np.meshgrid(x_avg,y_avg,z_avg)

ts = yt.load(plotfilepath + '/plt*')    #Load all plotfiles
n = len(ts)    #Number of plotfiles (correspond to time slices)

#Create arrays for particle positions, momenta, and fields. Note: particle 1 is
#the "weightless" particle with no self fields and particle 2 is the "massive"
#particle which does not move.
x0 = [0]*n;     px0 = [0]*n;    Ex0 = [0]*n;
y0 = [0]*n;     py0 = [0]*n;    Ey0 = [0]*n;
z0 = [0]*n;     pz0 = [0]*n;    Ez0 = [0]*n;
x1 = [0]*n;     px1 = [0]*n;    Ex1 = [0]*n;
y1 = [0]*n;     py1 = [0]*n;    Ey1 = [0]*n;
z1 = [0]*n;     pz1 = [0]*n;    Ez1 = [0]*n;

#Create arrays for fields, charge density, and along a diagonal slice
Ex = np.zeros([nx+1, ny+1, nz+1, n])
Ey = np.zeros([nx+1, ny+1, nz+1, n])
Ez = np.zeros([nx+1, ny+1, nz+1, n])
rho = np.zeros([nx+1, ny+1, nz+1, n])

#Constants and initial conditions
eps_0 = 8.8541878128e-12
m_e = 9.10938356e-31
q_e = -1.60217662e-19
pi = np.pi
r_0 = 0.1    #Initial distance of "massive" particle (does not move)
dt = 2e-4    #Timestep
t = [i*dt for i in range(n)]

#Define functions for exact forms of v(r), t(r), Er(r) with r as the distance
#from the origin. Note, the "massive" particle is offset by -r_0 from the
#origin with the weightless particle opposite the origin with distance r. Thus,
#the initial distance between the particles is 2*r_0 (at r = r_0).
v_exact = lambda r: np.sqrt((q_e**2/(2*pi*m_e*eps_0))*(1/(2*r_0)-1/(r+r_0)))
t_exact = lambda r: np.sqrt(r_0**3*16*pi*m_e*eps_0/(q_e**2)) \
    * (np.sqrt((r+r_0)/(2*r_0)-1)*np.sqrt((r+r_0)/(2*r_0)) \
       + np.log(np.sqrt((r+r_0)/(2*r_0)-1) + np.sqrt((r+r_0)/(2*r_0))))
E_exact = lambda r: abs(q_e)/(4*pi*eps_0*abs(r+r_0)**2)

for i in range(len(ts)):    #Loop over plotfiles
    ds = ts[i]
    data = ds.all_data()
    
    #Read particle data for each plotfile
    [x0[i], x1[i]] = data['particle_position_x'].to_ndarray()
    [y0[i], y1[i]] = data['particle_position_y'].to_ndarray()
    [z0[i], z1[i]] = data['particle_position_z'].to_ndarray()
    [px0[i], px1[i]] = data['particle_momentum_x'].to_ndarray()
    [py0[i], py1[i]] = data['particle_momentum_y'].to_ndarray()
    [pz0[i], pz1[i]] = data['particle_momentum_z'].to_ndarray()
    [Ex0[i], Ex1[i]] = data['particle_Ex'].to_ndarray()
    [Ey0[i], Ey1[i]] = data['particle_Ey'].to_ndarray()
    [Ez0[i], Ez1[i]] = data['particle_Ez'].to_ndarray()
    
    #Read field data for each plotfile
    fdata = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, 
                             dims=ds.domain_dimensions)
    Ex[:,:,:,i] = fdata['Ex'].to_ndarray()
    Ey[:,:,:,i] = fdata['Ey'].to_ndarray()
    Ez[:,:,:,i] = fdata['Ez'].to_ndarray()
    rho[:,:,:,i] = fdata['rho'].to_ndarray()

#Compute particle distance from origin
r0 = np.sqrt(np.power(x0,2) + np.power(y0,2) + np.power(z0,2))
r0_avg = (r0[1::]+r0[0:-1])/2    #Approx. position at half time steps
r0_avg = np.hstack((r0[0],r0_avg))

#Compute total particle momentum from origin (classical limit: v << c)
v0 = np.sqrt(np.power(px0,2) + np.power(py0,2) + np.power(pz0,2))/m_e

#Compute radial electric field magnitude
E0 = np.sqrt(np.power(Ex0,2) + np.power(Ey0,2) + np.power(Ez0,2))

#Plotting routine for time, velocity, and electric field at particle location
fig1 = plt.figure(1)
plt.plot(r0,t_exact(r0))
plt.plot(r0,t,'o')
plt.xlabel('r [m]')
plt.ylabel('t [s]')
plt.title('Time vs displacement (from origin) of weightless particle')

fig2 = plt.figure(2)
plt.plot(r0_avg,v_exact(r0_avg))
plt.plot(r0_avg,v0,'o')
plt.xlabel('r [m]')
plt.ylabel('v [m/s]')
plt.title('Velocity vs displacement (from origin) of weightless particle')

fig3 = plt.figure(3)
plt.plot(r0,E_exact(r0))
plt.plot(r0,E0,'o')
plt.xlabel('r [m]')
plt.ylabel('E [V/m]')
plt.title('Electric field vs displacement (from origin) of weightless particle')

#Create volume plot of field magnitude at first timestep
Etot = np.sqrt(np.power(Ex[:,:,:,1],2) + np.power(Ey[:,:,:,1],2) 
             + np.power(Ez[:,:,:,1],2))
verts, faces, __, __ = measure.marching_cubes_lewiner\
    (Etot[::2,::2,::2],level=E0[-1], spacing=(dx*2, dy*2, dz*2))
verts[:,0] = verts[:,0]+x_avg[0]
verts[:,1] = verts[:,1]+y_avg[0]
verts[:,2] = verts[:,2]+z_avg[0]
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111, projection='3d')
ax4.plot_trisurf(verts[:,0], verts[:,1], faces, verts[:,2],
                cmap='Spectral', lw=1)
ax4.set_xlim((x_avg[0],x_avg[-1]))
ax4.set_ylim((y_avg[0],y_avg[-1]))
ax4.set_zlim((z_avg[0],z_avg[-1]))
ax4.set_aspect('equal')

#Create 3d quiver plot of field at first timestep
fig5 = plt.figure(5)
ax5 = fig5.add_subplot(111, projection='3d')
ax5.quiver(X3D[::8,::8,::8],Y3D[::8,::8,::8],Z3D[::8,::8,::8],
            Ex[::8,::8,::8,1]*3e6,Ey[::8,::8,::8,1]*3e6,Ez[::8,::8,::8,1]*3e6)



