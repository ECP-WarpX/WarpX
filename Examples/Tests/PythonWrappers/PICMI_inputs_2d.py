import numpy as np
import matplotlib.pyplot as plt
import pywarpx.fields as pwxf
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from pywarpx import picmi

# Number of time steps
max_steps = 100

# Grid
nx = 128
nz = 128

# Domain
xmin = 0.e-6
zmin = 0.e-6
xmax = 50.e-6
zmax = 50.e-6

# Cell size
dx = (xmax - xmin) / nx
dz = (zmax - zmin) / nz

# Domain decomposition
max_grid_size_x = 64
max_grid_size_z = 64

# PML
nxpml = 10
nzpml = 10
field_boundary = ['open', 'open']

# Spectral order
nox = 8
noz = 8

# Guard cells
nxg = 8
nzg = 8

# Initialize grid
grid = picmi.Cartesian2DGrid(number_of_cells = [nx,nz],
                             lower_bound = [xmin,zmin],
                             upper_bound = [xmax,zmax],
                             lower_boundary_conditions = field_boundary,
                             upper_boundary_conditions = field_boundary,
                             guard_cells = [nxg,nzg],
                             moving_window_velocity = [0.,0.,0],
                             warpx_max_grid_size_x = max_grid_size_x,
                             warpx_max_grid_size_y = max_grid_size_z)

# Initialize field solver
solver = picmi.ElectromagneticSolver(grid=grid, cfl=0.95, method='PSATD',
                                     stencil_order = [nox,noz],
                                     divE_cleaning = 1,
                                     divB_cleaning = 1,
                                     pml_divE_cleaning = 1,
                                     pml_divB_cleaning = 1,
                                     warpx_psatd_update_with_rho = True)

# Initialize diagnostics
diag_field_list = ["E", "B"]
field_diag = picmi.FieldDiagnostic(name = 'diag1',
                                   grid = grid,
                                   period = 10,
                                   write_dir = '.',
                                   warpx_file_prefix = 'Python_wrappers_plt',
                                   data_list = diag_field_list)

# Initialize simulation
sim = picmi.Simulation(solver = solver,
                       max_steps = max_steps,
                       verbose = 1,
                       particle_shape = 'cubic',
                       warpx_current_deposition_algo = 'direct',
                       warpx_particle_pusher_algo = 'boris',
                       warpx_field_gathering_algo = 'energy-conserving',
                       warpx_use_filter = 1)

# Add diagnostics to simulation
sim.add_diagnostic(field_diag)

# Write input file to run with compiled version
sim.write_input_file(file_name = 'inputs_2d')

# Whether to include guard cells in data returned by Python wrappers
include_ghosts = 1

# Compute min and max of fields data
def compute_minmax(data):
    vmin = -max(abs(data.min()), abs(data.max()))
    vmax =  max(abs(data.min()), abs(data.max()))
    return (vmin, vmax)

# Plot fields data either in valid domain or in PML
def plot_data(it, data, pml, title, name):
    fig, ax = plt.subplots(nrows = 1, ncols = 1, gridspec_kw = dict(wspace = 0.5), figsize = [6,5])
    cax = make_axes_locatable(ax).append_axes('right', size='5%', pad='5%')
    lw  = 0.8
    ls = '-'
    if (pml):
        # Draw PMLs and ghost regions
        ax.axvline(x = 0       , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axvline(x = 0+nxg   , linewidth = lw, color = 'black', linestyle = ls)
        ax.axvline(x = -nxpml  , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axvline(x = nx      , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axvline(x = nx-nxg  , linewidth = lw, color = 'black', linestyle = ls)
        ax.axvline(x = nx+nxpml, linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axhline(y = 0       , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axhline(y = 0+nzg   , linewidth = lw, color = 'black', linestyle = ls)
        ax.axhline(y = -nzpml  , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axhline(y = nz      , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axhline(y = nz-nzg  , linewidth = lw, color = 'black', linestyle = ls)
        ax.axhline(y = nz+nzpml, linewidth = lw, color = 'grey' , linestyle = ls)
        # Set extent and sliced data
        extent = np.array([-nxg-nxpml, nx+nxpml+nxg, -nzg-nzpml, nz+nzpml+nzg])
        X = data[-nxg-nxpml:nx+nxpml+nxg,-nzg-nzpml:nz+nzpml+nzg].transpose()
    else:
        # Draw ghost regions
        ax.axvline(x = 0       , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axvline(x = nx      , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axhline(y = 0       , linewidth = lw, color = 'grey' , linestyle = ls)
        ax.axhline(y = nz      , linewidth = lw, color = 'grey' , linestyle = ls)
        # Set extent and sliced data
        extent = np.array([-nxg, nx+nxg, -nzg, nz+nzg])
        X = data[-nxg:nx+nxg,-nzg:nz+nzg].transpose()
    # Min and max for colorbar
    vmin, vmax = compute_minmax(X)
    # Display data as image
    im = ax.imshow(X = X, origin = 'lower', extent = extent, vmin = vmin, vmax = vmax, cmap = 'seismic')
    # Add colorbar to plot
    fig.colorbar(im, cax = cax)
    # Set label for x- and y-axis, set title
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title(title)
    # Set plot title
    suptitle = 'PML in (x,z), 4 grids 64 x 64'
    plt.suptitle(suptitle)
    # Save figure
    figname = 'figure_' + name + '.png'
    fig.savefig(figname, dpi = 100)

# Initialize fields data (unit pulse) and apply smoothing
def init_data(data):
    data[nx//2,nz//2] = 1.
    data[2:nx-1,:] *= 0.5
    data[2:nx-1,:] += 0.5 * (data[1:nx-2,:] + data[3:nx,:])
    data[:,2:nz-1] *= 0.5
    data[:,2:nz-1] += 0.5 * (data[:,1:nz-2] + data[:,3:nz])

# Advance simulation for one time step
sim.step(1)

# Get fields data using Python wrappers
Ex = pwxf.ExFPWrapper(include_ghosts = include_ghosts)
Ey = pwxf.EyFPWrapper(include_ghosts = include_ghosts)
Ez = pwxf.EzFPWrapper(include_ghosts = include_ghosts)
Bx = pwxf.BxFPWrapper(include_ghosts = include_ghosts)
By = pwxf.ByFPWrapper(include_ghosts = include_ghosts)
Bz = pwxf.BzFPWrapper(include_ghosts = include_ghosts)
F  = pwxf.FFPWrapper(include_ghosts = include_ghosts)
G  = pwxf.GFPWrapper(include_ghosts = include_ghosts)
Expml = pwxf.ExFPPMLWrapper(include_ghosts = include_ghosts)
Eypml = pwxf.EyFPPMLWrapper(include_ghosts = include_ghosts)
Ezpml = pwxf.EzFPPMLWrapper(include_ghosts = include_ghosts)
Bxpml = pwxf.BxFPPMLWrapper(include_ghosts = include_ghosts)
Bypml = pwxf.ByFPPMLWrapper(include_ghosts = include_ghosts)
Bzpml = pwxf.BzFPPMLWrapper(include_ghosts = include_ghosts)
Fpml  = pwxf.FFPPMLWrapper(include_ghosts = include_ghosts)
Gpml  = pwxf.GFPPMLWrapper(include_ghosts = include_ghosts)

# Initialize fields data in valid domain
init_data(Ex)
init_data(Ey)
init_data(Ez)
init_data(Bx)
init_data(By)
init_data(Bz)
init_data(F)
init_data(G)

# Advance simulation until last time step
sim.step(max_steps-1)

# Plot E
plot_data(max_steps-1, Ex, edge, pml = False, title = 'Ex', name = 'Ex')
plot_data(max_steps-1, Ey, edge, pml = False, title = 'Ey', name = 'Ey')
plot_data(max_steps-1, Ez, edge, pml = False, title = 'Ez', name = 'Ez')

# Plot B
plot_data(max_steps-1, Bx, edge, pml = False, title = 'Bx', name = 'Bx')
plot_data(max_steps-1, By, edge, pml = False, title = 'By', name = 'By')
plot_data(max_steps-1, Bz, edge, pml = False, title = 'Bz', name = 'Bz')

# Plot E in PML
plot_data(max_steps-1, Expml[:,:,0], edge, pml = True, title = 'Exy in PML', name = 'Exy')
plot_data(max_steps-1, Expml[:,:,1], edge, pml = True, title = 'Exz in PML', name = 'Exz')
plot_data(max_steps-1, Expml[:,:,2], edge, pml = True, title = 'Exx in PML', name = 'Exx')
plot_data(max_steps-1, Eypml[:,:,0], edge, pml = True, title = 'Eyz in PML', name = 'Eyz')
plot_data(max_steps-1, Eypml[:,:,1], edge, pml = True, title = 'Eyx in PML', name = 'Eyx')
plot_data(max_steps-1, Eypml[:,:,2], edge, pml = True, title = 'Eyy in PML', name = 'Eyy') # zero
plot_data(max_steps-1, Ezpml[:,:,0], edge, pml = True, title = 'Ezx in PML', name = 'Ezx')
plot_data(max_steps-1, Ezpml[:,:,1], edge, pml = True, title = 'Ezy in PML', name = 'Ezy') # zero
plot_data(max_steps-1, Ezpml[:,:,2], edge, pml = True, title = 'Ezz in PML', name = 'Ezz')

# Plot B in PML
plot_data(max_steps-1, Bxpml[:,:,0], edge, pml = True, title = 'Bxy in PML', name = 'Bxy')
plot_data(max_steps-1, Bxpml[:,:,1], edge, pml = True, title = 'Bxz in PML', name = 'Bxz')
plot_data(max_steps-1, Bxpml[:,:,2], edge, pml = True, title = 'Bxx in PML', name = 'Bxx')
plot_data(max_steps-1, Bypml[:,:,0], edge, pml = True, title = 'Byz in PML', name = 'Byz')
plot_data(max_steps-1, Bypml[:,:,1], edge, pml = True, title = 'Byx in PML', name = 'Byx')
plot_data(max_steps-1, Bypml[:,:,2], edge, pml = True, title = 'Byy in PML', name = 'Byy') # zero
plot_data(max_steps-1, Bzpml[:,:,0], edge, pml = True, title = 'Bzx in PML', name = 'Bzx')
plot_data(max_steps-1, Bzpml[:,:,1], edge, pml = True, title = 'Bzy in PML', name = 'Bzy') # zero
plot_data(max_steps-1, Bzpml[:,:,2], edge, pml = True, title = 'Bzz in PML', name = 'Bzz')
