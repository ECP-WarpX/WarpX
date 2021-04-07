import yt
import os
from scipy.constants import mu_0
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

fig = plt.figure(1, (1,15))

grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 3),
                axes_pad = 1.0,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")


path, dirs, files = next(os.walk("../diags/"))
N_dirs = len(dirs)

for i in range(0, N_dirs):
    filename = '../diags/diag'+str(int(i)).zfill(5)+'/'
    print(filename)
    ds = yt.load(filename)
    comp = 'Ez'
    p = yt.SlicePlot( ds, 2, comp, origin='native' )
    p.set_log(comp, False)
    p.set_zlim(field = comp, zmin = -5000, zmax = 5000)
    plot = p.plots[comp]
    plot.figure = fig
    plot.axes = grid[0].axes
    plot.cax = grid.cbar_axes[0]
    p._setup_plots()

    comp = 'Bx'
    p = yt.SlicePlot( ds, 2, comp, origin='native' )
    p.set_log(comp, False)
    p.set_zlim(field = comp, zmin = -8*mu_0, zmax = 8*mu_0)
    plot = p.plots[comp]
    plot.figure = fig
    plot.axes = grid[1].axes
    plot.cax = grid.cbar_axes[1]
    p._setup_plots()

    comp = 'By'
    p = yt.SlicePlot( ds, 2, comp, origin='native' )
    p.set_log(comp, False)
    p.set_zlim(field = comp, zmin = -8*mu_0, zmax = 8*mu_0)
    plot = p.plots[comp]
    plot.figure = fig
    plot.axes = grid[2].axes
    plot.cax = grid.cbar_axes[2]
    p._setup_plots()

    fig.set_size_inches(13,5)
    p.save(str(i))
