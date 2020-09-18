#!/usr/bin/env python3

###
### Python script to plot particles transitioning between
### 3 different plasma stages - used for milestone report
###
### License: BSD-3-Clause-LBNL
###
### You can run it with the command:
### python YtPlotfileStages.py <Dimensions> <Path scan> <Path diags>
###


# Import statements
import os, glob, sys
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.constants as scc
import yt ; yt.funcs.mylog.setLevel(50)


# Command line input parameters
print('*', sys.argv[0], 'started running *\n')
if len(sys.argv) < 3:
    print('ABORT:')
    print(' Missing command line input parameters:\n'
          ' <number of dimensions of simulation>\n'
          ' <path to folder containing scan>\n'
          ' <path to folder containing beam lab frame data>\n')
    exit
else:
    if sys.argv[1] == '3':
        if3d = 1
    else:
        if3d = 0
    if sys.argv[2][-1] != '/':
        path = sys.argv[2] + '/'
    else:
        path = sys.argv[2]
    print('Going into ', path, 'directory')
    f_names = sys.argv[3:]
    lnames=len(f_names)

# Species list to plot
species_name = ['beam','electrons','electrons2','electrons3']
spc_colors=['r','g','y','k']
# Conversion factors
zfac=100 # To have z axis in cm
xfac=1e6 # in um
yfac=1e6
efac=1e-12 # To have fields in TV/m


species_colors = {}
for spci in range(len(species_name)):
    species_colors[species_name[spci]]=spc_colors[spci]

plot = {}
lplots={}
for f in f_names:
# Double check this address and the number of ? needed
    plot[f]=glob.glob(path+f+'/diags/diag?????')
#    plot[f]=glob.glob(path+f+'/diags/plotfiles/plt?????')
    plot[f].sort()
    lplots[f]=len(plot[f])
    print(f,lplots[f])
    if lplots[f]>1:
        for j in range(lplots[f]):
            plot[f][j]=plot[f][j][-9:]
#            plot[f][j]=plot[f][j][-8:]
        plot[f].remove(plot[f][0])
        plot[f].remove(plot[f][0])
        for i in range(5):
            plot[f].remove(plot[f][0])
            plot[f].remove(plot[f][-1])

    lplots[f]=len(plot[f])
print(plot)

# Pots settings
STAND_SIZE = 30
plt.rc('font', size=STAND_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=STAND_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=STAND_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=STAND_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=STAND_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=STAND_SIZE)    # legend fontsize
plt.rc('figure', titlesize=STAND_SIZE)   # fontsize of the figure title

# plot_path='/diags/plotfiles/'
plot_path='/diags/'

grids={}
g_num={}
extent={}
f_dirs={}
species={}
xp, yp, zp, indp = [{} for i in range(4)]
uxp, uyp, uzp, wp = [{} for i in range(4)]
for f in f_names:
    grids[f]={}
    g_num[f]={}
    extent[f]={}
    f_dirs[f]={}
    species[f]={}
    xp[f],yp[f],zp[f] = [{} for i in range(3)]
    uxp[f], uyp[f], uzp[f], wp[f] = [{} for i in range(4)]
    for p in plot[f]:
        f_dirs[f][p]=path+f+plot_path+p
        grids[f][p]={}
        g_num[f][p]={}
        species[f][p]=[]
        xp[f][p],yp[f][p],zp[f][p] = [{} for i in range(3)]
        uxp[f][p], uyp[f][p], uzp[f][p], wp[f][p] = [{} for i in range(4)]

def read_data(idir,ploti):
    f=f_names[idir]
    p=plot[f][ploti]

    print('Now reading ',f_dirs[f][p])
    ds = yt.load( f_dirs[f][p] )
#     print(ds.field_list)

    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge,
                                            dims=ds.domain_dimensions)
    ad = ds.all_data()

    if if3d == 1:
        dx_o2=25.0
        print("Slicing in y around ",dx_o2)

    if if3d == 1:
        extent[f][p]=([ds.domain_left_edge[2]*zfac,
                       ds.domain_right_edge[2]*zfac,
                       ds.domain_left_edge[0]*yfac,
                       ds.domain_right_edge[0]*yfac])
    else:
        extent[f][p]=([ds.domain_left_edge[1]*zfac,
               ds.domain_right_edge[1]*zfac,
               ds.domain_left_edge[0]*yfac,
               ds.domain_right_edge[0]*yfac])

    for sn in species_name:
        if sn in [i[0] for i in ds.field_list]:
            species[f][p].append(sn)
            if if3d == 1:
                yp0=ad[sn, 'particle_position_y'].v*yfac
                xp0=ad[sn, 'particle_position_x'].v*xfac
                if sn != 'beam':
                    select = xp0**2<(dx_o2)**2
                    yp[f][p][sn]=yp0[select]
                    del xp0, yp0
                    zp0=ad[sn, 'particle_position_z'].v*zfac
                    zp[f][p][sn]=zp0[select]
                    del zp0
                else:
                    wp[f][p][sn] = ad[sn, 'particle_weight'].v
                    xp[f][p][sn]=xp0
                    yp[f][p][sn]=yp0
                    del xp0, yp0
                    zp[f][p][sn]=ad[sn, 'particle_position_z'].v*zfac
            else:
                wp[f][p][sn] = ad[sn, 'particle_weight'].v
                yp[f][p][sn] = ad[sn, 'particle_position_x'].v*xfac
                zp[f][p][sn] = ad[sn, 'particle_position_y'].v*zfac
            uxp[f][p][sn] = ad['beam', 'particle_momentum_x'].v/scc.m_e/scc.c
            uyp[f][p][sn] = ad['beam', 'particle_momentum_y'].v/scc.m_e/scc.c
            uzp[f][p][sn] = ad['beam', 'particle_momentum_z'].v/scc.m_e/scc.c

    hasgrid=0
    del ds, ad, all_data_level_0
    return hasgrid, species

for idir in range(lnames):
    print(f_names[idir])
    for ploti in range(lplots[f_names[idir]]):
        hasgrid, species=read_data(idir,ploti)

def new_dir(npd):
    # Creating output folder
    if os.path.isdir(npd) == False:
        print('Created directory: ', npd)
        os.mkdir(npd)
    return npd

def z_lens(f, p):
    zl = [0.34*100, 0.69*100]
    gb = 60
#    bb = math.sqrt(1.-1./gb**2)
    time = extent[f][p][1].v
    zl_boosted = [(zl[i])/gb - time for i in range(2)]
    #[(zl[i]/gb-corner)/(scc.c*(1.+bb))*scc.c for i in range(2)]
    return zl_boosted

def plot_content(fig,ax,f,p,ifcb=1):

    for isn in range(len(species_name)):
        sn=species_name[isn]
        if sn in species[f][p] and sn != 'beam':
            ax.scatter(zp[f][p][sn][::10],yp[f][p][sn][::10],c=species_colors[sn],
                       s=1.0,marker='.',rasterized=True)

    snb='beam'
    if snb in species[f][p]:
        ax.scatter(zp[f][p][snb],yp[f][p][snb],c=species_colors[snb],
                   s=1.0,marker='.',rasterized=True)

    plt.xlim(extent[f][p][0],extent[f][p][1])
    plt.ylim(extent[f][p][2],extent[f][p][3])

    plt.xlabel('Propagation direction (z) [cm]')

    ax.set_ylim(-100,100)


def plot_data_2(f,plist,fdir):
    if len(plist)!=2:
        print("This function is to plot two plotfiles")
        return

    fig=plt.figure(figsize=(20, 7))
    ax1=fig.add_subplot(1,2,1)
    plot_content(fig,ax1,f,plot[f][plist[0]])
    ax1.set_ylabel('Transverse direction (x) [um]')
    plt.text(0.6, 80.0, 'a)',fontsize=STAND_SIZE, bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

    ax2=fig.add_subplot(1,2,2)
    plot_content(fig,ax2,f,plot[f][plist[1]])
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.text(0.9, 80.0, 'b)',fontsize=STAND_SIZE)

    p_name=fdir+'/BeamPlasmaStages'
    print('Saving plot ',p_name)

    lw = 2.0
    zlens  = z_lens(f,plot[f][plist[0]])
    ax1.axvline(x=zlens[0], linestyle='--', color='b', linewidth=lw)
    ax1.axvline(x=zlens[1], linestyle='--', color='b', linewidth=lw)
    zlens  = z_lens(f,plot[f][plist[1]])
    ax2.axvline(x=zlens[0], linestyle='--', color='b', linewidth=lw)
    ax2.axvline(x=zlens[1], linestyle='--', color='b', linewidth=lw)

    plt.savefig(p_name+'.png', dpi=300, bbox_inches = 'tight',pad_inches = 0.1)
    plt.savefig(p_name+'.pdf', dpi=300, bbox_inches = 'tight',pad_inches = 0.1)

    plt.show()
    plt.close()


for idir in range(lnames):
    npd=path+f_names[idir]+'/Plots-NoFld'
    print('Created output folder '+new_dir(npd))
    plot_data_2(f,[2,5],npd)

exit
