#!/usr/bin/env python3

###
### Python script to get the beam properties from
### .h5 files of a scan and make comparison plots
###
### L.D.Amorim @ LBNL
### 3 Sep 2020
###
### You can run it with the command:
### python PropPlotsFromHDF5.py <Dimensions> <Path scan> <Paths lab data>
###


# Import statements
import os, sys
import numpy as np
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import h5py

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
        path_scan = sys.argv[2] + '/'
    else:
        path_scan = sys.argv[2]
    print('Going into ', path_scan, 'directory')
    sim_list = sys.argv[3:]
    nsims=len(sim_list)

# Input subfolder where .h5 files are stored
path_hdf5 = 'Beam-properties'
# Output folder where plots are stored
dir_output = 'ConvScan-Plots'
# Beam STD width in m
beam_w = 6e-7
# Laser wavelength in m
laser_w = 8.e-7
# Boosted frame Lorentz factor
gb = 60
bb = math.sqrt(1.-1./gb**2)

# Getting resolution information from runs names
dx_list = []
dz_list = []
ppw = {}
ppl = {}
for sim in sim_list:
    if sim[3] == 'p':
        decimal = len(sim.split('um')[0])
        dx_list.append((float(sim[2])+float(sim[4:decimal])/10**(decimal-4))*1e-6) # in m
    else:
        dx_list.append((float(sim[2]))*1e-6) # in m
    dz_list.append(dx_list[-1]/(1+bb)/gb)
    ppw[sim] = beam_w/dx_list[-1]
    ppl[sim] = laser_w/dz_list[-1]
print(dx_list)
print(dz_list)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
#print(colors)

fstyles=['-' for i in range(nsims)]
fmarkers=['P' for i in range(nsims)]
for i in range(nsims):
    if 'FDTD' in sim_list[i]:
        fstyles[i]='--'
        fmarkers[i]='o'
        for j in range(nsims):
            if sim_list[i] != sim_list[j]:
                if sim_list[i][:5] == sim_list[j][:5]:
                    colors[i]=colors[j]
f_colors = {}
f_styles = {}
f_markers = {}
lbl = {}
for isim in range(nsims):
    f_colors[sim_list[isim]] = colors[isim]
    f_styles[sim_list[isim]] = fstyles[isim]
    f_markers[sim_list[isim]] = fmarkers[isim]
    if isim == 0:
        lbl[sim_list[isim]] = 'PSATD'
    elif isim == 1:
        lbl[sim_list[isim]] = 'FDTD'
    else:
        lbl[sim_list[isim]] = None
print(f_colors)

# Creating output directory for the plots
new_dir=path_scan + dir_output
if os.path.isdir(new_dir) == False:
            print('Created directory: ', new_dir)
            os.mkdir(new_dir)

# Function to load data from .h5 files
def read_sim_hdf5(sim):
    path_file = (path_scan + sim + '/lab_frame_data/' + path_hdf5
                  + '/' + path_hdf5 + '-' + sim + '.hdf5')
    print('Now downloading data from: ', path_file)
    fin=h5py.File(path_file,'r')
    datasets=list(fin.keys())
    print(datasets)
    d_data={}
    for ds in datasets:
        d_data[ds]={}
        d_data[ds]['units'] = fin[ds].attrs['units']
        d_data[ds]['label'] = fin[ds].attrs['label']
        d_data[ds]['range'] = fin[ds].attrs['range']
        d_data[ds]['data']=fin[ds].value
    fin.close()

    return d_data

d_z      = {}
d_totq   = {}
d_totp   = {}
d_emean  = {}
d_estd   = {}
d_emitpx = {}
d_emitx  = {}
d_emitpy = {}
d_emity  = {}
d_bwdx   = {}
d_bwdy   = {}
d_bwdux   = {}
d_bwduy   = {}

for sim in sim_list:
    d_data = read_sim_hdf5(sim)
    d_z     [sim] = d_data['z']['data']
    d_totq  [sim] = d_data['totq']['data']
    d_totp  [sim] = d_data['totp']['data']
    d_emean [sim] = d_data['emean']['data']
    d_estd  [sim] = d_data['estd']['data']
    d_emitpx[sim] = d_data['emitpx']['data']
    d_emitx [sim] = d_data['emitx']['data']
    d_bwdx  [sim] = d_data['xsig']['data']
    d_bwdux [sim] = d_data['uxsig']['data']
    if (if3d == 1):
        d_emitpy[sim] = d_data['emitpy']['data']
        d_emity [sim] = d_data['emity']['data']
        d_bwdy  [sim] = d_data['ysig']['data']
        d_bwduy [sim] = d_data['uysig']['data']

# Plasma profile
upramp = .02
plateau = .295
gap = 0.03
downramp = 0.005
start = [0.0, gap+upramp+plateau+downramp, (gap+upramp+plateau+downramp)*2]
taxis = [start[0],start[0]+upramp, start[0]+upramp+plateau, start[0]+upramp+plateau+downramp,
         start[1],start[1]+upramp, start[1]+upramp+plateau, start[1]+upramp+plateau+downramp,
         start[2],start[2]+upramp, start[2]+upramp+plateau, start[2]+upramp+plateau+downramp]
taxis = [x for x in taxis]
zeros = [0 for x in taxis]
f_plasma1=np.asarray([0,1.0,1.0,0,0,0,0,0,0,0,0,0])
f_plasma2=np.asarray([0,0,0,0,0,1.0,1.0,0,0,0,0,0])
f_plasma3=np.asarray([0,0,0,0,0,0,0,0,0,1.0,1.0,0])
zlens= [taxis[3]+0.02, taxis[7]+0.02]
print(zlens)
print(taxis)

# Plots settings
STAND_SIZE=30
plt.rc('font', size=STAND_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=STAND_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=STAND_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=STAND_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=STAND_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=STAND_SIZE)    # legend fontsize
plt.rc('figure', titlesize=STAND_SIZE)   # fontsize of the figure title

def p_subplot_inset(ind,x,y,xl,yl,ax):
    kl=[sim_list[ind[i]] for i in range(len(ind))]
    fp = 0.0
    for sim in kl:
        ax.plot(x[sim],y[sim],f_styles[sim],label=lbl[sim],c=f_colors[sim])
        fp = max(fp,max(y[sim]))
    fp=fp*0.2
    ax.fill_between(taxis, zeros, f_plasma1*fp, facecolor='brown', alpha=.1)
    ax.fill_between(taxis, zeros, f_plasma2*fp, facecolor='brown', alpha=.1)
    ax.fill_between(taxis, zeros, f_plasma3*fp, facecolor='brown', alpha=.1)
    plt.xlabel(xl)
    plt.ylabel(yl)
# Plotting lenses position
    plt.axvline(x=zlens[0], color='k',linewidth=0.2)
    plt.axvline(x=zlens[1], color='k',linewidth=0.2)
    return 0


def plot_en_res(ind,fname):
    plt.figure(figsize=(24,6))

    ax0=plt.subplot(1,3,1)
    p_subplot_inset(ind,d_z,d_emean,
                    'Propagation direction [m]','Mean energy [GeV]',ax0)
    plt.ylim(0,22)
    plt.legend(loc='lower right',bbox_to_anchor=(0.65,0.65))
    ax1=plt.subplot(1,3,2)
    p_subplot_inset(ind,d_z,d_estd,
                    'Propagation direction [m]','Energy spread [%]', ax1)
    plt.ylim(0,13.0)
    ax2 = plt.subplot(1,3,3)
    size = 300
    for sim in sim_list:
        if 'FDTD' in sim:
            size = 200
        ax2.scatter(ppl[sim],ppw[sim],marker=f_markers[sim],c=f_colors[sim],s=size)
    ax2.set_xlabel('Points per laser wavelength')
    ax2.set_ylabel('Points per beam width')

    plt.subplots_adjust(hspace=0.4,wspace=0.35)

    plt.savefig(new_dir+'/'+fname+'-all.pdf', bbox_inches='tight')
    plt.savefig(new_dir+'/'+fname+'-all.png', bbox_inches='tight')
    print('Saving file ',new_dir+'/'+fname+'-all.pdf')

    plt.close()
    return 0


def plot_all_props(ind,fname):
    plt.figure(figsize=(24,18))

    ax0=plt.subplot(3,3,1)
    p_subplot_inset(ind,d_z,d_emean,
                    'Propagation direction [m]','Mean energy [GeV]',ax0)
    ax1=plt.subplot(3,3,2)
    p_subplot_inset(ind,d_z,d_emitx, 'Propagation direction [m]',
                    '<slice $\epsilon$ in x> (um)', ax1)
    plt.legend(loc='lower right',bbox_to_anchor=(0.65,0.5))
    ax2=plt.subplot(3,3,3)
    p_subplot_inset(ind,d_z,d_emity, 'Propagation direction [m]',
                    '<slice $\epsilon$ in y> [um]', ax2)
    ax0.get_xaxis().set_visible(False)
    ax1.get_xaxis().set_visible(False)
    ax2.get_xaxis().set_visible(False)

    ax3=plt.subplot(3,3,4)
    p_subplot_inset(ind,d_z,d_estd,
                     'Propagation direction [m]','Energy spread [%]', ax3)
    ax4=plt.subplot(3,3,5)
    p_subplot_inset(ind,d_z,d_bwdx,'Propagation direction [m]',
                    'Beam width in x [um]',  ax4)
    ax5=plt.subplot(3,3,6)
    p_subplot_inset(ind,d_z,d_bwdy, 'Propagation direction [m]',
                    'Beam width in y [um]', ax5)
    ax6=plt.subplot(3,3,7)
    p_subplot_inset(ind,d_z,d_totq, 'Propagation direction [m]',
                    'Beam charge [nC]', ax6)

    ax7=plt.subplot(3,3,8)
    size = 300
    for sim in sim_list:
        if 'FDTD' in sim:
            size = 200
        ax7.scatter(ppl[sim],ppw[sim],marker=f_markers[sim],c=f_colors[sim],s=size)
    ax7.set_xlabel('Points per laser wavelength')
    ax7.set_ylabel('Points per beam width')

    plt.subplots_adjust(hspace=0.4,wspace=0.35)

    plt.savefig(new_dir+'/'+fname+'-all.pdf', bbox_inches='tight')
    plt.savefig(new_dir+'/'+fname+'-all.png', bbox_inches='tight')
    print('Saving file ',new_dir+'/'+fname+'-all.pdf')

    plt.close()
    return 0

# Plotting the beam properties of all the runs
indtot=list(range(nsims))
plot_en_res(indtot,'all_en_scan')
plot_all_props(indtot,'all_props_scan')

exit
