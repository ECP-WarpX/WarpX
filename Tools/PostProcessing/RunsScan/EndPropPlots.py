#!/usr/bin/env python3

###
### Python script to get the beam properties from
### .h5 files of a scan and make comparison plots
###
### L.D.Amorim @ LBNL
### 3 Sep 2020
###
### You can run it with the command:
### python EndPropPlots.py <Path read_raw_data> <Path scan> <Paths lab data>
###

# Import statements
import os, glob, sys
import numpy as np
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.constants as scc
import math
import h5py

# Command line input parameters
print('*', sys.argv[0], 'started running *\n')
if len(sys.argv) < 3:
    print('ABORT:')
    print(' Missing command line input parameters:\n'
          ' <path to folder containing scan>\n'
          ' <path to folder containing beam lab frame data>\n')
    exit
else:
    if sys.argv[1][-1] != '/':
        path_scan = sys.argv[1] + '/'
    else:
        path_scan = sys.argv[1]
    print('Going into ', path_scan, 'directory')
    sim_list = sys.argv[2:]
    nsims=len(sim_list)

# Input subfolder where .h5 files are stored
path_hdf5 = 'Beam-properties'
# Output folder where plots are stored
dir_output = 'ConvScan-Plots'
# If the run is in 2D geometry
if2d = 0
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


# Creating output directory for the plots
new_dir=path_scan + dir_output
if os.path.isdir(new_dir) == False:
            print('Created directory: ', new_dir)
            os.mkdir(new_dir)


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
    if (if2d == 0):
        d_emitpy[sim] = d_data['emitpy']['data']
        d_emity [sim] = d_data['emity']['data']
        d_bwdy  [sim] = d_data['ysig']['data']
        d_bwduy [sim] = d_data['uysig']['data']



# Plots settings
STAND_SIZE=30
plt.rc('font', size=STAND_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=STAND_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=STAND_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=STAND_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=STAND_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=STAND_SIZE)    # legend fontsize
plt.rc('figure', titlesize=STAND_SIZE)   # fontsize of the figure title

def p_par_evol(fname):
    plot=plt.figure(figsize=(24,6)) #12))

    ax0=plt.subplot(1,3,1)
    ax0.set_ylabel('Mean energy (GeV)')
    ax1=plt.subplot(1,3,2)
    ax1.set_ylabel('Energy spread (%)')
    ax2=plt.subplot(1,3,3)
    ax2.set_ylabel('Points per beam width')
#     ax3=plt.subplot(2,3,4)
#     ax3.set_ylabel('Emittance in x (um)')
#     ax4=plt.subplot(2,3,5)
#     ax4.set_ylabel('Emittance in y (um)')
#     ax5=plt.subplot(2,3,6)
#     ax5.set_ylabel('Charge [pC]') #(R < 3RMS)
    
    xins=ppl
    xinsl='Points per laser wavelength'
    
    arr_fdtd = np.asarray([i for i in range(nsims) if 'FDTD' in sim_list[i]])
    c_fdtd = '#1f77b4'
    s_fdtd = 'o'
    size_fdtd = 300
    arr_psatd = np.asarray([i for i in range(nsims) if 'FDTD' not in sim_list[i]])
    c_psatd = '#ff7f0e'
    s_psatd = 'P'
    size_psatd = 250
    
    x=np.asarray([xins[sim] for sim in sim_list])
    y0=np.asarray([d_emean[sim][-1] for sim in sim_list])
    y1=np.asarray([d_estd[sim][-1] for sim in sim_list])
    y2=np.asarray([ppw[sim] for sim in sim_list])
#     y3=np.asarray([d_emitx[sim][-1] for sim in sim_list])
#     y4=np.asarray([d_emity[sim][-1] for sim in sim_list])
#     y5=np.asarray([d_totq[sim][-1] for sim in sim_list])
    
    ax0.plot(x[arr_fdtd[:]],y0[arr_fdtd[:]],color=c_fdtd)
    ax0.plot(x[arr_psatd[:]],y0[arr_psatd[:]],color=c_psatd)
    ax1.plot(x[arr_fdtd[:]],y1[arr_fdtd[:]],color=c_fdtd)
    ax1.plot(x[arr_psatd[:]],y1[arr_psatd[:]],color=c_psatd)
    ax2.plot(x[arr_fdtd[:]],y2[arr_fdtd[:]],color=c_fdtd)
    ax2.plot(x[arr_psatd[:]],y2[arr_psatd[:]],color=c_psatd)
#     ax3.plot(x[arr_fdtd[:]],y3[arr_fdtd[:]],color=c_fdtd)
#     ax3.plot(x[arr_psatd[:]],y3[arr_psatd[:]],color=c_psatd)
#     ax4.plot(x[arr_fdtd[:]],y4[arr_fdtd[:]],color=c_fdtd)
#     ax4.plot(x[arr_psatd[:]],y4[arr_psatd[:]],color=c_psatd)
#     ax5.plot(x[arr_fdtd[:]],y5[arr_fdtd[:]],color=c_fdtd)
#     ax5.plot(x[arr_psatd[:]],y5[arr_psatd[:]],color=c_psatd)
    
    ax0.scatter(x[arr_fdtd[:]], y0[arr_fdtd[:]], c=c_fdtd,s=size_fdtd,marker=s_fdtd)
    ax0.scatter(x[arr_psatd[:]], y0[arr_psatd[:]], c=c_psatd,s=size_psatd,marker=s_psatd)
    ax1.scatter(x[arr_fdtd[:]], y1[arr_fdtd[:]], c=c_fdtd,s=size_fdtd,marker=s_fdtd)
    ax1.scatter(x[arr_psatd[:]], y1[arr_psatd[:]], c=c_psatd,s=size_psatd,marker=s_psatd)
    ax2.scatter(x[arr_fdtd[:]], y2[arr_fdtd[:]], c=c_fdtd,s=size_fdtd,marker=s_fdtd, label = 'FDTD')
    ax2.scatter(x[arr_psatd[:]], y2[arr_psatd[:]], c=c_psatd,s=size_psatd,marker=s_psatd, label = 'PSATD')
#     ax3.scatter(x[arr_fdtd[:]], y3[arr_fdtd[:]], c=c_fdtd,s=size_fdtd,marker=s_fdtd)
#     ax3.scatter(x[arr_psatd[:]], y3[arr_psatd[:]], c=c_psatd,s=size_psatd,marker=s_psatd)
#     ax4.scatter(x[arr_fdtd[:]], y4[arr_fdtd[:]], c=c_fdtd,s=size_fdtd,marker=s_fdtd)
#     ax4.scatter(x[arr_psatd[:]], y4[arr_psatd[:]], c=c_psatd,s=size_psatd,marker=s_psatd)
#     ax5.scatter(x[arr_fdtd[:]], y5[arr_fdtd[:]], c=c_fdtd,s=size_fdtd,marker=s_fdtd)
#     ax5.scatter(x[arr_psatd[:]], y5[arr_psatd[:]], c=c_psatd,s=size_psatd,marker=s_psatd)
    
    ax0.set_xlabel('Points per laser wavelength')
    ax1.set_xlabel('Points per laser wavelength')
    ax2.set_xlabel('Points per laser wavelength')

    leg=ax2.legend(loc='lower left',bbox_to_anchor=(0.0,0.55))
        
    plt.subplots_adjust(hspace=0.3,wspace=0.4)

    plt.savefig(new_dir+'/'+fname+'-all.pdf', bbox_inches='tight')
    plt.savefig(new_dir+'/'+fname+'-all.png', bbox_inches='tight')
    print('Saving file ',new_dir+'/'+fname+'-all.pdf')
    
    plt.close()
    
    return 0

p_par_evol('all_end_scan')

exit 
