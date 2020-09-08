#!/usr/bin/env python3

###
### Python script to evaluate lab frame particle data
### and save the beam (RMS) properties in an .h5 file
###
### L.D.Amorim @ LBNL
### 3 Sep 2020
###
### You can run it with the command:
### python BeamPropToHDF5.py <Dimensions> <Path scan> <Paths lab data>
###


# Import statements
import os, glob, sys
import numpy as np
import scipy.constants as scc
import h5py


# Command line input parameters
print('*', sys.argv[0], 'started running *\n')
if len(sys.argv) < 4:
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
    print('Using default location /lab_frame_data/snapshots/snapshot?????')
    sim_list = sys.argv[3:]
    nsims=len(sim_list)


# Output subfolder where .h5 file is stored
dir_hdf5 = 'Beam-properties'

# Name of WarpX output species
species = 'beam'
# Number of RMS radial distances from beam centroid
d_rms = 3


# Definition of the functions used to compute beam selected particles
# and beam properties along its propagation
def get_particle_field(snapshot, species, field):
    fn = snapshot + '/' + species
    files = glob.glob(os.path.join(fn, field + '_*'))
    files.sort()
    all_data = np.array([])
    for f in files:
        data = np.fromfile(f)
        all_data = np.concatenate((all_data, data))
    return all_data

def w_ave( a, weights ):
    if not np.any(weights) and not np.any(a):
        return np.nan
    else:
        average = np.average(a, weights=weights)
        return( average )

# Function to select the particles
def rem_part(w,x,z,ux,uy,uz,y=[0]):
    if not np.any(w):
        if (if3d == 1):
            return w,x,y,z,ux,uy,uz
        else:
            return w,x,z,ux,uy,uz
    else:
        npart=len(w)
        xm = w_ave( x, w )
        xsq = w_ave( (x - xm) ** 2, w )
        zm = w_ave( z, w )
        zsq = w_ave( (z - zm) ** 2, w )
        indx=np.array([i for i in range(npart) if ((x[i]-xm)**2<d_rms**2*xsq)])
        indz=np.array([i for i in range(npart) if ((z[i]-zm)**2<d_rms**2*zsq)])
        if (if3d == 1):
            ym = w_ave( y, w )
            ysq = w_ave( (y - ym) ** 2, w )
            indy=np.array([i for i in range(npart) if ((y[i]-ym)**2<d_rms**2*ysq)])
            indexes0=np.intersect1d(indx,indz,return_indices=True)[0]
            indexes=np.intersect1d(indexes0,indy,return_indices=True)[0]
        else:
            indexes=np.intersect1d(indx,indz,return_indices=True)[0]
        if len(indexes)>0:
            w=w[indexes[:]]
            x=x[indexes[:]]
            z=z[indexes[:]]
            ux=ux[indexes[:]]
            uy=uy[indexes[:]]
            uz=uz[indexes[:]]
            if (if3d == 1):
                y=y[indexes[:]]
                return w,x,y,z,ux,uy,uz
            else:
                return w,x,z,ux,uy,uz
        else:
            return np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1),np.zeros(1)

# General emittance calculations
def emittance_from_coord(x, ux, uy, w,y=[0]):
    xm = w_ave( x, w )
    xsq = w_ave( (x-xm) ** 2, w )
    if (if3d == 1):
        ym = w_ave( y, w )
        ysq = w_ave( (y-ym) ** 2, w )
    uxm = w_ave( ux, w )
    uxsq = w_ave( (ux-uxm) ** 2, w )
    uym = w_ave( uy, w )
    uysq = w_ave( (uy-uym) ** 2, w )
    xux = w_ave( (x-xm) * (ux-uxm), w )
    emit_x = ( abs(xsq * uxsq - xux ** 2) )**.5
    if (if3d == 1):
        yuy = w_ave( (y-ym) * (uy-uym), w )
        emit_y = ( abs(ysq * uysq - yuy ** 2) )**.5
        return emit_x, emit_y
    else:
        return emit_x

# Function to compute beam emittance
def my_emittance(x, z, px, py, pz, w=None, kind='norm', sliced=False,
                nslices=20, beam_length=None, y=[0]):
    if w is None:
        w = np.ones_like(x)
    if sliced == False:
        xm = np.average( x, weights=w )
        xsq = np.average( (x-xm) ** 2, weights=w )
        if (if3d == 1):
            ym = np.average( y, weights=w )
            ysq = np.average( (y-ym) ** 2, weights=w )
        if kind == 'norm':
            ux = px
            uy = py
        elif kind == 'trace':
            ux = px/pz
            uy = py/pz
        uxm    = np.average( ux, weights=w )
        uxsq   = np.average( (ux-uxm) ** 2, weights=w )
        uym    = np.average( uy, weights=w )
        uysq   = np.average( (uy-uym) ** 2, weights=w )
        xux    = np.average( (x-xm) * (ux-uxm), weights=w )
        gamma = np.average(np.sqrt( 1 + px**2 + py**2 + pz**2 ), weights=w)
        emit_x = np.sqrt( abs(xsq * uxsq - xux ** 2) )
        if (if3d == 1):
            yuy    = np.average( (y-ym) * (uy-uym), weights=w )
            emit_y = np.sqrt( abs(ysq * uysq - yuy ** 2) )
            return emit_x, emit_y, gamma
        else:
            return emit_x, gamma
    if sliced == True:
        zavg = np.average(z, weights=w)
        z = z - zavg
        if beam_length is None:
            zmin = np.amin(z)
            zmax = np.amax(z)
            bins = np.linspace(zmin, zmax, nslices)
        else:
            bins = np.linspace(-beam_length/2, beam_length/2, nslices)
        #print()
        EMIT_X = np.zeros(len(bins[:-1]))
        if (if3d == 1):
            EMIT_Y = np.zeros(len(bins[:-1]))
        GAMMA = np.zeros(len(bins[:-1]))
        W_slice = np.zeros(len(bins[:-1]))
        for count, leftedge in enumerate(bins[:-1]):
            binsmean =  .5*bins[count]+.5*bins[count+1]
            binswidth = .5*bins[count+1]-.5*bins[count]
            select = ( np.abs(z-binsmean)<=binswidth )
            W_slice[count] = np.sum(w[select])
            if W_slice[count]>0:
                if (if3d == 1):
                    EMIT_X[count], EMIT_Y[count],GAMMA[count] = my_emittance(x[select], z[select],
                                                                px[select], py[select], pz[select],
                                                                w[select], kind=kind, y=y[select])
                else:
                    EMIT_X[count], GAMMA[count] = my_emittance(x[select], z[select],
                                                                px[select], py[select], pz[select],
                                                                w[select], kind=kind)
        if (if3d == 1):
            return EMIT_X, EMIT_Y, W_slice, GAMMA
        else:
            return EMIT_X, W_slice, GAMMA

# Function used to read the beam data and call previous functions
def read_sim(sim):
    print('Now downloading data from: ', sim)
    path_sim = path_scan + sim + '/'
    file_list = glob.glob(path_sim + '/lab_frame_data/snapshots/snapshot?????')
    file_list.sort()
    nfiles = len(file_list)#-1
    iteration_list = range(nfiles)
    zz         = np.zeros(nfiles)
    TotP       = np.zeros(nfiles)
    TotQ       = np.zeros(nfiles)
    Emean      = np.zeros(nfiles)
    Estd       = np.zeros(nfiles)
    beam_widthx= np.zeros(nfiles)
    beam_widthy= np.zeros(nfiles)
    beam_widthux= np.zeros(nfiles)
    beam_widthuy= np.zeros(nfiles)
    beam_widthy_center= np.zeros(nfiles)
    beam_div   = np.zeros(nfiles)
    emittancex = np.zeros(nfiles)
    emittancepx= np.zeros(nfiles)
    emittancey = np.zeros(nfiles)
    emittancepy= np.zeros(nfiles)
    nslices = 20
    EMITX      =  np.zeros((nfiles, nslices-1))
    EMITY      =  np.zeros((nfiles, nslices-1))
    W_slice    =  np.zeros((nfiles, nslices-1))
    for count, iteration in enumerate(iteration_list):
        snapshot = path_sim + '/lab_frame_data/snapshots/' + 'snapshot' + str(iteration).zfill(5)
        w = get_particle_field(snapshot, species, 'w')
        x = get_particle_field(snapshot, species, 'x')
        if (if3d == 1):
            y = get_particle_field(snapshot, species, 'y')
        z = get_particle_field(snapshot, species, 'z')
        ux = get_particle_field(snapshot, species, 'ux')/scc.c
        uy = get_particle_field(snapshot, species, 'uy')/scc.c
        uz = get_particle_field(snapshot, species, 'uz')/scc.c

        TotP      [count] = w.shape[0]
        # Removing electrons that are no longer in the beam (>d_rms RMS radius)
        if (if3d == 1):
            w,x,y,z,ux,uy,uz = rem_part(w,x,z,ux,uy,uz,y=y)
        else:
            w,x,z,ux,uy,uz = rem_part(w,x,z,ux,uy,uz)

        Ep = .511*(uz-1.)
        zz        [count] = np.mean(z)
        TotQ      [count] = np.sum(w)
        Emean     [count] = np.mean(Ep)
        Estd      [count] = np.std(Ep)
        beam_widthx[count] = np.std(x)
        beam_widthux[count] = np.std(ux)
        if (if3d == 1):
            beam_widthy[count] = np.std(y)
            beam_widthuy[count] = np.std(uy)
        if TotQ[count]!=0.0:
            cent_slice = np.abs(z-np.average(z,weights=w))<np.std(z)/30.
            beam_div  [count] = np.std(ux/uz)
            if (if3d == 1):
                beam_widthy_center[count] = np.std(y[cent_slice])
                lEMIT_X, lEMIT_Y, lW_slice, lGAMMA = my_emittance(x, z, ux,
                                uy, uz, sliced=True, nslices=nslices, y=y)
                EMITY[count,:] = lEMIT_Y
            else:
                lEMIT_X, lW_slice, lGAMMA = my_emittance(x, z, ux,
                                uy, uz, sliced=True, nslices=nslices)
            EMITX[count,:] = lEMIT_X
            W_slice[count,:] = lW_slice
            emittancex [count] = np.average(lEMIT_X, weights=lW_slice)
            if (if3d == 1):
                emittancepx[count]= emittance_from_coord(x, ux, uy, w, y)[0]
                emittancey [count] = np.average(lEMIT_Y, weights=lW_slice)
                emittancepy[count]= emittance_from_coord(x, ux, uy, w, y)[1]
            else:
                emittancepx[count]= emittance_from_coord(x, ux, uy, w)
        Emean[Emean==0.] = np.nan
        Estd[Estd==0.] = np.nan
        beam_widthx[beam_widthx==0.] = np.nan
        beam_widthux[beam_widthux==0.] = np.nan
        EMITX[EMITX==0.] = np.nan
        emittancex[emittancex==0.] = np.nan
        emittancepx[emittancepx==0.] = np.nan
        beam_div[beam_div==0.] = np.nan
        if (if3d == 1):
            beam_widthy[beam_widthy==0.] = np.nan
            beam_widthuy[beam_widthuy==0.] = np.nan
            EMITY[EMITY==0.] = np.nan
            emittancey[emittancey==0.] = np.nan
            emittancepy[emittancepy==0.] = np.nan
    if (if3d == 1):
        results = (zz, TotQ, TotP, Emean, Estd, emittancepx, emittancepy,
                emittancex, emittancey, beam_widthx, beam_widthy, beam_widthux,
                beam_widthuy)
    else:
        results = (zz, TotQ, TotP, Emean, Estd, emittancepx, emittancex,
                beam_widthx, beam_widthux)
    return results

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
    results = read_sim(sim)
    if (len(results) == 13): #same as if3d == 1, but without issue in github
         (zz, TotQ, TotP, Emean, Estd, emittancepx, emittancepy, emittancex,
          emittancey, beam_widthx, beam_widthy, beam_widthux,
          beam_widthuy) = results
    else:
        (zz, TotQ, TotP, Emean, Estd, emittancepx, emittancex, beam_widthx,
         beam_widthux) = results
    d_z     [sim] = zz
    d_totq  [sim] = TotQ*scc.e*1e9 # to be in nC
    d_totp  [sim] = TotP
    d_emean [sim] = Emean*1.e-3 # to be in GeV
    d_estd  [sim] = Estd/Emean*1.e2 # to be in %
    d_emitpx[sim] = emittancepx*1.e6 # to be in um
    d_emitx [sim] = emittancex*1.e6 # to be in um
    d_bwdx  [sim] = beam_widthx*1.e6 # to be in um
    d_bwdux  [sim] = beam_widthux
    if (if3d == 1):
        d_emitpy[sim] = emittancepy*1.e6 # to be in um
        d_emity [sim] = emittancey*1.e6 # to be in um
        d_bwdy  [sim] = beam_widthy*1.e6 # to be in um
        d_bwduy  [sim] = beam_widthuy

# Properties as they will be stored into the HDF5 file
a_list=['z','totq','totp','emean','estd', 'emitpx', 'emitx', 'xsig','uxsig',
        'emitpy','emity','ysig','uysig']
l_list=['z','Total charge in '+str(d_rms)+' RMS','Total particles in box',
        'Mean energy',
        'Relative energy spread','Projected $\epsilon$ in x',
        'Average slice $\epsilon$ in x',
        'Beam width in x','Beam width in ux','Projected $\epsilon$ in y',
        'Average slice $\epsilon$ in y','Beam width in y','Beam width in uy']
u_list=['m','nC','# macro-particles','GeV','%','$\mu$m','$\mu$m','$\mu$m','m/s',
        '$\mu$m','$\mu$m','$\mu$m','m/s']

d_att={}
for sim in sim_list:
    d_att[sim]={}
    for i in range(len(a_list)):
        d_att[sim][a_list[i]]={}
        d_att[sim][a_list[i]]['label']=l_list[i]
        d_att[sim][a_list[i]]['units']=u_list[i]
    d_att[sim][a_list[0]]['data']=d_z[sim]
    d_att[sim][a_list[1]]['data']=d_totq[sim]
    d_att[sim][a_list[2]]['data']=d_totp[sim]
    d_att[sim][a_list[3]]['data']=d_emean[sim]
    d_att[sim][a_list[4]]['data']=d_estd[sim]
    d_att[sim][a_list[5]]['data']=d_emitpx[sim]
    d_att[sim][a_list[6]]['data']=d_emitx[sim]
    d_att[sim][a_list[7]]['data']=d_bwdx[sim]
    d_att[sim][a_list[8]]['data']=d_bwdux[sim]
    if (if3d == 1):
        d_att[sim][a_list[9]]['data']=d_emitpy[sim]
        d_att[sim][a_list[10]]['data']=d_emity[sim]
        d_att[sim][a_list[11]]['data']=d_bwdy[sim]
        d_att[sim][a_list[12]]['data']=d_bwduy[sim]

def new_files(pd,dh5):
    # Creating output files names
    newf_name=[]
    for isims in sim_list:
        new_sub_dir=pd+isims+'/lab_frame_data/'+dh5
        if os.path.isdir(new_sub_dir) == False:
            print('Created directory: ', new_sub_dir)
            os.mkdir(new_sub_dir)
        newf_name.append(new_sub_dir+'/'+dh5+'-'+isims+'.hdf5')
        print(newf_name[-1])
    return newf_name

fhdf5_list= new_files(path_scan,dir_hdf5)

# Function to save the information in files
def save_hdf5(fn,att,fatt):
    f = h5py.File(fn, "w")
    for a in att:
        dset = f.create_dataset(a, tuple(fatt[a]['data'].shape), dtype='f')
        dset[:] = fatt[a]['data'][:]
        dset.attrs['label']=fatt[a]['label']
        dset.attrs['units']=fatt[a]['units']
        dset.attrs['range']=[np.amin(fatt[a]['data']),np.amax(fatt[a]['data'])]
    f.close()
    file_created='Stored dataset in file: '+fn
    return file_created


for isims in range(nsims):
    if (if3d == 1):
        print(save_hdf5(fhdf5_list[isims],a_list,d_att[sim_list[isims]]))
    else:
        print(save_hdf5(fhdf5_list[isims],a_list[:-4],d_att[sim_list[isims]]))

exit
