from matplotlib import cm, rcParams, use
import matplotlib.colors
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import openpmd_api as io
from scipy.constants import (c, centi, e, eV, femto, m_e, micron, milli,
                             physical_constants, pi, pico)

use('AGG')
rcParams.update({'font.size': 8})
MeV=1e6*eV
r_e = physical_constants["classical electron radius"][0]
my_dpi = 300.


def print_full():
    warpx_dir='./diags/particles/'


    for ts in range(400):
        series = io.Series(warpx_dir+"/openpmd_%06d.bp" % ts,io.Access.read_only)
        j = series.iterations[ts]
        time = j.time

        ps = j.particles["photonA"]
        x = ps["position"]["x"].load_chunk()
        z = ps["position"]["z"].load_chunk()
        #px = ps["momentum"]["x"].load_chunk()
        #py = ps["momentum"]["y"].load_chunk()
        #pz = ps["momentum"]["z"].load_chunk()
        w = ps["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk()
        series.flush()
        print(w,x,z)

        ps = j.particles["photonB"]
        x = ps["position"]["x"].load_chunk()
        z = ps["position"]["z"].load_chunk()
        #px = ps["momentum"]["x"].load_chunk()
        #py = ps["momentum"]["y"].load_chunk()
        #pz = ps["momentum"]["z"].load_chunk()
        w = ps["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk()
        #charge = electrons["charge"][io.Mesh_Record_Component.SCALAR]
        series.flush()
        print(w,x,z)

        ps = j.particles["electron"]
        x = ps["position"]["x"].load_chunk()
        z = ps["position"]["z"].load_chunk()
        #px = ps["momentum"]["x"].load_chunk()
        #py = ps["momentum"]["y"].load_chunk()
        #pz = ps["momentum"]["z"].load_chunk()
        w = ps["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk()
        #charge = electrons["charge"][io.Mesh_Record_Component.SCALAR]
        series.flush()
        print(w,x,z)

        ps = j.particles["positron"]
        x = ps["position"]["x"].load_chunk()
        z = ps["position"]["z"].load_chunk()
        #px = ps["momentum"]["x"].load_chunk()
        #py = ps["momentum"]["y"].load_chunk()
        #pz = ps["momentum"]["z"].load_chunk()
        w = ps["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk()
        #charge = electrons["charge"][io.Mesh_Record_Component.SCALAR]
        series.flush()
        print(w,x,z)

        print('ooooooooooooooooooooooooooooooooooooooooo')
        del series


def plot_one_panel(ax,data,cmap,title):
    im=ax.imshow(data, cmap=cmap, aspect='equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.001)
    fig.colorbar(im, cax=cax, orientation='vertical')
    ax.set_title(title)


def cross_section(E1_lab, E2_lab, theta):
    s = E1_lab*E2_lab/(2.*m_e**2*c**4)*(1.-np.cos(theta))
    beta = np.sqrt(1.-1./s)
    factor1 = 0.5*pi*r_e**2*(1.-beta**2)
    term1 = (3.-beta**4)*np.log((1.+beta)/(1.-beta))
    term2 = 2*beta*(beta**2-2.)
    factor2 = term1 + term2
    sigma = factor1*factor2
    return sigma

def plot_cross_section():
    emin = 0.5*m_e*c**2
    emax = 6*m_e*c**2
    energy = np.linspace(emin, emax, 1000)
    E1_lab, E2_lab = np.meshgrid(energy, energy)
    theta = pi
    CS = cross_section(E1_lab, E2_lab, theta)
    CS2 = CS[~np.isnan(CS)]
    Elab2 = E1_lab[~np.isnan(CS)]
    amax = np.unravel_index(CS2.argmax(), CS2.shape)
    print('energy for max cross section', Elab2[amax]/(m_e*c**2))
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(1000./my_dpi, 1000./my_dpi), dpi=my_dpi)
    divider = make_axes_locatable(ax)
    im=ax.imshow(CS/r_e**2, aspect='equal', cmap = 'rainbow', origin='lower', extent=[emin/(m_e*c**2), emax/(m_e*c**2), emin/(m_e*c**2), emax/(m_e*c**2)])
    cax = divider.append_axes('right', size='2%', pad=0.1)
    fig.colorbar(im, cax=cax, orientation='vertical', label=r'cross section')
    plt.tight_layout()
    plt.savefig('cs.png',dpi=my_dpi)
    plt.close("all")

def plot_reduced():
    fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(4000./my_dpi, 2000./my_dpi), dpi=my_dpi, sharex=True)

    mydir= './diags/reducedfiles/'
    lw = 4

    # number of photonA
    x = np.loadtxt(mydir+'ParticleNumber.txt')
    ax[0][0].plot(x[:,1], x[:,3], lw=lw)
    ax[0][0].set_title(r'macro photonA number')

    # number of photonB
    ax[0][1].plot(x[:,1], x[:,4], lw=lw)
    ax[0][1].set_title(r'macro photonB number')

    # number of electrons
    ax[0][2].plot(x[:,1], x[:,5], lw=lw)
    ax[0][2].set_title(r'macro electrons')

    # number of positrons
    x = np.loadtxt(mydir+'ParticleNumber.txt')
    ax[0][3].plot(x[:,1], x[:,6], lw=lw)
    ax[0][3].set_title('macro positrons')

    # max weight photonA
    x = np.loadtxt(mydir+'ParticleExtrema_photonA.txt')
    ax[1][0].plot(x[:,1], x[:,17], lw=lw)
    ax[1][0].set_title('max weight photonA')

    # max weight photonB
    x = np.loadtxt(mydir+'ParticleExtrema_photonB.txt')
    ax[1][1].plot(x[:,1], x[:,17], lw=lw)
    ax[1][1].set_title('max weight photonB')

    # max weight electron
    x = np.loadtxt(mydir+'ParticleExtrema_electron.txt')
    ax[1][2].plot(x[:,1], x[:,17], lw=lw)
    ax[1][2].set_title('max weight electron')

    # max weight positron
    x = np.loadtxt(mydir+'ParticleExtrema_positron.txt')
    ax[1][3].plot(x[:,1], x[:,17], lw=lw)
    ax[1][3].set_title('max weight positron')

    plt.tight_layout()
    plt.savefig('reduced.png', dpi=my_dpi, bbox_inches='tight')
    plt.close("all")


def main():
    warpx_used_inputs = open('./warpx_used_inputs', 'r').read()
    #if re.search('photonA.momentum_function_uz', warpx_used_inputs):
    uz1_lab, uz2_lab = 2., -2.
    p1_lab, p2_lab = uz1_lab * m_e * c, uz2_lab * m_e * c
    E1_lab, E2_lab = np.abs(p1_lab) * c, np.abs(p2_lab) * c
    print(p1_lab, p2_lab)
    print(E1_lab*E2_lab > m_e**2 * c**4)
    plot_cross_section()
    #with open('./warpx_used_inputs', 'rt') as f:
    #    lines = f.readlines()
    #    for line in lines:
    #        if 'photonA.momentum_function_uz' in line:
    #            print(line)
    #            E1_lab = re.findall(r'[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?', line)[0]
    #        if 'photonB.momentum_function_uz' in line:
    #            print(line)
    #            E2_lab = re.findall(r'[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?', line)[0]
    #print(E1_lab, E2_lab)
    print(cross_section(E1_lab, E2_lab, pi))
    plot_reduced()
    print_full()

if __name__ == "__main__":
    main()
