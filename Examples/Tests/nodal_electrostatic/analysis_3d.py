import numpy as np

# check that the maximum chi value is small 
fname = 'diags/reducedfiles/ParticleExtrema_beam_p.txt'
chi_max = np.loadtxt(fname)[:,19]
assert(np.sum(chi_max<1.e-8) == len(chi_max))

# check that no photons have been produced 
fname = 'diags/reducedfiles/ParticleNumber.txt'
pho_num = np.loadtxt(fname)[:,7]
assert(pho_num.all()==0.)

    
