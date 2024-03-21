
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron, c, pi, m_e, femto, e, milli, eV, physical_constants, alpha, nano
from matplotlib import use, cm
import matplotlib.colors
import pandas as pd 

r_e = physical_constants['classical electron radius'][0]
my_dpi=300
sigmaz = 10*nano

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(2000./my_dpi, 1000./my_dpi), dpi=my_dpi)

rdir= './diags/reducedfiles/'

df_cr = pd.read_csv(f'{rdir}'+'ColliderRelevant_beam1_beam2.txt', sep=" ", header=0)
df_pn = pd.read_csv(f'{rdir}'+'ParticleNumber.txt', sep=" ", header=0)


times = df_cr[[col for col in df_cr.columns if f']time' in col]].to_numpy()
steps = df_cr[[col for col in df_cr.columns if f']step' in col]].to_numpy()

x = df_cr[[col for col in df_cr.columns if f']dL_dt' in col]].to_numpy()
coll_index = np.argmax(x)
coll_time = times[coll_index]

# number of photons per beam particle 
np1 = df_pn[[col for col in df_pn.columns if f']pho1_weight' in col]].to_numpy()
np2 = df_pn[[col for col in df_pn.columns if f']pho2_weight' in col]].to_numpy()
Ne = df_pn[[col for col in df_pn.columns if f']beam1_weight' in col]].to_numpy()[0]
Np = df_pn[[col for col in df_pn.columns if f']beam2_weight' in col]].to_numpy()[0]

ax[0].plot((times-coll_time)/(sigmaz/c), (np1+np2)/(Ne+Np), lw=2)
ax[0].set_title(r'photon number/beam particle')

# number of NLBW particles per beam particle 
e1 = df_pn[[col for col in df_pn.columns if f']ele1_weight' in col]].to_numpy()
e2 = df_pn[[col for col in df_pn.columns if f']ele2_weight' in col]].to_numpy()

ax[1].plot((times-coll_time)/(sigmaz/c), (e1+e2)/(Ne+Np), lw=2)
ax[1].set_title(r'NLBW particles/beam particle')
   
for a in ax.reshape(-1):
    a.set_xlabel(r'time [$\sigma_z/c$]')
image_file_name ='reduced.png'
plt.tight_layout()
plt.savefig(image_file_name,dpi=300, bbox_inches='tight') 
plt.close("all")
