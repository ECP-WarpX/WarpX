"""File to keep constants not already in picmi.py that are needed in mewarpx."""
from pywarpx.picmi import constants

import numpy as np

# CONSTANTS - SI
h = constants.hbar * 2*np.pi # J s
e = constants.q_e # C
m_e = constants.m_e # kg
kb_J = constants.kb # J/K

# The theoretical value for the Richardson constant, ~120 A/cm^2/K^2
A0 = 4 * np.pi * m_e * kb_J**2 * e / h**3 * 1e-4

torr_SI = constants.torr_SI # 1 torr in Pa
erg_SI = 1e-7 # 1 erg in J

# CONSTANTS - CGS
kb_cgs = kb_J / erg_SI # erg/K
torr_cgs = torr_SI * 10 # 1 torr in dyne/cm^2

kb_eV = kb_J / e # eV/K
