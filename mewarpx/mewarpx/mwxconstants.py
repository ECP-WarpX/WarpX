"""File to keep constants not already in picmi.py that are needed in mewarpx."""
from pywarpx.picmi import constants

# CONSTANTS - SI
kb_J = constants.kb # J/K
torr_SI = constants.torr_SI # 1 torr in Pa
erg_SI = 1e-7 # 1 erg in J

# CONSTANTS - CGS
kb_cgs = kb_J / erg_SI # erg/K
torr_cgs = torr_SI * 10 # 1 torr in dyne/cm^2
