from pywarpx import picmi

def init_mcc(electrons, ions, N_INERT, T_INERT):
    """Initializes MCC objects."""
    cross_sec_direc = '../../../warpx-data/MCC_cross_sections/He/'
    mcc_electrons = picmi.MCCCollisions(
        name='coll_elec',
        species=electrons,
        background_density=N_INERT,
        background_temperature=T_INERT,
        background_mass=ions.mass,
        scattering_processes={
            'elastic' : {
                'cross_section' : cross_sec_direc+'electron_scattering.dat'
            },
            'excitation1' : {
                'cross_section': cross_sec_direc+'excitation_1.dat',
                'energy' : 19.82
            },
            'excitation2' : {
                'cross_section': cross_sec_direc+'excitation_2.dat',
                'energy' : 20.61
            },
            'ionization' : {
                'cross_section' : cross_sec_direc+'ionization.dat',
                'energy' : 24.55,
                'species' : ions
            },
        }
    )

    mcc_ions = picmi.MCCCollisions(
        name='coll_ion',
        species=ions,
        background_density=N_INERT,
        background_temperature=T_INERT,
        scattering_processes={
            'elastic' : {
                'cross_section' : cross_sec_direc+'ion_scattering.dat'
            },
            'back' : {
                'cross_section' : cross_sec_direc+'ion_back_scatter.dat'
            },
            # 'charge_exchange' : {
            #    'cross_section' : cross_sec_direc+'charge_exchange.dat'
            # }
        }
    )

    return mcc_electrons, mcc_ions