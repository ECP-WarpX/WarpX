#ifndef WARPX_qed_commons_h_
#define WARPX_qed_commons_h_

//Uses SI units
#define PXRMP_WITH_SI_UNITS
#include "qed_commons.h"

#include <AMReX_AmrCore.H>

//A constant used to convert normalized momenta into SI momenta
const amrex::Real qed_me_c =
  picsar::multi_physics::electron_mass * picsar::multi_physics::light_speed;


#endif //WARPX_qed_commons_h_
