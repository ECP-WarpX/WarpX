#ifndef WARPX_breit_wheeler_engine_wrapper_h_
#define WARPX_breit_wheeler_engine_wrapper_h_

//This file provides a wrapper aroud the breit_wheeler engine
//provided by the standard template library



//BW ENGINE
//#define PXRMP_GPU __host__ __device__
#include "warpx_qed_commons.h"
#include "breit_wheeler_engine.hpp"

#include "amrex_rng_wrapper.h"

class warpx_breit_wheeler_engine:
  public picsar::multi_physics::breit_wheeler_engine<amrex::Real, amrex_rng_wrapper>
{
  public:
    warpx_breit_wheeler_engine();
};

//___________________________________________

#endif //WARPX_breit_wheeler_engine_wrapper_H_
