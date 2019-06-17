#ifndef WARPX_schwinger_pair_prod_wrapper_h_
#define WARPX_schwinger_pair_prod_wrapper_h_

//This file provides a wrapper aroud the schwinger pair production engine
//provided by PICSAR QED

#include "warpx_qed_commons.h"
#include "schwinger_pair_engine.hpp"

#include "amrex_rng_wrapper.h"

class warpx_schwinger_pair_engine:
  public picsar::multi_physics::schwinger_pair_engine<amrex::Real, amrex_rng_wrapper>
{
  public:
    warpx_schwinger_pair_engine();
};

#endif //WARPX_schwinger_pair_prod_wrapper_h_
