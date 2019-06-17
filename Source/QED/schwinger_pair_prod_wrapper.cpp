//This file provides a wrapper aroud the schwinger pair production engine
//provided by PICSAR QED

#include "schwinger_pair_prod_wrapper.h"

warpx_schwinger_pair_engine::warpx_schwinger_pair_engine():
    schwinger_pair_engine{std::move(amrex_rng_wrapper{})}
{}
