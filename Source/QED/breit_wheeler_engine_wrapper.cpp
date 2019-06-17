//This file provides a wrapper aroud the breit_wheeler engine
//provided by PICSAR QED

#include "breit_wheeler_engine_wrapper.h"

warpx_breit_wheeler_engine::warpx_breit_wheeler_engine():
    breit_wheeler_engine{std::move(amrex_rng_wrapper{})}
{}
