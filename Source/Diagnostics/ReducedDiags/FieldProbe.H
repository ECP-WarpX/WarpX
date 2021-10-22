/* Copyright 2021 Lorenzo Giacomel, Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_

#include "ReducedDiags.H"
#include "FieldProbeParticleContainer.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_Vector.H>

#include <string>

/**
 *  This class mainly contains a function that computes the value of each component
 * of the EM field at a given point
 */
class FieldProbe : public ReducedDiags
{
public:

    /**
     * constructor
     * @param[in] rd_name reduced diags names
     */
    FieldProbe(std::string rd_name);

    /** 
     * This function assins test/data particles to constructed environemnt
     */
    void InitData() override final;

    /**
     * This function computes the value of Ex, Ey, Ez, Bx, By, Bz and at a given point
     *
     * @param[in] step current time step
     */
    void ComputeDiags(int step) override final;

    /*
     * Define constants used throughout FieldProbe
     */

    /**
     * np is number of particles to be used
     */
    static constexpr int np = 1;

    /**
     *noutputs is 7 (Ex, Ey, Ez, Bx, By, Bz, S)
     */
    static constexpr int noutputs = static_cast<int>(ParticleVal::nattribs);

    /**
     * this boolean is used to disqualify out of bounds particles
     */
    bool m_probe_in_domain;
    bool m_probe_in_domain_lev_0 = false;

    /**
    * initializes particle container m_probe where particles will be stored
    */
    FieldProbeParticleContainer m_probe;

private:
    amrex::Real x_probe, y_probe, z_probe;

    //! for field_probe_integrate = 0, will not integrate values
    int field_probe_integrate = 0;
    int interp_order = 1;
    bool raw_fields = false;

    //vacper is vacuum permeability
    amrex::Real vacper = PhysConst::mu0;

    /**
     * Built-in function in ReducedDiags to write out test data
     */
    virtual void WriteToFile (int step) const override;

};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
