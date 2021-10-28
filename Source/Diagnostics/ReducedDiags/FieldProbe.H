/* Copyright 2021 Lorenzo Giacomel, Neil Zaim, Yinjian Zhao
 *                Tiberius Rheaume, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_

#include "ReducedDiags.H"
#include "FieldProbeParticleContainer.H"

#include <AMReX.H>
#include <AMReX_Vector.H>

#include <unordered_map>
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

    //! noutputs is 7 (Ex, Ey, Ez, Bx, By, Bz, S)
    static constexpr int noutputs = FieldProbePIdx::nattribs;

    //! u_map maps FieldProbe observables to units
    std::unordered_map< int, std::string > u_map;

    //! this boolean is used to disqualify out of bounds particles
    bool m_probe_in_domain;
    bool m_probe_in_domain_lev_0 = false;

    //! this is the particle container in which probe particles are stored
    FieldProbeParticleContainer m_probe;

private:
    amrex::Real x_probe, y_probe, z_probe;

    //! if  true, integrate values over time instead of probing instantaneous values
    bool m_field_probe_integrate = false;
    int interp_order = 1;
    bool raw_fields = false;

    /**
     * Built-in function in ReducedDiags to write out test data
     */
    virtual void WriteToFile (int step) const override;

};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
