/* Copyright 2021 Lorenzo Giacomel, Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_

#include "ReducedDiags.H"
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
     * This function computes the value of Ex, Ey, Ez, Bx, By, Bz and at a given point
     *
     * @param[in] step current time step
     */
    virtual void ComputeDiags(int step) override final;

    static constexpr int noutputs = 6; // probing Ex,Ey,Ez,Bx,By,Bz
    static constexpr int index_Ex = 0;
    static constexpr int index_Ey = 1;
    static constexpr int index_Ez = 2;
    static constexpr int index_Bx = 3;
    static constexpr int index_By = 4;
    static constexpr int index_Bz = 5;

    bool m_probe_in_domain;
    bool m_probe_in_domain_lev_0 = false;

private:
#if (AMREX_SPACEDIM == 2)
    amrex::Real x_probe, z_probe;
#elif (AMREX_SPACEDIM == 3)
    amrex::Real x_probe, y_probe, z_probe;
#endif

    int interp_order = 1;
    bool raw_fields = false;

    virtual void WriteToFile (int step) const override;

};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
