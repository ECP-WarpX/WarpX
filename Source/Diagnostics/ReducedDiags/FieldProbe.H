/* Copyright 2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_

#include "ReducedDiags.H"

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
     * This function computes the value of Ex, Ey, Ez, |E|, Bx, By, Bz and |B| at a given point
     *
     * @param[in] step current time step
     */
    virtual void ComputeDiags(int step) override final;

private:
    amrex::ParticleReal x_probe, y_probe, z_probe;
    int interp_order;
    bool raw_fields;
};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDMAXIMUM_H_
