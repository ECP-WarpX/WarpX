/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDREDUCED_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDREDUCED_H_

#include "ReducedDiags.H"

/**
 *  This class mainly contains a function that
 *  computes the field energy.
 */
class FieldReduced : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    FieldReduced(const std::string& rd_name, const std::vector<std::string>& field_type);

    /** This funciton computes the field energy (EF).
     *  EF = E eps / 2 + B / mu / 2,
     *  where E is the electric field,
     *  B is the magnetic field,
     *  eps is the vacuum permittivity,
     *  mu is the vacuum permeability. */
    virtual void ComputeDiags(int step) override final;

private:

    /// This bool is set to true if field energy is activated
    bool m_fieldEnergy = false;
    /// This bool is set to true if maximum field is activated
    bool m_maxField = false;
    /// Offset indicating where to write maximum value of fields in m_data
    int m_offset_maxField;

};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDREDUCED_H_
