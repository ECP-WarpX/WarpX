#ifndef WARPX_MAGNETICMATERIAL_H_
#define WARPX_MAGNETICMATERIAL_H_

#include <AMReX_REAL.H>

class MagneticMaterial
{
public:

    virtual void UpdateHandM (
              amrex::Real& H_x     ,       amrex::Real& H_y     ,       amrex::Real& H_z    ,
              amrex::Real& M_x     ,       amrex::Real& M_y     ,       amrex::Real& M_z    ,
        const amrex::Real& B_next_x, const amrex::Real& B_next_y, const amrex::Real& B_next_z) const = 0;

};

#endif /* WARPX_MAGNETICMATERIAL_H_ */
