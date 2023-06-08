#Laser enevelope model

void WarpX::AllocateLaserEnvelope (int lev, const BoxArray& ba, const DistributionMapping& dm, const IntVect& ngA, const bool aux_is_nodal)
{
    // Declare nodal flags
    IntVect A_nodal_flag
    // Set nodal flags
    #if   defined(WARPX_DIM_1D_Z)
    A_nodal_flag = IntVect(0)
    #elif   defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    A_nodal_flag = IntVect(0,0)
    #elif defined(WARPX_DIM_3D)
    A_nodal_flag = IntVect(0,0,0)
    #endif

    if (electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic)
    {
    A_nodal_flag = IntVect::TheNodeVector();
    }

    const int A_ncomps = (WarpX::do_multi_J) ? ncomps : 2*ncomps;
    AllocInitMultiFab(A_potential_vector[lev], amrex::convert(ba, A_nodal_flag), dm, A_ncomps, ngA, tag("A_fp"), 0.0_rt);
     
}