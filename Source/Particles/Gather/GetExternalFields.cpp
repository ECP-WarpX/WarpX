#include "WarpX.H"
#include "Particles/Gather/GetExternalFields.H"

GetExternalEBField::GetExternalEBField (const WarpXParIter& a_pti, int a_offset) noexcept
{
    auto& warpx = WarpX::GetInstance();
    auto& mypc  = warpx.GetPartContainer();

    // E field
    if (mypc.m_E_ext_particle_s=="constant" || mypc.m_E_ext_particle_s=="default")
    {
        m_E_type = Constant;
        m_E_field_value[0] = mypc.m_E_external_particle[0];
        m_E_field_value[1] = mypc.m_E_external_particle[1];
        m_E_field_value[2] = mypc.m_E_external_particle[2];
    }
    else if (mypc.m_E_ext_particle_s=="parse_e_ext_particle_function")
    {
        m_E_type = Parser;
        m_time = warpx.gett_new(a_pti.GetLevel());
        m_get_position = GetParticlePosition(a_pti, a_offset);
        m_Ex_field_partparser = mypc.m_Ex_particle_parser.get();
        m_Ey_field_partparser = mypc.m_Ey_particle_parser.get();
        m_Ez_field_partparser = mypc.m_Ez_particle_parser.get();
    }

    // B field
    if (mypc.m_B_ext_particle_s=="constant" || mypc.m_B_ext_particle_s=="default")
    {
        m_B_type = Constant;
        m_B_field_value[0] = mypc.m_B_external_particle[0];
        m_B_field_value[1] = mypc.m_B_external_particle[1];
        m_B_field_value[2] = mypc.m_B_external_particle[2];
    }
    else if (mypc.m_B_ext_particle_s=="parse_b_ext_particle_function")
    {
        m_B_type = Parser;
        m_time = warpx.gett_new(a_pti.GetLevel());
        m_get_position = GetParticlePosition(a_pti, a_offset);
        m_Bx_field_partparser = mypc.m_Bx_particle_parser.get();
        m_By_field_partparser = mypc.m_By_particle_parser.get();
        m_Bz_field_partparser = mypc.m_Bz_particle_parser.get();
    }
}
