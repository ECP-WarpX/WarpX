#include "Particles/Gather/GetExternalFields.H"

#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "WarpX.H"

#include <AMReX_Vector.H>

#include <string>


GetExternalEField::GetExternalEField (const WarpXParIter& a_pti, int a_offset) noexcept
{
    auto& warpx = WarpX::GetInstance();
    auto& mypc = warpx.GetPartContainer();
    if (mypc.m_E_ext_particle_s=="constant" || mypc.m_E_ext_particle_s=="default")
    {
        m_type = Constant;
        m_field_value[0] = mypc.m_E_external_particle[0];
        m_field_value[1] = mypc.m_E_external_particle[1];
        m_field_value[2] = mypc.m_E_external_particle[2];
    }
    else if (mypc.m_E_ext_particle_s=="parse_e_ext_particle_function")
    {
        m_type = Parser;
        m_time = warpx.gett_new(a_pti.GetLevel());
        m_get_position = GetParticlePosition(a_pti, a_offset);
        m_xfield_partparser = getParser(mypc.m_Ex_particle_parser);
        m_yfield_partparser = getParser(mypc.m_Ey_particle_parser);
        m_zfield_partparser = getParser(mypc.m_Ez_particle_parser);
    }
    else if (mypc.m_E_ext_particle_s=="repeated_plasma_lens")
    {
        m_type = RepeatedPlasmaLens;
        m_dt = warpx.getdt(a_pti.GetLevel());
        m_get_position = GetParticlePosition(a_pti, a_offset);
        auto& attribs = a_pti.GetAttribs();
        m_ux = attribs[PIdx::uy].dataPtr() + a_offset;
        m_uy = attribs[PIdx::uz].dataPtr() + a_offset;
        m_uz = attribs[PIdx::uz].dataPtr() + a_offset;
        m_repeated_plasma_lens_period = mypc.m_repeated_plasma_lens_period;
        int const n_lenses = static_cast<int>(mypc.m_repeated_plasma_lens_starts.size());
        m_repeated_plasma_lens_starts.resize(n_lenses);
        m_repeated_plasma_lens_lengths.resize(n_lenses);
        m_repeated_plasma_lens_strengths.resize(n_lenses);
        for (int i=0 ; i < n_lenses ; i++) {
            m_repeated_plasma_lens_starts[i] = mypc.m_repeated_plasma_lens_starts[i];
            m_repeated_plasma_lens_lengths[i] = mypc.m_repeated_plasma_lens_lengths[i];
            m_repeated_plasma_lens_strengths[i] = mypc.m_repeated_plasma_lens_strengths[i];
        }
    }
}

GetExternalBField::GetExternalBField (const WarpXParIter& a_pti, int a_offset) noexcept
{
    auto& warpx = WarpX::GetInstance();
    auto& mypc = warpx.GetPartContainer();
    if (mypc.m_B_ext_particle_s=="constant" || mypc.m_B_ext_particle_s=="default")
    {
        m_type = Constant;
        m_field_value[0] = mypc.m_B_external_particle[0];
        m_field_value[1] = mypc.m_B_external_particle[1];
        m_field_value[2] = mypc.m_B_external_particle[2];
    }
    else if (mypc.m_B_ext_particle_s=="parse_b_ext_particle_function")
    {
        m_type = Parser;
        m_time = warpx.gett_new(a_pti.GetLevel());
        m_get_position = GetParticlePosition(a_pti, a_offset);
        m_xfield_partparser = getParser(mypc.m_Bx_particle_parser);
        m_yfield_partparser = getParser(mypc.m_By_particle_parser);
        m_zfield_partparser = getParser(mypc.m_Bz_particle_parser);
    }
}
