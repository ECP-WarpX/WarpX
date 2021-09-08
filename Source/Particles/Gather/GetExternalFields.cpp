#include "Particles/Gather/GetExternalFields.H"

#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "WarpX.H"

#include <AMReX_Vector.H>

#include <string>

using namespace amrex::literals;

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
        m_xfield_partparser = mypc.m_Ex_particle_parser->compile<4>();
        m_yfield_partparser = mypc.m_Ey_particle_parser->compile<4>();
        m_zfield_partparser = mypc.m_Ez_particle_parser->compile<4>();
    }
    else if (mypc.m_E_ext_particle_s=="repeated_plasma_lens")
    {
        m_type = RepeatedPlasmaLens;
        m_time = warpx.gett_new(a_pti.GetLevel());
        m_gamma_boost = WarpX::gamma_boost;
        m_uz_boost = std::sqrt(WarpX::gamma_boost*WarpX::gamma_boost - 1._rt)*PhysConst::c;
        m_lens_is_electric = 1;
        m_dt = warpx.getdt(a_pti.GetLevel());
        m_get_position = GetParticlePosition(a_pti, a_offset);
        auto& attribs = a_pti.GetAttribs();
        m_ux = attribs[PIdx::ux].dataPtr() + a_offset;
        m_uy = attribs[PIdx::uy].dataPtr() + a_offset;
        m_uz = attribs[PIdx::uz].dataPtr() + a_offset;
        m_repeated_plasma_lens_period = mypc.m_repeated_plasma_lens_period;
        m_n_lenses = static_cast<int>(mypc.h_repeated_plasma_lens_starts.size());
        m_repeated_plasma_lens_starts = mypc.d_repeated_plasma_lens_starts.data();
        m_repeated_plasma_lens_lengths = mypc.d_repeated_plasma_lens_lengths.data();
        m_repeated_plasma_lens_strengths_E = mypc.d_repeated_plasma_lens_strengths_E.data();
        m_repeated_plasma_lens_strengths_B = mypc.d_repeated_plasma_lens_strengths_B.data();
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
        m_xfield_partparser = mypc.m_Bx_particle_parser->compile<4>();
        m_yfield_partparser = mypc.m_By_particle_parser->compile<4>();
        m_zfield_partparser = mypc.m_Bz_particle_parser->compile<4>();
    }
    else if (mypc.m_B_ext_particle_s=="repeated_plasma_lens")
    {
        m_type = RepeatedPlasmaLens;
        m_time = warpx.gett_new(a_pti.GetLevel());
        m_gamma_boost = WarpX::gamma_boost;
        m_uz_boost = std::sqrt(WarpX::gamma_boost*WarpX::gamma_boost - 1._rt)*PhysConst::c;
        m_lens_is_electric = 0;
        m_dt = warpx.getdt(a_pti.GetLevel());
        m_get_position = GetParticlePosition(a_pti, a_offset);
        auto& attribs = a_pti.GetAttribs();
        m_ux = attribs[PIdx::ux].dataPtr() + a_offset;
        m_uy = attribs[PIdx::uy].dataPtr() + a_offset;
        m_uz = attribs[PIdx::uz].dataPtr() + a_offset;
        m_repeated_plasma_lens_period = mypc.m_repeated_plasma_lens_period;
        m_n_lenses = static_cast<int>(mypc.h_repeated_plasma_lens_starts.size());
        m_repeated_plasma_lens_starts = mypc.d_repeated_plasma_lens_starts.data();
        m_repeated_plasma_lens_lengths = mypc.d_repeated_plasma_lens_lengths.data();
        m_repeated_plasma_lens_strengths_E = mypc.d_repeated_plasma_lens_strengths_E.data();
        m_repeated_plasma_lens_strengths_B = mypc.d_repeated_plasma_lens_strengths_B.data();
    }
}
