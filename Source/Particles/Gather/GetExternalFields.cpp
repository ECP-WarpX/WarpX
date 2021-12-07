#include "Particles/Gather/GetExternalFields.H"

#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "WarpX.H"

#include <AMReX_Vector.H>

#include <string>

using namespace amrex::literals;

GetExternalEBField::GetExternalEBField (const WarpXParIter& a_pti, int a_offset) noexcept
{
    auto& warpx = WarpX::GetInstance();
    auto& mypc = warpx.GetPartContainer();

    m_gamma_boost = WarpX::gamma_boost;
    m_uz_boost = std::sqrt(WarpX::gamma_boost*WarpX::gamma_boost - 1._rt)*PhysConst::c;

    m_Etype = Unknown;
    m_Btype = Unknown;

    if (mypc.m_E_ext_particle_s == "none") m_Etype = None;
    if (mypc.m_B_ext_particle_s == "none") m_Btype = None;

    if (mypc.m_E_ext_particle_s == "constant")
    {
        m_Etype = Constant;
        m_Efield_value[0] = mypc.m_E_external_particle[0];
        m_Efield_value[1] = mypc.m_E_external_particle[1];
        m_Efield_value[2] = mypc.m_E_external_particle[2];
    }

    if (mypc.m_B_ext_particle_s == "constant")
    {
        m_Btype = Constant;
        m_Bfield_value[0] = mypc.m_B_external_particle[0];
        m_Bfield_value[1] = mypc.m_B_external_particle[1];
        m_Bfield_value[2] = mypc.m_B_external_particle[2];
    }

    if (mypc.m_E_ext_particle_s == "parse_e_ext_particle_function" ||
        mypc.m_B_ext_particle_s == "parse_b_ext_particle_function" ||
        mypc.m_E_ext_particle_s == "repeated_plasma_lens" ||
        mypc.m_B_ext_particle_s == "repeated_plasma_lens")
    {
        m_time = warpx.gett_new(a_pti.GetLevel());
        m_get_position = GetParticlePosition(a_pti, a_offset);
    }

    if (mypc.m_E_ext_particle_s == "parse_e_ext_particle_function")
    {
        m_Etype = Parser;
        m_Exfield_partparser = mypc.m_Ex_particle_parser->compile<4>();
        m_Eyfield_partparser = mypc.m_Ey_particle_parser->compile<4>();
        m_Ezfield_partparser = mypc.m_Ez_particle_parser->compile<4>();
    }

    if (mypc.m_B_ext_particle_s == "parse_b_ext_particle_function")
    {
        m_Btype = Parser;
        m_Bxfield_partparser = mypc.m_Bx_particle_parser->compile<4>();
        m_Byfield_partparser = mypc.m_By_particle_parser->compile<4>();
        m_Bzfield_partparser = mypc.m_Bz_particle_parser->compile<4>();
    }

    if (mypc.m_E_ext_particle_s == "repeated_plasma_lens" ||
        mypc.m_B_ext_particle_s == "repeated_plasma_lens")
    {
        if (mypc.m_E_ext_particle_s == "repeated_plasma_lens") m_Etype = RepeatedPlasmaLens;
        if (mypc.m_B_ext_particle_s == "repeated_plasma_lens") m_Btype = RepeatedPlasmaLens;
        m_dt = warpx.getdt(a_pti.GetLevel());
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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Etype != Unknown, "Unknown E_ext_particle_init_style");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Btype != Unknown, "Unknown B_ext_particle_init_style");

}
