/* Copyright 2019 Andrew Myers, Axel Huebl, David Grote
 * Luca Fedeli, Maxence Thevenet, Revathi Jambunathan
 * Weiqun Zhang, levinem
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <MultiParticleContainer.H>
#include <WarpX.H>
#include <AMReX_Random.H>
#include <cmath>

using namespace amrex;

void
RigidInjectedParticleContainer::ReadHeader (std::istream& is)
{
    is >> charge >> mass;
    WarpX::GotoNextLine(is);

    int nlevs;
    is >> nlevs;
    WarpX::GotoNextLine(is);

    AMREX_ASSERT(zinject_plane_levels.size() == 0);
    AMREX_ASSERT(done_injecting.size() == 0);

    for (int i = 0; i < nlevs; ++i)
    {
        int zinject_plane_tmp;
        is >> zinject_plane_tmp;
        zinject_plane_levels.push_back(zinject_plane_tmp);
        WarpX::GotoNextLine(is);
    }

    for (int i = 0; i < nlevs; ++i)
    {
        int done_injecting_tmp;
        is >> done_injecting_tmp;
        done_injecting.push_back(done_injecting_tmp);
        WarpX::GotoNextLine(is);
    }
}

void
RigidInjectedParticleContainer::WriteHeader (std::ostream& os) const
{
    // no need to write species_id
    os << charge << " " << mass << "\n";
    int nlevs = zinject_plane_levels.size();
    os << nlevs << "\n";
    for (int i = 0; i < nlevs; ++i)
    {
        os << zinject_plane_levels[i] << "\n";
    }
    for (int i = 0; i < nlevs; ++i)
    {
        os << done_injecting[i] << "\n";
    }
}

void
WarpXParticleContainer::ReadHeader (std::istream& is)
{
    is >> charge >> mass;
    WarpX::GotoNextLine(is);
}

void
WarpXParticleContainer::WriteHeader (std::ostream& os) const
{
    // no need to write species_id
    os << charge << " " << mass << "\n";
}

void
MultiParticleContainer::Checkpoint (const std::string& dir) const
{
    for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
        allcontainers[i]->Checkpoint(dir, species_names[i]);
    }
}

void
MultiParticleContainer::WritePlotFile (const std::string& dir) const
{

    for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
        auto& pc = allcontainers[i];
        if (pc->plot_species) {

            Vector<std::string> real_names;
            Vector<std::string> int_names;
            Vector<int> int_flags;

            real_names.push_back("weight");

            real_names.push_back("momentum_x");
            real_names.push_back("momentum_y");
            real_names.push_back("momentum_z");

            real_names.push_back("Ex");
            real_names.push_back("Ey");
            real_names.push_back("Ez");

            real_names.push_back("Bx");
            real_names.push_back("By");
            real_names.push_back("Bz");

#ifdef WARPX_DIM_RZ
            real_names.push_back("theta");
#endif

            if(pc->do_field_ionization){
                int_names.push_back("ionization_level");
                // int_flags specifies, for each integer attribs, whether it is
                // dumped to plotfiles. So far, ionization_level is the only
                // integer attribs, and it is automatically dumped to plotfiles
                // when ionization is on.
                int_flags.resize(1, 1);
            }

#ifdef WARPX_QED
                if(pc->m_do_qed){
                        real_names.push_back("tau");
                }
#endif

            // Convert momentum to SI
            pc->ConvertUnits(ConvertDirection::WarpX_to_SI);
            // Select particles to be written
            // (by setting IDs of unselected particles to negative numbers)
            pc->SelectParticlesForIO(species_names[i]);
            // real_names contains a list of all particle attributes.
            // pc->plot_flags is 1 or 0, whether quantity is dumped or not.
            pc->WritePlotFile(dir, species_names[i],
                              pc->plot_flags, int_flags,
                              real_names, int_names);
            // Reset negative IDs back to positive
            pc->UnSelectParticlesAfterIO();
            // Convert momentum back to WarpX units
            pc->ConvertUnits(ConvertDirection::SI_to_WarpX);
        }
    }
}

void
MultiParticleContainer::Restart (const std::string& dir)
{
    for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
        allcontainers[i]->Restart(dir, species_names[i]);
    }
}

void
MultiParticleContainer::ReadHeader (std::istream& is)
{
    for (auto& pc : allcontainers) {
        pc->ReadHeader(is);
    }
}

void
MultiParticleContainer::WriteHeader (std::ostream& os) const
{
    for (const auto& pc : allcontainers) {
        pc->WriteHeader(os);
    }
}

// Particle momentum is defined as gamma*velocity, which is neither
// SI mass*gamma*velocity nor normalized gamma*velocity/c.
// This converts momentum to SI units (or vice-versa) to write SI data
// to file.
void
PhysicalParticleContainer::ConvertUnits(ConvertDirection convert_direction)
{
    BL_PROFILE("PPC::ConvertUnits()");

    // Compute conversion factor
    Real factor = 1;
    if (convert_direction == ConvertDirection::WarpX_to_SI){
        factor = mass;
    } else if (convert_direction == ConvertDirection::SI_to_WarpX){
        factor = 1./mass;
    }

    const int nLevels = finestLevel();
    for (int lev=0; lev<=nLevels; lev++){
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // - momenta are stored as a struct of array, in `attribs`
            auto& attribs = pti.GetAttribs();
            ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            ParticleReal* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            ParticleReal* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
            // Loop over the particles and convert momentum
            const long np = pti.numParticles();
            ParallelFor( np,
                [=] AMREX_GPU_DEVICE (long i) {
                    ux[i] *= factor;
                    uy[i] *= factor;
                    uz[i] *= factor;
                }
            );
        }
    }
}

// This function selects particles for IO
// by setting unselected particle IDs to be negative
void
PhysicalParticleContainer::SelectParticlesForIO(std::string species_name)
{
    BL_PROFILE("PPC::SelectParticlesForIO()");

    // plot downsampling type
    std::string type = "default";
    // plot downsampling fraction
    Real fraction = 1.0;
    // flags
    int is_type_given = 0;
    int is_fraction_given = 0;

    // read input
    ParmParse pp(species_name);
    is_type_given = pp.query("plot_downsampling_type",type);
    is_fraction_given = pp.query("plot_downsampling_fraction",fraction);
    int inv_fraction = int(floor(1.0/fraction));

    // assert
    if ( is_type_given == 1 && is_fraction_given == 0 )
    {
        amrex::Abort("Missing a <species>.plot_downdampling_fraction.");
    }
    if ( is_type_given == 0 && is_fraction_given == 1 )
    {
        amrex::Abort("Missing a <species>.plot_downdampling_type.");
    }

    // parser
    auto & mypc = WarpX::GetInstance().GetPartContainer();
    auto nSpecies = mypc.nSpecies();
    auto species_names = mypc.GetSpeciesNames();
    int i_s;
    for ( int i = 0; i < nSpecies; ++i )
    { if ( species_names[i] == species_name ) { i_s = i; } }
    auto & myspc = mypc.GetParticleContainer(i_s);
    ParserWrapper* function_partparser = myspc.m_particle_filter_parser.get();

    const int nLevels = finestLevel();
    for (int lev=0; lev<=nLevels; lev++){
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto pstructs = pti.GetArrayOfStructs()().dataPtr();
            const long np = pti.numParticles();
            // assert
            ParallelFor( np,
                [=] AMREX_GPU_DEVICE (long i) {
                    AMREX_ALWAYS_ASSERT( pstructs[i].id() > 0 );
                }
            );
            // filter
            if ( myspc.m_is_particle_filter )
            {
                ParallelFor( np,
                    [=] AMREX_GPU_DEVICE (long i) {
                        Real x = pstructs[i].pos(0);
                        Real y = pstructs[i].pos(1);
                        Real z = pstructs[i].pos(2);
                        Real t = WarpX::GetInstance().gett_new(lev);
                        if ( function_partparser->getField(x,y,z,t) == 0.0 )
                        { pstructs[i].id() *= -1; }
                    }
                );
            }
            // fraction
            if ( type == "uniform" )
            {
                ParallelFor( np,
                    [=] AMREX_GPU_DEVICE (long i) {
                        if ( i % inv_fraction != 0 && pstructs[i].id() > 0 )
                        {
                            pstructs[i].id() *= -1;
                        }
                    }
                );
            }
            else if ( type == "random" )
            {
                ParallelFor( np,
                    [=] AMREX_GPU_DEVICE (long i) {
                        if ( Random() > fraction && pstructs[i].id() > 0 )
                        {
                            pstructs[i].id() *= -1;
                        }
                    }
                );
            }
            else if ( type == "default" ) {} // do nothing
            else
            { amrex::Abort("Unknown plot downsampling type."); }
        }
    }
}

// This function reset unselect particle IDs back to positive
void
PhysicalParticleContainer::UnSelectParticlesAfterIO()
{
    BL_PROFILE("PPC::UnSelectParticlesAfterIO()");

    const int nLevels = finestLevel();
    for (int lev=0; lev<=nLevels; lev++){
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto pstructs = pti.GetArrayOfStructs()().dataPtr();
            // Loop over particles and change id back to positive
            const long np = pti.numParticles();
            ParallelFor( np,
                [=] AMREX_GPU_DEVICE (long i) {
                    if ( pstructs[i].id() < 0 )
                    { pstructs[i].id() *= -1; }
                }
            );
        }
    }
}
