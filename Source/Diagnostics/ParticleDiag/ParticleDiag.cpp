#include "ParticleDiag.H"

#include "Diagnostics/ParticleDiag/ParticleDiag.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <AMReX_ParmParse.H>

#include <map>
#include <vector>

using namespace amrex;

ParticleDiag::ParticleDiag(std::string diag_name, std::string name, WarpXParticleContainer* pc, PinnedMemoryParticleContainer* pinned_pc)
    : m_diag_name(diag_name), m_name(name), m_pc(pc), m_pinned_pc(pinned_pc)
{
    //variable to set m_plot_flags size
    const int plot_flag_size = pc->NumRealComps();

    // By default output all attributes
    m_plot_flags.resize(plot_flag_size, 1);

    const ParmParse pp_diag_name_species_name(diag_name + "." + name);
    amrex::Vector<std::string> variables;
    const int variables_specified = pp_diag_name_species_name.queryarr("variables", variables);

    if (variables_specified){
        // If only specific variables have been specified, fill m_plot_flags with zero and only set
        // requested variables to one
        std::fill(m_plot_flags.begin(), m_plot_flags.end(), 0);
        if (variables[0] != "none"){
            const std::map<std::string, int> existing_variable_names = pc->getParticleComps();
            for (const auto& var : variables){
                const auto search = existing_variable_names.find(var);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    search != existing_variable_names.end(),
                    "variables argument '" + var
                    +"' is not an existing attribute for this species");
                m_plot_flags[existing_variable_names.at(var)] = 1;
            }
        }
    }

#ifdef WARPX_DIM_RZ
    // Always write out theta, whether or not it's requested,
    // to be consistent with always writing out r and z.
    // TODO: openPMD does a reconstruction to Cartesian, so we can now skip force-writing this
    m_plot_flags[pc->getParticleComps().at("theta")] = 1;
#endif

    // build filter functors
    m_do_random_filter = utils::parser::queryWithParser(
        pp_diag_name_species_name, "random_fraction", m_random_fraction);
    m_do_uniform_filter = utils::parser::queryWithParser(
        pp_diag_name_species_name, "uniform_stride",m_uniform_stride);
    std::string buf;
    m_do_parser_filter = pp_diag_name_species_name.query("plot_filter_function(t,x,y,z,ux,uy,uz)",
                                                         buf);

    if (m_do_parser_filter) {
        std::string function_string = "";
        utils::parser::Store_parserString(
            pp_diag_name_species_name,"plot_filter_function(t,x,y,z,ux,uy,uz)", function_string);
        m_particle_filter_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(function_string,{"t","x","y","z","ux","uy","uz"}));
    }
}
