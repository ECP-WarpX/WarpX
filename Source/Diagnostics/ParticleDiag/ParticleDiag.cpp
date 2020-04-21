#include "WarpX.H"
#include "ParticleDiag.H"
// #include "Utils/WarpXUtil.H"

#include <AMReX_ParmParse.H>

using namespace amrex;

ParticleDiag::ParticleDiag(std::string diag_name, std::string name, WarpXParticleContainer* pc)
    : m_diag_name(diag_name), m_name(name), m_pc(pc)
{
    ParmParse pp(diag_name + "." + name);
    if (!pp.queryarr("variables", variables)){
        variables = {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "ux", "uy", "uz", "w"};
    }

    //variable to set plot_flags size
    int plot_flag_size = PIdx::nattribs;
    if(WarpX::do_back_transformed_diagnostics && m_pc->doBackTransformedDiagnostics())
        plot_flag_size += 6;

#ifdef WARPX_QED
    if(m_pc->DoQED()){
        // plot_flag will have an entry for the optical depth
        plot_flag_size++;
    }
#endif

    // Set plot_flag to 0 for all attribs
    plot_flags.resize(plot_flag_size, 0);

    // If not none, set plot_flags values to 1 for elements in plot_vars.
    if (variables[0] != "none"){
        for (const auto& var : variables){
            // Return error if var not in PIdx.
            WarpXUtilMsg::AlwaysAssert(
                ParticleStringNames::to_index.count(var),
                "ERROR: variables argument '" + var +
                "' not in ParticleStringNames"
                );
            plot_flags[ParticleStringNames::to_index.at(var)] = 1;
        }
    }

#ifdef WARPX_DIM_RZ
    // Always write out theta, whether or not it's requested,
    // to be consistent with always writing out r and z.
    plot_flags[ParticleStringNames::to_index.at("theta")] = 1;
#endif

#ifdef WARPX_QED
    if(m_pc->DoQED()){
        //Optical depths is always plotted if QED is on
        plot_flags[plot_flag_size-1] = 1;
    }
#endif

    // build filter functors
    m_do_random_filter = pp.query("random_fraction", m_random_fraction);
    m_do_uniform_filter = pp.query("uniform_stride",  m_uniform_stride);
    std::string buf;
    m_do_parser_filter = pp.query("plot_filter_function(t,x,y,z,ux,uy,uz)", buf);

    if (m_do_parser_filter) {
        std::string function_string = "";
        Store_parserString(pp,"plot_filter_function(t,x,y,z,ux,uy,uz)",
                           function_string);
        m_particle_filter_parser.reset(new ParserWrapper<7>(
            makeParser(function_string,{"t","x","y","z","ux","uy","uz"})));
    }
}
