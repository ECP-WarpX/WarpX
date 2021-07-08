/* Copyright 2019 Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXParser.H"

#include <algorithm>

WarpXParser::WarpXParser (std::string const& func_body)
{
    define(func_body);
}

void
WarpXParser::define (std::string const& func_body)
{
    clear();

    m_expression = func_body;
    m_expression.erase(std::remove(m_expression.begin(),m_expression.end(),'\n'),
                       m_expression.end());
    std::string f = m_expression + "\n";

#ifdef AMREX_USE_OMP

    int nthreads = omp_get_max_threads();
    m_variables.resize(nthreads);
    m_varnames.resize(nthreads);
    m_parser.resize(nthreads);
    m_parser[0] = wp_c_parser_new(f.c_str());
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        if (tid > 0) {
            m_parser[tid] = wp_parser_dup(m_parser[0]);
        }
    }

#else

    m_parser = wp_c_parser_new(f.c_str());

#endif
}

WarpXParser::~WarpXParser ()
{
    clear();
}

void
WarpXParser::clear ()
{
    m_expression.clear();
    m_varnames.clear();

#ifdef AMREX_USE_OMP

    if (!m_parser.empty())
    {
#pragma omp parallel
        {
            wp_parser_delete(m_parser[omp_get_thread_num()]);
        }
    }
    m_parser.clear();
    m_variables.clear();

#else

    if (m_parser) wp_parser_delete(m_parser);
    m_parser = nullptr;

#endif
}

void
WarpXParser::registerVariable (std::string const& name, amrex::Real& var)
{
    // We assume this is called inside OMP parallel region
#ifdef AMREX_USE_OMP
    wp_parser_regvar(m_parser[omp_get_thread_num()], name.c_str(), &var);
    m_varnames[omp_get_thread_num()].push_back(name);
#else
    wp_parser_regvar(m_parser, name.c_str(), &var);
    m_varnames.push_back(name);
#endif
}

void
WarpXParser::registerVariables (std::vector<std::string> const& names)
{
#ifdef AMREX_USE_OMP

// This must be called outside OpenMP parallel region.
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        struct wp_parser* p = m_parser[tid];
        auto& v = m_variables[tid];
        for (int j = 0; j < static_cast<int>(names.size()); ++j) {
            wp_parser_regvar(p, names[j].c_str(), &(v[j]));
            m_varnames[tid].push_back(names[j]);
        }
    }

#else

    for (auto j = 0u; j < names.size(); ++j) {
        wp_parser_regvar(m_parser, names[j].c_str(), &(m_variables[j]));
        m_varnames.push_back(names[j]);
    }

#endif
}

void
WarpXParser::setConstant (std::string const& name, amrex::Real c)
{
#ifdef AMREX_USE_OMP

    bool in_parallel = omp_in_parallel();
    // We don't know if this is inside OMP parallel region or not
#pragma omp parallel if (!in_parallel)
    {
        wp_parser_setconst(m_parser[omp_get_thread_num()], name.c_str(), c);
    }

#else

    wp_parser_setconst(m_parser, name.c_str(), c);

#endif
}

void
WarpXParser::print () const
{
#ifdef AMREX_USE_OMP
#pragma omp critical(warpx_parser_pint)
    wp_ast_print(m_parser[omp_get_thread_num()]->ast);
#else
    wp_ast_print(m_parser->ast);
#endif
}

int
WarpXParser::depth () const
{
    int n = 0;
#ifdef AMREX_USE_OMP
    wp_ast_depth(m_parser[omp_get_thread_num()]->ast, &n);
#else
    wp_ast_depth(m_parser->ast, &n);
#endif
    return n;
}

std::string const&
WarpXParser::expr () const
{
    return m_expression;
}

std::set<std::string>
WarpXParser::symbols () const
{
    std::set<std::string> results;
#ifdef AMREX_USE_OMP
    wp_ast_get_symbols(m_parser[omp_get_thread_num()]->ast, results);
#else
    wp_ast_get_symbols(m_parser->ast, results);
#endif
    return results;
}
