/* Copyright 2026 The WarpX Community
 *
 * Authors: Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "Python/pyWarpX.H"

#include <FieldSolver/ImplicitSolvers/ImplicitSolver.H>
#include <WarpX.H>

#include <ablastr/fields/MultiFabRegister.H>

#include <memory>
#include <string>


void init_ImplicitSolver (py::module& m)
{
    // The implicit solver is owned by WarpX (WarpX::m_implicit_solver) and is
    // destroyed with it in WarpX::Finalize, so the Python object must never
    // delete it: py::nodelete makes the returned pointer non-owning.
    py::class_<ImplicitSolver, std::unique_ptr<ImplicitSolver, py::nodelete>>(m, "ImplicitSolver")
        .def("finish_mass_matrices", &ImplicitSolver::FinishMassMatrices,
            R"pbdoc(Fill the second half of the symmetric diagonal mass matrices

The deposition only fills half of the band of the diagonal blocks (``Sxx``,
``Syy``, ``Szz``), exploiting their symmetry. This mirrors the other half.)pbdoc"
        )
        .def("apply_mass_matrices",
            [](ImplicitSolver& solver,
               std::string const& out_name,
               std::string const& in_name,
               bool zero_out_first)
            {
                auto& fields = WarpX::GetInstance().GetMultiFabRegister();
                // the mass matrices only exist on the levels the solver was set up for
                const int finest_level = solver.numAMRLevels() - 1;

                ablastr::fields::MultiLevelVectorField out =
                    fields.get_mr_levels_alldirs(out_name, finest_level);
                const ablastr::fields::MultiLevelVectorField in =
                    fields.get_mr_levels_alldirs(in_name, finest_level);

                solver.ApplyMassMatrices(out, in, nullptr, nullptr, 1.0, zero_out_first);
            },
            py::arg("out_name"), py::arg("in_name"), py::arg("zero_out_first") = false,
            R"pbdoc(Apply the mass matrices to a vector field

Computes ``out += S * in``, where ``S`` are the mass matrices, the linear
response of the deposited current density to the electric field
(``dJ = S dE``). Both fields are looked up by name in the MultiFab register.
``in`` must have its guard cells filled: the stencil reads the neighbors of
every point it writes.

Parameters
----------
out_name: str
  Name of the vector field to accumulate into, with the staggering of the current density
in_name: str
  Name of the vector field to apply the mass matrices to, with the staggering of the electric field
zero_out_first: bool, optional
  Whether to zero ``out`` first, so that it holds the result of this operation only;
  otherwise the result is added to its existing contents (defaults to False))pbdoc"
        )
    ;
}
