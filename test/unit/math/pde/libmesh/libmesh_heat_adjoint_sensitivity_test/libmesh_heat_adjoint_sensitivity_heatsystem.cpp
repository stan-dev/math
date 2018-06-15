// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "libmesh/getpot.h"

#include "libmesh_heat_adjoint_sensitivity_heatsystem.h"

#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"

void HeatSystem::init_data ()
{
  T_var = this->add_variable ("T", static_cast<Order>(_fe_order),
                              Utility::string_to_enum<FEFamily>(_fe_family));

  _k = 1.0;

  parameters.push_back(_k);

  // Set equation system parameters _k and theta, so they can be read by the exact solution
  this->get_equation_systems().parameters.set<Real>("_k") = _k;

  // The temperature is evolving, with a first order time derivative
  this->time_evolving(T_var, 1);

  const boundary_id_type all_ids[4] = {0, 1, 2, 3};
  std::set<boundary_id_type> all_bdys(all_ids, all_ids+(2*2));

  std::vector<unsigned int> T_only(1, T_var);

  ZeroFunction<Number> zero;

  // Most DirichletBoundary users will want to supply a "locally
  // indexed" functor
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary (all_bdys, T_only, zero,
                        LOCAL_VARIABLE_ORDER));

  FEMSystem::init_data();
}



void HeatSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);

  // Now make sure we have requested all the data
  // we need to build the linear system.
  elem_fe->get_JxW();
  elem_fe->get_dphi();
  elem_fe->get_xyz();

  // We'll have a more automatic solution to preparing adjoint
  // solutions for time integration, eventually...
  if (c.is_adjoint())
    {
      // A reference to the system context is built with
      const System & sys = c.get_system();

      // Get a pointer to the adjoint solution vector
      NumericVector<Number> & adjoint_solution =
        const_cast<System &>(sys).get_adjoint_solution(0);

      // Add this adjoint solution to the vectors that diff context should localize
      c.add_localized_vector(adjoint_solution, sys);
    }

  FEMSystem::init_context(context);
}

#define optassert(X) {if (!(X)) libmesh_error_msg("Assertion " #X " failed.");}
//#define optassert(X) libmesh_assert(X);

bool HeatSystem::element_time_derivative (bool request_jacobian,
                                          DiffContext & context)
{
  bool compute_jacobian = request_jacobian && _analytic_jacobians;

  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // Element basis functions
  const std::vector<std::vector<RealGradient>> & dphi = elem_fe->get_dphi();

  // Workaround for weird FC6 bug
  optassert(c.get_dof_indices().size() > 0);

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(0).size();

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> & K = c.get_elem_jacobian(0, 0);
  DenseSubVector<Number> & F = c.get_elem_residual(0);

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution gradient at the Newton iterate
      Gradient grad_T = c.interior_gradient(0, qp);

      for (unsigned int i=0; i != n_u_dofs; i++)
        F(i) += JxW[qp] * -parameters[0] * (grad_T * dphi[i][qp]);
      if (compute_jacobian)
        for (unsigned int i=0; i != n_u_dofs; i++)
          for (unsigned int j=0; j != n_u_dofs; ++j)
            K(i,j) += JxW[qp] * -parameters[0] * (dphi[i][qp] * dphi[j][qp]);
    } // end of the quadrature point qp-loop

  return compute_jacobian;
}

// Perturb and accumulate dual weighted residuals
void HeatSystem::perturb_accumulate_residuals(ParameterVector & parameters_in)
{
  const unsigned int Np = parameters_in.size();

  this->update();

  for (unsigned int j=0; j != Np; ++j)
    {
      Number old_parameter = *parameters_in[j];

      *parameters_in[j] = old_parameter - dp;

      this->assembly(true, false);

      this->rhs->close();

      std::unique_ptr<NumericVector<Number>> R_minus = this->rhs->clone();

      // The contribution at a single time step would be [f(z;p+dp) - <partialu/partialt, z>(p+dp) - <g(u),z>(p+dp)] * dt
      // But since we compute the residual already scaled by dt, there is no need for the * dt
      R_minus_dp += -R_minus->dot(this->get_adjoint_solution(0));

      *parameters_in[j] = old_parameter + dp;

      this->assembly(true, false);

      this->rhs->close();

      std::unique_ptr<NumericVector<Number>> R_plus = this->rhs->clone();

      R_plus_dp += -R_plus->dot(this->get_adjoint_solution(0));

      *parameters_in[j] = old_parameter;
    }
}

// We only have one QoI, so we don't bother checking the qois argument
// to see if it was requested from us
void HeatSystem::element_qoi_derivative (DiffContext & context,
                                         const QoISet & /* qois */)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // The basis functions for the element
  const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.get_dof_indices(0).size();
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Fill the QoI RHS corresponding to this QoI. ince this is the 0th QoI
  // we fill in the [0][i] subderivatives, i corresponding to the variable index.
  // Our system has only one variable, so we only have to fill the [0][0] subderivative
  DenseSubVector<Number> & Q = c.get_qoi_derivatives(0, 0);

  // A reference to the system context is built with
  const System & sys = c.get_system();

  // Get a pointer to the adjoint solution vector
  NumericVector<Number> & adjoint_solution = const_cast<System &>(sys).get_adjoint_solution(0);

  // Get the previous adjoint solution values at all the qps

  std::vector<Number> old_adjoint (n_qpoints, 0);

  c.interior_values<Number>(0, adjoint_solution, old_adjoint);

  // Our QoI depends solely on the final time, so there are no QoI contributions.
  // However, there is a contribution from the adjoint solution timestep, for the
  // time part of the residual of the adjoint problem
  // Loop over the qps
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    for (unsigned int i=0; i != n_T_dofs; i++)
      Q(i) += -JxW[qp] * old_adjoint[qp] * phi[i][qp];
}

// Define the postprocess function to compute QoI 0
void HeatSystem::element_postprocess (DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe( 0, elem_fe );

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  const std::vector<Point> & xyz = elem_fe->get_xyz();

  // The number of local degrees of freedom in each variable

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  Number dQoI_0 = 0.;

  // Loop over quadrature points

  for (unsigned int qp = 0; qp != n_qpoints; qp++)
    {
      // Get co-ordinate locations of the current quadrature point
      const Real x = xyz[qp](0);
      const Real y = xyz[qp](1);

      Number T = c.interior_value(0, qp);

      // Flux with weight s = R(u^h, s) = int ( f*s - alpha*(grad_u*grad_s) ) dx
      dQoI_0 += JxW[qp] * ( T * sin(M_PI * x) * sin(M_PI * y) );
    }

  // Update the computed value of the global functional R, by adding the contribution from this element

  computed_QoI[0] = computed_QoI[0] + dQoI_0;

}
