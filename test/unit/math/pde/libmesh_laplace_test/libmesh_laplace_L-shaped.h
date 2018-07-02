#ifndef STAN_MATH_PDE_LAPALCE_TEST_L_SHAPED_H
#define STAN_MATH_PDE_LAPALCE_TEST_L_SHAPED_H

#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_system.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/qoi_set.h"
#include "libmesh/system.h"
#include "libmesh/parameter_pointer.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

#include "libmesh/getpot.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parallel.h"
#include "libmesh/fem_context.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class LaplaceSystem : public FEMSystem
{
public:
  // Constructor
  LaplaceSystem(EquationSystems & es,
                const std::string & name_in,
                const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
      _fe_family("LAGRANGE"),
      _fe_order(1),
      _analytic_jacobians(true)
  { qoi.resize(1); }

  Number computed_QoI[1];

  // Parameters associated with the system
  std::vector<Number> parameters;

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  // Postprocessing function which we are going to override for this application

  virtual void postprocess();

  Number & get_parameter_value(unsigned int parameter_index)
  {
    return parameters[parameter_index];
  }

  ParameterVector & get_parameter_vector()
  {
    parameter_vector.clear();
    for (std::size_t i = 0; i != parameters.size(); ++i)
      parameter_vector.push_back(libmesh_make_unique<ParameterPointer<Number>>(&parameters[i]));

    return parameter_vector;
  }

protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context);

  // Constraint parts
  virtual bool side_constraint (bool request_jacobian,
                                DiffContext & context);

  Number exact_solution (const Point &);

  // The ParameterVector object that will contain pointers to
  // the system parameters
  ParameterVector parameter_vector;

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;

  virtual void element_postprocess(DiffContext & context);
};

// Bring in everything from the libMesh namespace
using namespace libMesh;

void LaplaceSystem::init_data ()
{
  unsigned int T_var =
    this->add_variable ("T", static_cast<Order>(_fe_order),
                        Utility::string_to_enum<FEFamily>(_fe_family));

  GetPot infile("libmesh_laplace_l-shaped.in");
  parameters.push_back(infile("alpha_1", 1.0));
  parameters.push_back(infile("alpha_2", 1.0));

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // The temperature is evolving, with a first order time derivative
  this->time_evolving(T_var, 1);
}

void LaplaceSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // Now make sure we have requested all the data
  // we need to build the linear system.
  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);
  elem_fe->get_JxW();
  elem_fe->get_phi();
  elem_fe->get_dphi();
  elem_fe->get_xyz();

  FEBase * side_fe = libmesh_nullptr;
  c.get_side_fe(0, side_fe);
  side_fe->get_JxW();
  side_fe->get_phi();
  side_fe->get_xyz();
}

#define optassert(X) {if (!(X)) libmesh_error_msg("Assertion " #X " failed.");}

// Assemble the element contributions to the stiffness matrix
bool LaplaceSystem::element_time_derivative (bool request_jacobian,
                                             DiffContext & context)
{
  // Are the jacobians specified analytically ?
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

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.get_dof_indices(0).size();

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

      // The residual contribution from this element
      for (unsigned int i=0; i != n_T_dofs; i++)
        // F(i) += JxW[qp] * (grad_T * dphi[i][qp]);
        F(i) += JxW[qp] * (parameters[0] + (2.*parameters[1])) * (grad_T * dphi[i][qp]);
      if (compute_jacobian)
        for (unsigned int i=0; i != n_T_dofs; i++)
          for (unsigned int j=0; j != n_T_dofs; ++j)
            // The analytic jacobian
            // K(i,j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
            K(i,j) += JxW[qp] * (parameters[0] + (2.*parameters[1])) * (dphi[i][qp] * dphi[j][qp]);
    } // end of the quadrature point qp-loop

  return compute_jacobian;
}

// Set Dirichlet bcs, side contributions to global stiffness matrix
bool LaplaceSystem::side_constraint (bool request_jacobian,
                                     DiffContext & context)
{
  // Are the jacobians specified analytically ?
  bool compute_jacobian = request_jacobian && _analytic_jacobians;

  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * side_fe = libmesh_nullptr;
  c.get_side_fe(0, side_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = side_fe->get_JxW();

  // Side basis functions
  const std::vector<std::vector<Real>> & phi = side_fe->get_phi();

  // Side Quadrature points
  const std::vector<Point > & qside_point = side_fe->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.get_dof_indices(0).size();

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> & K = c.get_elem_jacobian(0, 0);
  DenseSubVector<Number> & F = c.get_elem_residual(0);

  unsigned int n_qpoints = c.get_side_qrule().n_points();

  const Real penalty = 1./(TOLERANCE*TOLERANCE);

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution at the old Newton iterate
      Number T = c.side_value(0, qp);

      // We get the Dirichlet bcs from the exact solution
      Number u_dirichlet = exact_solution (qside_point[qp]);

      // The residual from the boundary terms, penalize non-zero temperature
      for (unsigned int i=0; i != n_T_dofs; i++)
        F(i) += JxW[qp] * penalty * (T - u_dirichlet) * phi[i][qp];
      if (compute_jacobian)
        for (unsigned int i=0; i != n_T_dofs; i++)
          for (unsigned int j=0; j != n_T_dofs; ++j)
            // The analytic jacobian
            K(i,j) += JxW[qp] * penalty * phi[i][qp] * phi[j][qp];

    } // end of the quadrature point qp-loop

  return compute_jacobian;
}

// The exact solution to the singular problem,
// u_exact = r^(2/3)*sin(2*theta/3). We use this to set the Dirichlet boundary conditions
Number LaplaceSystem::exact_solution(const Point & p)
{
  const Real x1 = p(0);
  const Real x2 = p(1);

  Real theta = atan2(x2, x1);

  // Make sure 0 <= theta <= 2*pi
  if (theta < 0)
    theta += 2. * libMesh::pi;

  return (parameters[0] + (2.*parameters[1])) * pow(x1*x1 + x2*x2, 1./3.)*sin(2./3.*theta);
  // return pow(x1*x1 + x2*x2, 1./3.)*sin(2./3.*theta);
}

void LaplaceSystem::postprocess()
{
  // Reset the array holding the computed QoIs
  computed_QoI[0] = 0.0;

  FEMSystem::postprocess();

  this->comm().sum(computed_QoI[0]);
}

void LaplaceSystem::element_postprocess (DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();
  const std::vector<Point> & xyz = elem_fe->get_xyz();

  // The number of local degrees of freedom in each variable
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // The function R = int_{omega} T dR
  // omega is a subset of Omega (the whole domain), omega = [0.75, 1.0] x [0.0, 0.25]
  Number dQoI_0 = 0.;

  // Loop over quadrature points
  for (unsigned int qp = 0; qp != n_qpoints; qp++)
    {
      // Get co-ordinate locations of the current quadrature point
      const Real x = xyz[qp](0);
      const Real y = xyz[qp](1);

      // If in the sub-domain omega, add the contribution to the integral R
      if (std::abs(x - 0.875) <= 0.125 && std::abs(y - 0.125) <= 0.125)
        {
          // Get the solution value at the quadrature point
          Number T = c.interior_value(0, qp);

          // Update the elemental increment dR for each qp
          dQoI_0 += JxW[qp] * T;
        }
    }

  // Update the computed value of the global functional R, by adding the contribution from this element
  computed_QoI[0] = computed_QoI[0] + dQoI_0;
}

#endif
