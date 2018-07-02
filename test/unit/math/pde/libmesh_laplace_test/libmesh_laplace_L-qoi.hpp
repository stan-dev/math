#ifndef STAN_MATH_PDE_LAPALCE_TEST_L_QOI_H
#define STAN_MATH_PDE_LAPALCE_TEST_L_QOI_H

#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/diff_qoi.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// Bring in everything from the libMesh namespace
using namespace libMesh;

class LaplaceQoI : public DifferentiableQoI
{
public:
  LaplaceQoI(){}
  virtual ~LaplaceQoI(){}

  virtual void init_qoi(std::vector<Number> & sys_qoi);

  // Context initialization
  virtual void init_context (DiffContext & context);

  virtual void postprocess() {}

  virtual void element_qoi_derivative(DiffContext & context, const QoISet & qois);

  virtual void element_qoi (DiffContext & context, const QoISet & qois);

  virtual std::unique_ptr<DifferentiableQoI> clone()
  {
    return libmesh_make_unique<LaplaceQoI>(*this);
  }

};

using namespace libMesh;

void LaplaceQoI::init_qoi(std::vector<Number> & sys_qoi)
{
  // Only 1 qoi to worry about
  sys_qoi.resize(1);
}

void LaplaceQoI::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // Now make sure we have requested all the data
  // we need to build the linear system.
  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);
  elem_fe->get_JxW();
  elem_fe->get_phi();
  elem_fe->get_xyz();
}


// We only have one QoI, so we don't bother checking the qois argument
// to see if it was requested from us
void LaplaceQoI::element_qoi (DiffContext & context,
                              const QoISet & /* qois */)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  const std::vector<Point> & xyz = elem_fe->get_xyz();

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  Number dQoI_0 = 0.;

  // Loop over quadrature points
  for (unsigned int qp = 0; qp != n_qpoints; qp++)
    {
      // Get co-ordinate locations of the current quadrature point
      const Real xf = xyz[qp](0);
      const Real yf = xyz[qp](1);

      // If in the sub-domain omega, add the contribution to the integral R
      if (std::abs(xf - 0.875) <= 0.125 && std::abs(yf - 0.125) <= 0.125)
        {
          // Get the solution value at the quadrature point
          Number T = c.interior_value(0, qp);

          // Update the elemental increment dR for each qp
          dQoI_0 += JxW[qp] * T;
        }
    }

  // Update the computed value of the global functional R, by adding the contribution from this element
  c.get_qois()[0] = c.get_qois()[0] + dQoI_0;
}

// We only have one QoI, so we don't bother checking the qois argument
// to see if it was requested from us
void LaplaceQoI::element_qoi_derivative (DiffContext & context,
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

  // The element quadrature points
  const std::vector<Point > & q_point = elem_fe->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.get_dof_indices(0).size();
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
  // we fill in the [0][i] subderivatives, i corresponding to the variable index.
  // Our system has only one variable, so we only have to fill the [0][0] subderivative
  DenseSubVector<Number> & Q = c.get_qoi_derivatives(0, 0);

  // Loop over the qps
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);

      // If in the sub-domain over which QoI 0 is supported, add contributions
      // to the adjoint rhs
      if (std::abs(x - 0.875) <= 0.125 && std::abs(y - 0.125) <= 0.125)
        for (unsigned int i=0; i != n_T_dofs; i++)
          Q(i) += JxW[qp]*phi[i][qp];
    } // end of the quadrature point qp-loop
}

#endif // L_QOI_H
