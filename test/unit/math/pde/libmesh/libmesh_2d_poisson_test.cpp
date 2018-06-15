#ifdef STAN_LIBMESH

// solve a simple
// Poisson system.  This example also introduces the notion
// of customized matrix assembly functions, working with an
// exact solution, and using element iterators.
// We will not comment on things that
// were already explained in the second example.

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/enum_solver_package.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void assemble_poisson(EquationSystems & es,
                      const std::string & system_name);

/**
 * This is the exact solution that
 * we are trying to obtain.  We will solve
 *
 * - (u_xx + u_yy) = f
 *
 * and take a finite difference approximation using this
 * function to get f.  This is the well-known "method of
 * manufactured solutions".
 */
Real exact_solution (const Real x,
                     const Real y,
                     const Real z = 0.)
{
  static const Real pi = acos(-1.);

  return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
}


int libmesh_2d_poisson(int argc, char ** argv, double& sol2)
{
  // Initialize libraries, like in example 2.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
  // elements instead of the default QUAD4's we used in example 2
  // allow us to use higher-order approximation.
  MeshTools::Generation::build_square (mesh,
                                       15, 15,
                                       -1., 1.,
                                       -1., 1.,
                                       QUAD9);

  // Print information about the mesh to the screen.
  // Note that 5x5 QUAD9 elements actually has 11x11 nodes,
  // so this mesh is significantly larger than the one in example 2.
  // mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the Poisson system and its variables.
  // The Poisson system is another example of a steady system.
  equation_systems.add_system<LinearImplicitSystem> ("Poisson");

  // Adds the variable "u" to "Poisson".  "u"
  // will be approximated using second-order approximation.
  equation_systems.get_system("Poisson").add_variable("u", SECOND);

  // Give the system a pointer to the matrix assembly
  // function.  This will be called when needed by the
  // library.
  equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // // Prints information about the system to the screen.
  // equation_systems.print_info();

  // Solve the system "Poisson".  Note that calling this
  // member will assemble the linear system and invoke
  // the default numerical solver.  With PETSc the solver can be
  // controlled from the command line.  For example,
  // you can invoke conjugate gradient with:
  //
  // ./introduction_ex3 -ksp_type cg
  //
  // You can also get a nice X-window that monitors the solver
  // convergence with:
  //
  // ./introduction-ex3 -ksp_xmonitor
  //
  // if you linked against the appropriate X libraries when you
  // built PETSc.
  equation_systems.get_system("Poisson").solve();
  
  // build sol
  std::vector<double> u;
  equation_systems.build_solution_vector(u);
  sol2 = stan::math::dot_self(u);

#if defined(LIBMESH_HAVE_VTK) && !defined(LIBMESH_ENABLE_PARMESH)

  // After solving the system write the solution
  // to a VTK-formatted plot file.
  VTKIO (mesh).write_equation_systems ("out.pvtu", equation_systems);

#endif // #ifdef LIBMESH_HAVE_VTK

  // All done.
  return 0;
}


// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_poisson(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Poisson");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Poisson");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.  Introduction Example 4
  // describes some advantages of  std::unique_ptr's in the context of
  // quadrature rules.
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // Boundary integration requires one quadrature rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim-1, FIFTH);

  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point> & q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".  These datatypes are templated on
  //  Number, which allows the same code to work for real
  // or complex numbers.
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element iterators are a nice way to iterate through all the
  // elements, or all the elements that have some property.  The
  // iterator el will iterate from the first to the last element on
  // the local processor.  The iterator end_el tells us when to stop.
  // It is smart to make this one const so that we don't accidentally
  // mess it up!  In case users later modify this program to include
  // refinement, we will be safe and will only consider the active
  // elements; hence we use a variant of the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).

      // The  DenseMatrix::resize() and the  DenseVector::resize()
      // members will automatically zero out the matrix  and vector.
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          // Now we will build the element matrix.  This involves
          // a double loop to integrate the test functions (i) against
          // the trial functions (j).
          for (std::size_t i=0; i<phi.size(); i++)
            for (std::size_t j=0; j<phi.size(); j++)
              {
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
              }

          // This is the end of the matrix summation loop
          // Now we build the element right-hand-side contribution.
          // This involves a single loop in which we integrate the
          // "forcing function" in the PDE against the test functions.
          {
            const Real x = q_point[qp](0);
            const Real y = q_point[qp](1);
            const Real eps = 1.e-3;


            // "fxy" is the forcing function for the Poisson equation.
            // In this case we set fxy to be a finite difference
            // Laplacian approximation to the (known) exact solution.
            //
            // We will use the second-order accurate FD Laplacian
            // approximation, which in 2D is
            //
            // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
            //                u(i-1,j) + u(i+1,j) +
            //                -4*u(i,j))/h^2
            //
            // Since the value of the forcing function depends only
            // on the location of the quadrature point (q_point[qp])
            // we will compute it here, outside of the i-loop
            const Real fxy = -(exact_solution(x, y-eps) +
                               exact_solution(x, y+eps) +
                               exact_solution(x-eps, y) +
                               exact_solution(x+eps, y) -
                               4.*exact_solution(x, y))/eps/eps;

            for (std::size_t i=0; i<phi.size(); i++)
              Fe(i) += JxW[qp]*fxy*phi[i][qp];
          }
        }

      // We have now reached the end of the RHS summation,
      // and the end of quadrature point loop, so
      // the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions.
      //
      // There are several ways Dirichlet boundary conditions
      // can be imposed.  A simple approach, which works for
      // interpolary bases like the standard Lagrange polynomials,
      // is to assign function values to the
      // degrees of freedom living on the domain boundary. This
      // works well for interpolary bases, but is more difficult
      // when non-interpolary (e.g Legendre or Hierarchic) bases
      // are used.
      //
      // Dirichlet boundary conditions can also be imposed with a
      // "penalty" method.  In this case essentially the L2 projection
      // of the boundary values are added to the matrix. The
      // projection is multiplied by some large factor so that, in
      // floating point arithmetic, the existing (smaller) entries
      // in the matrix and right-hand-side are effectively ignored.
      //
      // This amounts to adding a term of the form (in latex notation)
      //
      // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i
      //
      // where
      //
      // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1
      {

        // The following loop is over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (auto side : elem->side_index_range())
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              // The value of the shape functions at the quadrature
              // points.
              const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();

              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              // The XYZ locations (in physical space) of the
              // quadrature points on the face.  This is where
              // we will interpolate the boundary value function.
              const std::vector<Point> & qface_point = fe_face->get_xyz();

              // Compute the shape function values on the element
              // face.
              fe_face->reinit(elem, side);

              // Loop over the face quadrature points for integration.
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  // The location on the boundary of the current
                  // face quadrature point.
                  const Real xf = qface_point[qp](0);
                  const Real yf = qface_point[qp](1);

                  // The penalty value.  \frac{1}{\epsilon}
                  // in the discussion above.
                  const Real penalty = 1.e10;

                  // The boundary value.
                  const Real value = exact_solution(xf, yf);

                  // Matrix contribution of the L2 projection.
                  for (std::size_t i=0; i<phi_face.size(); i++)
                    for (std::size_t j=0; j<phi_face.size(); j++)
                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

                  // Right-hand-side contribution of the L2
                  // projection.
                  for (std::size_t i=0; i<phi_face.size(); i++)
                    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                }
            }
      }

      // We have now finished the quadrature point loop,
      // and have therefore applied all the boundary conditions.

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The  SparseMatrix::add_matrix()
      // and  NumericVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

  // All done!
}

TEST(pde_solvers, libmesh_2d_poisson) {
  char arg0[] = "libmesh_setup";
  char arg1[] = "arg1";
  char arg2[] = "arg2";
  char *argv[] = {arg0, arg1, arg2, NULL};
  int argc = sizeof(argv) / sizeof(char*) - 1;
  double sol2;
  libmesh_2d_poisson(argc, argv, sol2);
  const double sol2_numeric = 240.001;
  ASSERT_FLOAT_EQ(sol2_numeric, sol2);
}

#endif
