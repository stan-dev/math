#ifdef STAN_LIBMESH

// <h1>Based on libmesh's Introduction Example 2 - Defining a Simple System</h1>
// by Benjamin S. Kirk
//
// This is the second example program.  It demonstrates how to
// create an equation system for a simple scalar system.  This
// example will also introduce some of the issues involved with using PETSc
// in your application.
//
// The test also indirectly
// uses the PETSc library.  By default equation data is stored
// in PETSc vectors, which may span multiple processors.  Before
// PETSc is used it must be initialized via libMesh::init().  Note that
// by passing argc and argv to PETSc you may specify
// command line arguments to PETSc.  For example, you might
// try running this example as:
//
// ./introduction_ex2 -log_info
//
// to see what PETSc is doing behind the scenes or
//
// ./introduction_ex2 -log_summary
//
// to get a summary of what PETSc did.
// Among other things, libMesh::init() initializes the MPI
// communications library and PETSc numeric library on your system if
// you haven't already done so.

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

// C++ include files that we need
#include <iostream>
//Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/enum_xdr_mode.h"
// Include file that defines various mesh generation utilities
#include "libmesh/mesh_generation.h"
// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/enum_solver_package.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;



int libmesh_setup_test(int argc, char ** argv)
{
  LibMeshInit init (argc, argv);

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // // A brief message to the user to inform her of the
  // // exact name of the program being run, and its command line.
  // libMesh::out << "Running " << argv[0];
  // for (int i=1; i<argc; i++)
  //   libMesh::out << " " << argv[i];
  // libMesh::out << std::endl << std::endl;

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the unit square.  By default a mesh of QUAD4
  // elements will be created.  We instruct the mesh generator
  // to build a mesh of 5x5 elements.
  MeshTools::Generation::build_square (mesh, 5, 5);

  // Create an equation systems object. This object can
  // contain multiple systems of different
  // flavors for solving loosely coupled physics.  Each system can
  // contain multiple variables of different approximation orders.
  // Here we will simply create a single system with one variable.
  // Later on, other flavors of systems will be introduced.  For the
  // moment, we use the general system.
  // The EquationSystems object needs a reference to the mesh
  // object, so the order of construction here is important.
  EquationSystems equation_systems (mesh);

  // Add a flag "test" that is visible for all systems.  This
  // helps in inter-system communication.
  equation_systems.parameters.set<bool> ("test") = true;

  // Set a simulation-specific parameter visible for all systems.
  // This helps in inter-system-communication.
  equation_systems.parameters.set<Real> ("dummy") = 42.;

  // Set another simulation-specific parameter
  equation_systems.parameters.set<Real> ("nobody") = 0.;

  // Now we declare the system and its variables.
  // We begin by adding a "TransientLinearImplicitSystem" to the
  // EquationSystems object, and we give it the name
  // "Simple System".
  equation_systems.add_system<TransientLinearImplicitSystem> ("Simple System");

  // Adds the variable "u" to "Simple System".  "u"
  // will be approximated using first-order approximation.
  equation_systems.get_system("Simple System").add_variable("u", FIRST);

  // add an "ExplicitSystem" to the
  // EquationSystems object, and we give it the name
  // "Complex System".
  equation_systems.add_system<ExplicitSystem> ("Complex System");

  // Give "Complex System" three variables -- each with a different approximation
  // order.  Variables "c" and "T" will use first-order Lagrange approximation,
  // while variable "dv" will use a second-order discontinuous
  // approximation space.
  equation_systems.get_system("Complex System").add_variable("c", FIRST);
  equation_systems.get_system("Complex System").add_variable("T", FIRST);
  equation_systems.get_system("Complex System").add_variable("dv", SECOND, MONOMIAL);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // unit test is based on the ndofs of the system
  return equation_systems.n_dofs();
}

TEST(pde_solvers, libmesh_setup) {
  char arg0[] = "libmesh_setup";
  char arg1[] = "arg1";
  char arg2[] = "arg2";
  char *argv[] = {arg0, arg1, arg2, NULL};
  int argc = sizeof(argv) / sizeof(char*) - 1;
  ASSERT_EQ(258, libmesh_setup_test(argc, argv));
}

#endif
