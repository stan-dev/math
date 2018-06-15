#ifdef STAN_LIBMESH

// General Localized Vectors and Unsteady Adjoints.
// \author Vikram Garg
// \date 2013
//
// This example showcases three new capabilities in libMesh. The
// primary motivation is adjoint sensitivity analysis for unsteady
// problems. The PDE we are interested in is the simple 2-d heat
// equation:
// partial(T)/partial(t) - K Laplacian(T) = 0
// with initial condition:
// T(x,y;0) = sin(pi*x) sin(pi*y) and boundary conditions:
// T(boundary;t) = 0

// For these initial and boundary conditions, the exact solution
// u = exp(-K pi^2 t) * sin(pi*x) * sin(pi*y)

// We specify our Quantity of Interest (QoI) as
// Q(u) = int_{domain} u(x,y;1) sin(pi*x) sin(pi*y) dx dy, and
// are interested in computing the sensitivity dQ/dK

// For this QoI, the continuous adjoint problem reads,
// -partial(z)/partial(t) - K Laplacian(z) = 0
// with initial condition:
// T(x,y;1) = sin(pi*x) sin(pi*y)
// and boundary condition:
// T(boundary;t) = 0

// which has the exact solution,
// z = exp(-K pi^2 (1 - t)) * sin(pi*x) * sin(pi*y)
// which is the mirror image in time of the forward solution

// For an adjoint consistent space-time formulation, the discrete
// adjoint can be obtained by marching backwards from the adjoint
// initial condition and solving the transpose of the discrete primal
// problem at the last nonlinear solve of the corresponding primal
// timestep. This necessitates the storage of the primal solution at
// all timesteps, which is accomplished here using a
// MemorySolutionHistory object. As the name suggests, this object
// simply stores the primal solution (and other vectors we may choose
// to save), so that we can retrieve them later, whenever necessary.

// The discrete adjoint system for implicit time steppers requires the
// localization of vectors other than system.solution, which is
// accomplished using the localize_vectors method. In this particular
// example, we use the localized adjoint solution to assemble the
// residual contribution for the current adjoint timestep from the last
// computed adjoint timestep.

// Finally, The adjoint_advance_timestep method, the backwards time
// analog of advance_timestep prepares the time solver for solving the
// adjoint system, while the retrieve_timestep method retrieves the
// saved solutions at the current system.time, so that the adjoint
// sensitivity contribution for the current time can be computed.

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

#include "libmesh_heat_adjoint_sensitivity_femparameters.h"
#include "libmesh_heat_adjoint_sensitivity_heatsystem.h"

// Libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/system_norm.h"
#include "libmesh/numeric_vector.h"
// #include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/enum_solver_package.h"

#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/petsc_diff_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/twostep_time_solver.h"

#include "libmesh/getpot.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exodusII_io.h"

// SolutionHistory Includes
#include "libmesh/solution_history.h"
#include "libmesh/memory_solution_history.h"

// function factory
#include "libmesh/dense_vector.h"
#include "libmesh/factory.h"
#include "libmesh/function_base.h"
// #include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ includes
#include <iostream>
#include <sys/time.h>
#include <iomanip>

// adjoint initial
void adjoint_read_initial_parameters() {}

void adjoint_finish_initialization() {}


// Initial conditions
Number adjoint_initial_value(const Point & p,
                             const Parameters &,
                             const std::string &,
                             const std::string &)
{
  Real x = p(0), y = p(1);

  return (sin(M_PI * x) * sin(M_PI * y));
}


Gradient adjoint_initial_grad(const Point & p,
                              const Parameters &,
                              const std::string &,
                              const std::string &)
{
  Real x = p(0), y = p(1);

  return Gradient(Number(M_PI*cos(M_PI * x) * sin(M_PI * y)),
                  Number(M_PI*sin(M_PI * x) * cos(M_PI * y)));
}


void read_initial_parameters() {}

void finish_initialization() {}

// Initial conditions
Number initial_value(const Point & p,
                     const Parameters &,
                     const std::string &,
                     const std::string &)
{
  Real x = p(0), y = p(1);

  return sin(M_PI * x) * sin(M_PI * y);
}

Gradient initial_grad(const Point & p,
                      const Parameters &,
                      const std::string &,
                      const std::string &)
{
  Real x = p(0), y = p(1);

  return Gradient(Number(M_PI*cos(M_PI * x) * sin(M_PI * y)),
                  Number(M_PI*sin(M_PI * x) * cos(M_PI * y)));
}


// An example of a hand-coded function
class ExampleOneFunction : public FunctionBase<Number>
{
  virtual Number operator() (const Point & /*p*/,
                             const Real /*time*/)
  {
    return 1;
  }

  virtual void operator() (const Point & /*p*/,
                           const Real /*time*/,
                           DenseVector<Number> & output)
  {
    for (unsigned int i=0; i != output.size(); ++i)
      output(i) = 1;
  }

  virtual void init() {}
  virtual void clear() {}
  virtual std::unique_ptr<FunctionBase<Number>> clone() const
  {
    return libmesh_make_unique<ExampleOneFunction>();
  }
};

#ifdef LIBMESH_USE_SEPARATE_NAMESPACE
namespace libMesh {
#endif

//-------------------------------------------------
// Full specialization for the Factory<FunctionBase<Number>>
// So we can look up hand-coded functions by name string
template<>
std::map<std::string, Factory<FunctionBase<Number>> *> &
Factory<FunctionBase<Number>>::factory_map()
{
  static std::map<std::string, Factory<FunctionBase<Number>> *> _map;
  return _map;
}

FactoryImp<ExampleOneFunction, FunctionBase<Number>> example_one_factory ("example_one");

#ifdef LIBMESH_USE_SEPARATE_NAMESPACE
} // namespace libMesh
#endif


void write_output(EquationSystems & es,
                  unsigned int t_step,       // The current time step count
                  std::string solution_type, // primal or adjoint solve
                  FEMParameters & param)
{
  // Ignore parameters when there are no output formats available.
  libmesh_ignore(es);
  libmesh_ignore(t_step);
  libmesh_ignore(solution_type);
  libmesh_ignore(param);

#ifdef LIBMESH_HAVE_GMV
  if (param.output_gmv)
    {
      MeshBase & mesh = es.get_mesh();

      std::ostringstream file_name_gmv;
      file_name_gmv << solution_type
                    << ".out.gmv."
                    << std::setw(2)
                    << std::setfill('0')
                    << std::right
                    << t_step;

      GMVIO(mesh).write_equation_systems(file_name_gmv.str(), es);
    }
#endif

#ifdef LIBMESH_HAVE_EXODUS_API
  if (param.output_exodus)
    {
      MeshBase & mesh = es.get_mesh();

      // We write out one file per timestep. The files are named in
      // the following way:
      // foo.e
      // foo.e-s002
      // foo.e-s003
      // ...
      // so that, if you open the first one with Paraview, it actually
      // opens the entire sequence of adapted files.
      std::ostringstream file_name_exodus;

      file_name_exodus << solution_type << ".e";
      if (t_step > 0)
        file_name_exodus << "-s"
                         << std::setw(3)
                         << std::setfill('0')
                         << std::right
                         << t_step + 1;

      // TODO: Get the current time from the System...
      ExodusII_IO(mesh).write_timestep(file_name_exodus.str(),
                                       es,
                                       1,
                                       /*time=*/t_step + 1);
    }
#endif
}

void set_system_parameters(HeatSystem & system,
                           FEMParameters & param)
{
  // Use the prescribed FE type
  system.fe_family() = param.fe_family[0];
  system.fe_order() = param.fe_order[0];

  // Use analytical jacobians?
  system.analytic_jacobians() = param.analytic_jacobians;

  // Verify analytic jacobians against numerical ones?
  system.verify_analytic_jacobians = param.verify_analytic_jacobians;
  system.numerical_jacobian_h = param.numerical_jacobian_h;

  // More desperate debugging options
  system.print_solution_norms    = param.print_solution_norms;
  system.print_solutions         = param.print_solutions;
  system.print_residual_norms    = param.print_residual_norms;
  system.print_residuals         = param.print_residuals;
  system.print_jacobian_norms    = param.print_jacobian_norms;
  system.print_jacobians         = param.print_jacobians;
  system.print_element_jacobians = param.print_element_jacobians;
  system.print_element_residuals = param.print_element_residuals;

  // Solve this as a time-dependent or steady system
  if (param.transient)
    {
      UnsteadySolver *innersolver;
      if (param.timesolver_core == "euler")
        {
          EulerSolver *eulersolver =
            new EulerSolver(system);

          eulersolver->theta = param.timesolver_theta;
          innersolver = eulersolver;
        }
      else
        libmesh_error_msg("This example (and unsteady adjoints in libMesh) only support Backward Euler and explicit methods.");

      system.time_solver =
        std::unique_ptr<TimeSolver>(innersolver);
    }
  else
    system.time_solver = libmesh_make_unique<SteadySolver>(system);

  // The Memory Solution History object we will set the system SolutionHistory object to
  MemorySolutionHistory heatsystem_solution_history(system);
  system.time_solver->set_solution_history(heatsystem_solution_history);

  system.time_solver->reduce_deltat_on_diffsolver_failure =
    param.deltat_reductions;
  system.time_solver->quiet           = param.time_solver_quiet;

  // Create any Dirichlet boundary conditions
  typedef
    std::map<boundary_id_type, FunctionBase<Number> *>::
    const_iterator Iter;

  for (Iter i = param.dirichlet_conditions.begin();
       i != param.dirichlet_conditions.end(); ++i)
    {
      boundary_id_type b = i->first;
      FunctionBase<Number> *f = i->second;
      std::set<boundary_id_type> bdys; bdys.insert(b);

      system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(bdys,
                                                                    param.dirichlet_condition_variables[b],
                                                                    f));

      libMesh::out << "Added Dirichlet boundary " << b << " for variables ";
      for (std::size_t vi=0; vi != param.dirichlet_condition_variables[b].size(); ++vi)
        libMesh::out << param.dirichlet_condition_variables[b][vi];
      libMesh::out << std::endl;
    }

  // Set the time stepping options
  system.deltat = param.deltat;

  // And the integration options
  system.extra_quadrature_order = param.extra_quadrature_order;

  // And the nonlinear solver options
  if (param.use_petsc_snes)
    {
#ifdef LIBMESH_HAVE_PETSC
      PetscDiffSolver *solver = new PetscDiffSolver(system);
      system.time_solver->diff_solver() = std::unique_ptr<DiffSolver>(solver);
#else
      libmesh_error_msg("This example requires libMesh to be compiled with PETSc support.");
#endif
    }
  else
    {
      NewtonSolver *solver = new NewtonSolver(system);
      system.time_solver->diff_solver() = std::unique_ptr<DiffSolver>(solver);

      solver->quiet                       = param.solver_quiet;
      solver->verbose                     = param.solver_verbose;
      solver->max_nonlinear_iterations    = param.max_nonlinear_iterations;
      solver->minsteplength               = param.min_step_length;
      solver->relative_step_tolerance     = param.relative_step_tolerance;
      solver->relative_residual_tolerance = param.relative_residual_tolerance;
      solver->require_residual_reduction  = param.require_residual_reduction;
      solver->linear_tolerance_multiplier = param.linear_tolerance_multiplier;
      if (system.time_solver->reduce_deltat_on_diffsolver_failure)
        {
          solver->continue_after_max_iterations = true;
          solver->continue_after_backtrack_failure = true;
        }

      // And the linear solver options
      solver->max_linear_iterations       = param.max_linear_iterations;
      solver->initial_linear_tolerance    = param.initial_linear_tolerance;
      solver->minimum_linear_tolerance    = param.minimum_linear_tolerance;
    }
}

// stan-libmesh API. Users provide this function to
// calculate primal & sentisivity. The output is
// a vector of vectors:
// [QoI0, dQoI0/dtheta0, dQoI0/dtheta1, dQoI0/dtheta2...]
// [QoI1, dQoI1/dtheta0, dQoI1/dtheta1, dQoI1/dtheta2...]
// ...
//
// methods provided by this API struct:
// init(): envionment setup
// finalize(): envionment takedown
// solve(...): solve PDE and return QoI
// solve_with_sensitivity(...): solve PDE and return QoI &
// QoI sensitivity
class LibMeshHeatModel {
  LibMeshInit init;

public:
  LibMeshHeatModel(int argc, const char * const * argv) : init(argc, argv) {}

  inline std::vector<std::vector<double> >
  solve_with_sensitivity(const std::vector<double>& theta,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i,
                         std::ostream* msgs = nullptr) const {
  // Make sure the general input file exists, and parse it
  {
    std::ifstream i("libmesh_heat_adjoint_sensitivity_general.in");
    if (!i)
      libmesh_error_msg('[' << init.comm().rank() << "] Can't find libmesh_heat_adjoint_sensitivity_general.in; exiting early.");
  }
  GetPot infile("libmesh_heat_adjoint_sensitivity_general.in");

  // Read in parameters from the input file
  FEMParameters param(init.comm());
  param.read(infile);

  // Create a mesh with the given dimension, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm(), param.dimension);

  // And an object to refine it
  auto mesh_refinement = libmesh_make_unique<MeshRefinement>(mesh);

  // And an EquationSystems to run on it
  EquationSystems equation_systems (mesh);

  libMesh::out << "Building mesh" << std::endl;

  // Build a unit square
  ElemType elemtype;

  if (param.elementtype == "tri" ||
      param.elementtype == "unstructured")
    elemtype = TRI3;
  else
    elemtype = QUAD4;

  MeshTools::Generation::build_square (mesh, param.coarsegridx, param.coarsegridy,
                                       param.domain_xmin, param.domain_xmin + param.domain_edge_width,
                                       param.domain_ymin, param.domain_ymin + param.domain_edge_length,
                                       elemtype);

  libMesh::out << "Building system" << std::endl;

  HeatSystem & system = equation_systems.add_system<HeatSystem> ("HeatSystem");

  set_system_parameters(system, param);

  libMesh::out << "Initializing systems" << std::endl;

  // Initialize the system
  equation_systems.init ();

  // use theta from Stan Math
  // std::cout << "taki test: " << "set k" << "\n";
  system.set_k(theta[0]);
  system.get_equation_systems().parameters.set<Real>("_k") = theta[0];
  system.set_parameters(theta);

  // Refine the grid again if requested
  for (unsigned int i=0; i != param.extrarefinements; ++i)
    {
      mesh_refinement->uniformly_refine(1);
      equation_systems.reinit();
    }

  libMesh::out << "Setting primal initial conditions" << std::endl;

  read_initial_parameters();

  system.project_solution(initial_value, initial_grad,
                          equation_systems.parameters);

  // Output the H1 norm of the initial conditions
  libMesh::out << "|U("
               << system.time
               << ")|= "
               << system.calculate_norm(*system.solution, 0, H1)
               << std::endl
               << std::endl;

  // Add an adjoint vector, this will be computed after the forward
  // time stepping is complete
  //
  // Tell the library not to save adjoint solutions during the forward
  // solve
  //
  // Tell the library not to project this vector, and hence, memory
  // solution history to not save it.
  //
  // Make this vector ghosted so we can localize it to each element
  // later.
  const std::string & adjoint_solution_name = "adjoint_solution0";
  system.add_vector("adjoint_solution0", false, GHOSTED);

  // Close up any resources initial.C needed
  finish_initialization();

  // Plot the initial conditions
  write_output(equation_systems, 0, "primal", param);

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

  // we have one QoI and one parameter
  std::vector<std::vector<double> > result(1);

  // In optimized mode we catch any solver errors, so that we can
  // write the proper footers before closing.  In debug mode we just
  // let the exception throw so that gdb can grab it.
#ifdef NDEBUG
  try
    {
#endif
      // Now we begin the timestep loop to compute the time-accurate
      // solution of the equations.
      for (unsigned int t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
          // A pretty update message
          libMesh::out << " Solving time step "
                       << t_step
                       << ", time = "
                       << system.time
                       << std::endl;

          // Solve the forward problem at time t, to obtain the solution at time t + dt
          system.solve();

          // Output the H1 norm of the computed solution
          libMesh::out << "|U("
                       << system.time + system.deltat
                       << ")|= "
                       << system.calculate_norm(*system.solution, 0, H1)
                       << std::endl;

          // Advance to the next timestep in a transient problem
          libMesh::out << "Advancing timestep" << std::endl << std::endl;
          system.time_solver->advance_timestep();

          // Write out this timestep
          write_output(equation_systems, t_step+1, "primal", param);
        }
      // End timestep loop

      ///////////////// Now for the Adjoint Solution //////////////////////////////////////

      // Now we will solve the backwards in time adjoint problem
      libMesh::out << std::endl << "Solving the adjoint problem" << std::endl;

      // We need to tell the library that it needs to project the adjoint, so
      // MemorySolutionHistory knows it has to save it

      // Tell the library to project the adjoint vector, and hence, memory solution history to
      // save it
      system.set_vector_preservation(adjoint_solution_name, true);

      libMesh::out << "Setting adjoint initial conditions Z("
                   << system.time
                   << ")"
                   <<std::endl;

      // Need to call adjoint_advance_timestep once for the initial condition setup
      libMesh::out<<"Retrieving solutions at time t="<<system.time<<std::endl;
      system.time_solver->adjoint_advance_timestep();

      // Output the H1 norm of the retrieved solutions (u^i and u^i+1)
      libMesh::out << "|U("
                   << system.time + system.deltat
                   << ")|= "
                   << system.calculate_norm(*system.solution, 0, H1)
                   << std::endl;

      libMesh::out << "|U("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                   << std::endl;

      // The first thing we have to do is to apply the adjoint initial
      // condition. The user should supply these. Here they are specified
      // in the functions adjoint_initial_value and adjoint_initial_gradient
      system.project_vector(adjoint_initial_value,
                            adjoint_initial_grad,
                            equation_systems.parameters,
                            system.get_adjoint_solution(0));

      // Since we have specified an adjoint solution for the current
      // time (T), set the adjoint_already_solved boolean to true, so
      // we dont solve unnecessarily in the adjoint sensitivity method
      system.set_adjoint_already_solved(true);

      libMesh::out << "|Z("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(), 0, H1)
                   << std::endl
                   << std::endl;

      write_output(equation_systems, param.n_timesteps, "dual", param);

      // Now that the adjoint initial condition is set, we will start the
      // backwards in time adjoint integration

      // For loop stepping backwards in time
      for (unsigned int t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
          //A pretty update message
          libMesh::out << " Solving adjoint time step "
                       << t_step
                       << ", time = "
                       << system.time
                       << std::endl;

          // The adjoint_advance_timestep function calls the retrieve
          // function of the memory_solution_history class via the
          // memory_solution_history object we declared earlier.  The
          // retrieve function sets the system primal vectors to their
          // values at the current timestep.
          libMesh::out << "Retrieving solutions at time t=" << system.time << std::endl;
          system.time_solver->adjoint_advance_timestep();

          // Output the H1 norm of the retrieved solution
          libMesh::out << "|U("
                       << system.time + system.deltat
                       << ")|= "
                       << system.calculate_norm(*system.solution, 0, H1)
                       << std::endl;

          libMesh::out << "|U("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                       << std::endl;

          system.set_adjoint_already_solved(false);

          system.adjoint_solve();

          // Now that we have solved the adjoint, set the
          // adjoint_already_solved boolean to true, so we dont solve
          // unnecessarily in the error estimator
          system.set_adjoint_already_solved(true);

          libMesh::out << "|Z("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_adjoint_solution(), 0, H1)
                       << std::endl
                       << std::endl;

          // Get a pointer to the primal solution vector
          NumericVector<Number> & primal_solution = *system.solution;

          // Get a pointer to the solution vector of the adjoint problem for QoI 0
          NumericVector<Number> & dual_solution_0 = system.get_adjoint_solution(0);

          // Swap the primal and dual solutions so we can write out the adjoint solution
          primal_solution.swap(dual_solution_0);

          write_output(equation_systems, param.n_timesteps - (t_step + 1), "dual", param);

          // Swap back
          primal_solution.swap(dual_solution_0);

          // postprocessing: calc QoI
          system.postprocess_sides = true;
          system.postprocess();
        }
      // End adjoint timestep loop

      // Now that we have computed both the primal and adjoint solutions, we compute the sensitivities to the parameter p
      // dQ/dp = partialQ/partialp - partialR/partialp
      // partialQ/partialp = (Q(p+dp) - Q(p-dp))/(2*dp), this is not supported by the library yet
      // partialR/partialp = (R(u,z;p+dp) - R(u,z;p-dp))/(2*dp), where
      // R(u,z;p+dp) = int_{0}^{T} f(z;p+dp) - <partialu/partialt, z>(p+dp) - <g(u),z>(p+dp)
      // To do this we need to step forward in time, and compute the perturbed R at each time step and accumulate it
      // Then once all time steps are over, we can compute (R(u,z;p+dp) - R(u,z;p-dp))/(2*dp)

      // Now we begin the timestep loop to compute the time-accurate
      // adjoint sensitivities
      for (unsigned int t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
          // A pretty update message
          libMesh::out << "Retrieving "
                       << t_step
                       << ", time = "
                       << system.time
                       << std::endl;

          // Retrieve the primal and adjoint solutions at the current timestep
          system.time_solver->retrieve_timestep();

          libMesh::out << "|U("
                       << system.time + system.deltat
                       << ")|= "
                       << system.calculate_norm(*system.solution, 0, H1)
                       << std::endl;

          libMesh::out << "|U("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                       << std::endl;

          libMesh::out << "|Z("
                       << system.time
                       << ")|= "
                       << system.calculate_norm(system.get_adjoint_solution(0), 0, H1)
                       << std::endl
                       << std::endl;

          // Call the postprocess function which we have overloaded to compute
          // accumulate the perturbed residuals
          dynamic_cast<HeatSystem &>(system).perturb_accumulate_residuals(dynamic_cast<HeatSystem &>(system).get_parameter_vector());

          // Move the system time forward (retrieve_timestep does not do this)
          system.time += system.deltat;
        }

      // A pretty update message
      libMesh::out << "Retrieving final time = "
                   << system.time
                   << std::endl;

      // Retrieve the primal and adjoint solutions at the current timestep
      system.time_solver->retrieve_timestep();

      libMesh::out << "|U("
                   << system.time + system.deltat
                   << ")|= "
                   << system.calculate_norm(*system.solution, 0, H1)
                   << std::endl;

      libMesh::out << "|U("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1)
                   << std::endl;

      libMesh::out << "|Z("
                   << system.time
                   << ")|= "
                   << system.calculate_norm(system.get_adjoint_solution(0), 0, H1)
                   << std::endl
                   << std::endl;

      // Call the postprocess function which we have overloaded to compute
      // accumulate the perturbed residuals
      dynamic_cast<HeatSystem &>(system).perturb_accumulate_residuals(dynamic_cast<HeatSystem &>(system).get_parameter_vector());

      // Now that we computed the accumulated, perturbed residuals, we can compute the
      // approximate sensitivity
      Number sensitivity_0_0 = (dynamic_cast<HeatSystem &>(system)).compute_final_sensitivity();

      // Print it out
      libMesh::out << "Sensitivity of QoI 0 w.r.t parameter 0 is: "
                   << sensitivity_0_0
                   << std::endl;

      // QoI and QoI sensitivity
      result[0].push_back(system.get_QoI_value(0));
      result[0].push_back(sensitivity_0_0);

#ifdef NDEBUG
    }
  catch (...)
    {
      libMesh::err << '[' << mesh.processor_id()
                   << "] Caught exception; exiting early." << std::endl;
    }
#endif

  libMesh::err << '[' << mesh.processor_id()
               << "] Completing output."
               << std::endl;

  return result;
}

  // no sensitivity version is simply a wrapper
  inline std::vector<double>
  solve(const std::vector<double>& theta,
        const std::vector<double>& x_r,
        const std::vector<int>& x_i,
        std::ostream* msgs = nullptr) const {
    std::vector<std::vector<double> > sen_res =
      solve_with_sensitivity(theta, x_r, x_i, msgs);
    std::vector<double> res(sen_res.size());
    std::transform(sen_res.begin(), sen_res.end(),
                   res.begin(), [](std::vector<double>& v)
                   -> double { return v[0]; } );
    return res;
  }
};

TEST(LibMeshheatAdjointSensitivity, PDE_struct) {
  char arg0[] = "libmesh_heat_adjoint_sensitivity";
  char arg1[] = "arg1";
  char arg2[] = "arg2";
  char *argv[] = {arg0, arg1, arg2, NULL};
  int argc = sizeof(argv) / sizeof(char*) - 1;

  LibMeshHeatModel pde(argc, argv);

  const std::vector<double> theta{1.e-3};
  double res;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::vector<std::vector<double> > res_sen =
    pde.solve_with_sensitivity(theta, x_r, x_i);
  std::vector<double> res_prm =
    pde.solve(theta, x_r, x_i);
  res = 2.472073;
  ASSERT_FLOAT_EQ(res_sen[0][0], res);
  ASSERT_FLOAT_EQ(res_prm[0], res);
  res = -5.37173;
  ASSERT_FLOAT_EQ(res_sen[0][1], res);

  // test run MPI again
  std::vector<std::vector<double> > res2 =
    pde.solve_with_sensitivity(theta, x_r, x_i);
  ASSERT_FLOAT_EQ(res2[0][0], res_sen[0][0]);
  ASSERT_FLOAT_EQ(res2[0][1], res_sen[0][1]);
}

#endif
