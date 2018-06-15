
#ifndef FEMPARAMETERS_H
#define FEMPARAMETERS_H

#include "libmesh/libmesh_common.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/function_base.h"
#include "libmesh/getpot.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"

#include <limits>
#include <map>
#include <string>
#include <vector>



class FEMParameters : public libMesh::ParallelObject
{
public:
  FEMParameters(const libMesh::Parallel::Communicator & comm_in);

  ~FEMParameters();

  void read(GetPot & input,
            const std::vector<std::string> * other_variable_names = libmesh_nullptr);

  // Parameters applicable to entire EquationSystems:

  unsigned int initial_timestep, n_timesteps;
  bool transient;
  unsigned int deltat_reductions;
  std::string timesolver_core;
  libMesh::Real end_time, deltat, timesolver_theta,
    timesolver_maxgrowth, timesolver_tolerance,
    timesolver_upper_tolerance, steadystate_tolerance;
  std::vector<libMesh::FEMNormType> timesolver_norm;

  //   Mesh generation

  unsigned int dimension;
  std::string domaintype, domainfile, elementtype;
  libMesh::Real elementorder;
  libMesh::Real domain_xmin, domain_ymin, domain_zmin;
  libMesh::Real domain_edge_width, domain_edge_length, domain_edge_height;
  unsigned int coarsegridx, coarsegridy, coarsegridz;
  unsigned int coarserefinements, extrarefinements;
  std::string mesh_redistribute_func;

  //   Mesh partitioning
  std::string mesh_partitioner_type;

  //   Mesh refinement

  unsigned int nelem_target;
  libMesh::Real global_tolerance;
  libMesh::Real refine_fraction, coarsen_fraction, coarsen_threshold;
  unsigned int max_adaptivesteps;
  unsigned int initial_adaptivesteps;

  //   Output

  unsigned int write_interval;
  bool write_gmv_error, write_tecplot_error,
    write_exodus_error,
    output_xda, output_xdr,
    output_bz2, output_gz,
    output_gmv, output_tecplot,
    output_exodus, output_nemesis;

  // Types of Systems to create

  std::vector<std::string> system_types;

  // Parameters applicable to each system:

  //   Boundary and initial conditions

#ifdef LIBMESH_ENABLE_PERIODIC
  std::vector<libMesh::PeriodicBoundary> periodic_boundaries;
#endif

  std::map<libMesh::subdomain_id_type, libMesh::FunctionBase<libMesh::Number> *>
  initial_conditions;
  std::map<libMesh::boundary_id_type, libMesh::FunctionBase<libMesh::Number> *>
  dirichlet_conditions,
    neumann_conditions;
  std::map<libMesh::boundary_id_type, std::vector<unsigned int>>
  dirichlet_condition_variables,
    neumann_condition_variables;
  std::map<int, std::map<libMesh::subdomain_id_type,
                         libMesh::FunctionBase<libMesh::Number> *>>
  other_interior_functions;
  std::map<int, std::map<libMesh::boundary_id_type,
                         libMesh::FunctionBase<libMesh::Number> *>>
  other_boundary_functions;

  //   Execution type

  bool run_simulation, run_postprocess;

  //   Approximation type

  std::vector<std::string> fe_family;
  std::vector<unsigned int> fe_order;
  int extra_quadrature_order;

  //   Assembly options

  bool analytic_jacobians;
  libMesh::Real verify_analytic_jacobians;
  libMesh::Real numerical_jacobian_h;

  bool print_solution_norms, print_solutions,
    print_residual_norms, print_residuals,
    print_jacobian_norms, print_jacobians,
    print_element_solutions,
    print_element_residuals,
    print_element_jacobians;

  //   Solver options

  bool use_petsc_snes;
  bool time_solver_quiet, solver_quiet, solver_verbose,
    reuse_preconditioner, require_residual_reduction;
  libMesh::Real min_step_length;
  unsigned int max_linear_iterations, max_nonlinear_iterations;
  libMesh::Real relative_step_tolerance, relative_residual_tolerance,
    absolute_residual_tolerance,
    initial_linear_tolerance, minimum_linear_tolerance,
    linear_tolerance_multiplier;

  // Initialization

  unsigned int initial_sobolev_order;
  unsigned int initial_extra_quadrature;

  //   Error indicators

  bool refine_uniformly;
  std::string indicator_type;
  bool patch_reuse;
  unsigned int sobolev_order;

  //   System-specific parameters file

  std::string system_config_file;
};

#endif // FEMPARAMETERS_H
