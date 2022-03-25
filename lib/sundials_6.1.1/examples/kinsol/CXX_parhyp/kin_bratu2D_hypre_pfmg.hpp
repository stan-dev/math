/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------*/

#ifndef KIN_BRATU2D_NONLIN_HYPRE_PFMG_HPP
#define KIN_BRATU2D_NONLIN_HYPRE_PFMG_HPP

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>

#include "kinsol/kinsol.h"                  // access to KINSOL
#include "nvector/nvector_parallel.h"       // access to the MPI N_Vector
#include "sunlinsol/sunlinsol_pcg.h"        // access to PCG SUNLinearSolver
#include "HYPRE_struct_ls.h"                // HYPRE structured grid solver interface
#include "mpi.h"                            // MPI header file

// Macros for problem constants
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define SIX   RCONST(6.0)
#define EIGHT RCONST(8.0)

// Macro to access (x,y) location in 1D NVector array
#define IDX(x,y,n) ((n)*(y)+(x))

using namespace std;

// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------

struct UserData
{
  // Exponential term coefficient
  realtype C;

  // Upper bounds in x and y directions
  realtype xu;
  realtype yu;

  // Global number of nodes in the x and y directions
  sunindextype nx;
  sunindextype ny;

  // Global total number of nodes
  sunindextype nodes;

  // Mesh spacing in the x and y directions
  realtype dx;
  realtype dy;

  // Local number of nodes in the x and y directions
  sunindextype nx_loc;
  sunindextype ny_loc;

  // Overall number of local nodes
  sunindextype nodes_loc;

  // Global x and y indices of this subdomain
  sunindextype is;  // x starting index
  sunindextype ie;  // x ending index
  sunindextype js;  // y starting index
  sunindextype je;  // y ending index

  // MPI variables
  MPI_Comm comm_c; // Cartesian communicator in space

  int nprocs_w; // total number of MPI processes in Comm world
  int npx;      // number of MPI processes in the x-direction
  int npy;      // number of MPI processes in the y-direction

  int myid_c; // process ID in Cartesian communicator

  // Flags denoting if this process has a neighbor
  bool HaveNbrW;
  bool HaveNbrE;
  bool HaveNbrS;
  bool HaveNbrN;

  // Neighbor IDs for exchange
  int ipW;
  int ipE;
  int ipS;
  int ipN;

  // Fixed Point Solver settings
  realtype rtol;        // relative tolerance
  int      maa;         // m for Anderson Acceleration
  realtype damping;     // daming for Anderson Acceleration
  int      orthaa;      // orthogonalization routine for AA
  int      maxits;      // max number of fixed point iterations

  // Linear solver and preconditioner settings
  bool     lsinfo;    // output residual history
  int      liniters;  // number of linear iterations
  int      msbp;      // max number of steps between preconditioner setups
  realtype epslin;    // linear solver tolerance factor

  // Linear solver object
  SUNLinearSolver LS;  // linear solver memory structure

  // hypre objects
  HYPRE_StructGrid    grid;
  HYPRE_StructStencil stencil;
  HYPRE_StructMatrix  Jmatrix;
  HYPRE_StructVector  bvec;
  HYPRE_StructVector  xvec;
  HYPRE_StructVector  vvec;
  HYPRE_StructVector  Jvvec;
  HYPRE_StructSolver  precond;

  // hypre grid extents
  HYPRE_Int ilower[2];
  HYPRE_Int iupper[2];

  // hypre workspace
  HYPRE_Int   nwork;
  HYPRE_Real *work;

  // hypre counters
  HYPRE_Int pfmg_its;

  // hypre PFMG settings (hypre defaults)
  HYPRE_Int pfmg_relax;  // type of relaxation:
                         //   0 - Jacobi
                         //   1 - Weighted Jacobi
                         //   2 - symmetric R/B Gauss-Seidel (*)
                         //   3 - nonsymmetric R/B Gauss-Seidel
  HYPRE_Int pfmg_nrelax; // number of pre and post relaxation sweeps (2)

  // Ouput variables
  int      output; // output level
  ofstream uout;   // output file stream
  ofstream rout;   // output residual file stream
  N_Vector e;      // error vector

  // Timing variables
  bool   timing;     // print timings
  double totaltime;
  double fevaltime;
  double matfilltime;
  double jvtime;
  double psetuptime;
  double psolvetime;
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS iterator
// -----------------------------------------------------------------------------

// Nonlinear fixed point function
static int FPFunction(N_Vector u, N_Vector f, void *user_data);

// Jacobian-vector product function
static int JTimes(void *user_data, N_Vector v, N_Vector Jv);

// Preconditioner setup and solve functions
static int PSetup(void *user_data);

static int PSolve(void *user_data, N_Vector r, N_Vector z,
                  realtype tol, int lr);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Setup the parallel decomposition
static int SetupDecomp(MPI_Comm comm_w, UserData *udata);

// Create hypre objects
static int SetupHypre(UserData *udata);

// Create PCG Linear Solver and attach hypre
static int SetupLS(N_Vector u, void *user_data, SUNContext sunctx);

// Fill Jacobian and A = I - gamma * J
static int Jac(UserData *udata);

// Initial Guess
static int InitialGuess(N_Vector u, UserData *udata);

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Set the default values in the UserData structure
static int InitUserData(UserData *udata);

// Free memory allocated within UserData
static int FreeUserData(UserData *udata);

// Read the command line inputs and set UserData values
static int ReadInputs(int *argc, char ***argv, UserData *udata, bool outproc);

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Print the command line options
static void InputHelp();

// Print some UserData information
static int PrintUserData(UserData *udata);

// Print solver statistics
static int OutputStats(void *kinsol_mem, UserData *udata);

// Print solver timing
static int OutputTiming(UserData *udata);

// Write solution to a file
static int WriteSolution(N_Vector u, UserData *udata);

// Functions for outputting residual history to file
static int OpenResOutput(UserData *udata);
static int WriteResOutput(UserData *udata);
static int CloseResOutput(UserData *udata);

// Check function return values
static int check_retval(void *flagvalue, const string funcname, int opt);

#endif
