/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * -----------------------------------------------------------------------------
 * Shared header file for 2D diffusion benchmark problem
 * ---------------------------------------------------------------------------*/

#ifndef DIFFUSION_2D_HPP
#define DIFFUSION_2D_HPP

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>

#include "mpi.h"

#if defined(USE_HIP)
#include "nvector/nvector_mpiplusx.h"
#include "nvector/nvector_hip.h"
#elif defined(USE_CUDA)
#include "nvector/nvector_mpiplusx.h"
#include "nvector/nvector_cuda.h"
#else
#include "nvector/nvector_parallel.h"
#endif

#include "sunlinsol/sunlinsol_pcg.h"
#include "sunlinsol/sunlinsol_spgmr.h"

// Macros for problem constants
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define EIGHT RCONST(8.0)

// Macro to access (x,y) location in 1D NVector array
#define IDX(x,y,n) ((n)*(y)+(x))

// Ceiling for integers ceil(a/b) = ceil((a + b - 1) / b)
#define ICEIL(a,b) (((a) + (b) - 1) / (b))

using namespace std;

// -----------------------------------------------------------------------------
// UserData structure
// -----------------------------------------------------------------------------

struct UserData
{
  SUNProfiler prof = NULL;

  // Diffusion coefficients in the x and y directions
  realtype kx = ONE;
  realtype ky = ONE;

  // Enable/disable forcing
  bool forcing = true;

  // Final time
  realtype tf = ONE;

  // Lower and Upper bounds in x and y directions
  realtype xl = ZERO;
  realtype yl = ZERO;
  realtype xu = ONE;
  realtype yu = ONE;

  // Global number of nodes in the x and y directions
  sunindextype nx = 32;
  sunindextype ny = 32;

  // Global total number of nodes
  sunindextype nodes = nx * ny;

  // Mesh spacing in the x and y directions
  realtype dx = xu / (nx - 1);
  realtype dy = yu / (ny - 1);

  // Local number of nodes in the x and y directions
  sunindextype nx_loc = 0;
  sunindextype ny_loc = 0;

  // Overall number of local nodes
  sunindextype nodes_loc = 0;

  // Global x and y indices of this subdomain
  sunindextype is = 0;  // x starting index
  sunindextype ie = 0;  // x ending index
  sunindextype js = 0;  // y starting index
  sunindextype je = 0;  // y ending index

  // MPI variables
  MPI_Comm comm_c = MPI_COMM_NULL; // Cartesian communicator in space

  int np  = 1; // total number of MPI processes in Comm world
  int npx = 0; // number of MPI processes in the x-direction
  int npy = 0; // number of MPI processes in the y-direction

  int myid_c = 0; // process ID in Cartesian communicator

  // Flags denoting if this process has a neighbor
  bool HaveNbrW = true;
  bool HaveNbrE = true;
  bool HaveNbrS = true;
  bool HaveNbrN = true;

  // Neighbor IDs for exchange
  int ipW = -1;
  int ipE = -1;
  int ipS = -1;
  int ipN = -1;

  // Receive buffers for neighbor exchange
  realtype *Wrecv = NULL;
  realtype *Erecv = NULL;
  realtype *Srecv = NULL;
  realtype *Nrecv = NULL;

  // Receive requests for neighbor exchange
  MPI_Request reqRW;
  MPI_Request reqRE;
  MPI_Request reqRS;
  MPI_Request reqRN;

  // Send buffers for neighbor exchange
  realtype *Wsend = NULL;
  realtype *Esend = NULL;
  realtype *Ssend = NULL;
  realtype *Nsend = NULL;

  // Send requests for neighor exchange
  MPI_Request reqSW;
  MPI_Request reqSE;
  MPI_Request reqSS;
  MPI_Request reqSN;

  // Inverse of Jacobian diagonal for preconditioner
  N_Vector diag = NULL;

  UserData(SUNProfiler prof) : prof(prof) {}
  ~UserData();

  // Helper functions
  int parse_args(vector<string> &args, bool outproc);
  void help();
  void print();
  int setup();
  int start_exchange(const N_Vector u);
  int end_exchange();

private:
  int allocate_buffers();
  int pack_buffers(const N_Vector u);
  int free_buffers();
};

// -----------------------------------------------------------------------------
// UserOutput structure
// -----------------------------------------------------------------------------

struct UserOutput
{
  // Ouput variables
  int      output = 1;    // 0 = no output, 1 = stats output, 2 = output to disk
  int      nout   = 20;   // number of output times
  N_Vector error  = NULL; // error vector
  ofstream uoutstream;    // output file stream
  ofstream eoutstream;    // error file stream

  // Helper functions
  int parse_args(vector<string> &args, bool outproc);
  void help();
  void print();

  // Output functions
  int open(UserData* udata);
  int write(realtype t, N_Vector u, UserData *udata);
  int close(UserData* udata);
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// Common function for computing the Laplacian
int laplacian(realtype t, N_Vector u, N_Vector f, void *user_data);

#if defined(BENCHMARK_ODE)

// ODE right hand side function
int diffusion(realtype t, N_Vector u, N_Vector f, void *user_data);

// Preconditioner setup and solve functions
int PSetup(realtype t, N_Vector u, N_Vector f, booleantype jok,
           booleantype *jcurPtr, realtype gamma, void *user_data);

int PSolve(realtype t, N_Vector u, N_Vector f, N_Vector r,
           N_Vector z, realtype gamma, realtype delta, int lr,
           void *user_data);

#elif defined(BENCHMARK_DAE)

// DAE residual function
int diffusion(realtype t, N_Vector u, N_Vector up, N_Vector res,
              void *user_data);

// Preconditioner setup and solve functions
int PSetup(realtype t, N_Vector u, N_Vector up, N_Vector res, realtype cj,
           void *user_data);

int PSolve(realtype t, N_Vector u, N_Vector up, N_Vector res, N_Vector r,
           N_Vector z, realtype cj, realtype delta, void *user_data);

#else
#error "Missing ODE/DAE preprocessor directive"
#endif

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Synchronize devices
int DeviceSynchronize();

// Copy data from devices
int CopyDataFromDevice(N_Vector y);

// Compute the true solution and its derivative
int Solution(realtype t, N_Vector u, UserData *udata);
int SolutionDerivative(realtype t, N_Vector up, UserData *udata);

// Compute the solution error
int SolutionError(realtype t, N_Vector u, N_Vector e, UserData *udata);

// Check function return values
int check_flag(void *flagvalue, const string funcname, int opt);

#endif
