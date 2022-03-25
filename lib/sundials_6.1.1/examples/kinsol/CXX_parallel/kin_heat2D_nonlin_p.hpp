/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 *                David Gardner @ LLNL
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
 * -----------------------------------------------------------------------------*/

#ifndef KIN_HEAT2D_NONLIN_P_HPP
#define KIN_HEAT2D_NONLIN_P_HPP

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>

#include "kinsol/kinsol.h"             // access to KINSOL
#include "nvector/nvector_parallel.h"  // access to the MPI N_Vector
#include "mpi.h"                       // MPI header file

// Macros for problem constants
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define HALF  RCONST(0.5)
#define TWO   RCONST(2.0)
#define EIGHT RCONST(8.0)

// Macro to access (x,y) location in 1D NVector array
#define IDX(x,y,n) ((n)*(y)+(x))

// Define c function type
typedef realtype (*cFn)(realtype u_val);

using namespace std;

// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------

struct UserData
{
  // Diffusion coefficients in the x and y directions
  realtype kx;
  realtype ky;

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

  // Receive buffers for neighbor exchange
  realtype *Wrecv;
  realtype *Erecv;
  realtype *Srecv;
  realtype *Nrecv;

  // Receive requests for neighbor exchange
  MPI_Request reqRW;
  MPI_Request reqRE;
  MPI_Request reqRS;
  MPI_Request reqRN;

  // Send buffers for neighbor exchange
  realtype *Wsend;
  realtype *Esend;
  realtype *Ssend;
  realtype *Nsend;

  // Send requests for neighor exchange
  MPI_Request reqSW;
  MPI_Request reqSE;
  MPI_Request reqSS;
  MPI_Request reqSN;

  // Fixed Point Solver settings
  realtype rtol;        // relative tolerance
  int      maa;         // m for Anderson Acceleration
  double   damping;     // daming for Anderson Acceleration
  int      orthaa;      // orthogonalization routine for AA
  int      maxits;      // max number of fixed point iterations

  // c(u) Function and integer for help setting
  cFn c;
  int c_int;

  // Vectors to help with FPFunction definition and execution
  N_Vector b;      // defined using c(u_exact)
  N_Vector vtemp;  // temporary vector for function evaluation

  // Ouput variables
  int      output; // output level
  N_Vector e;      // error vector
  ofstream uout;   // output file stream
  ofstream rout;   // output residual file stream
  ofstream eout;   // error file stream

  // Timing variables
  bool   timing;     // print timings
  double totaltime;
  double fevaltime;
  double exchangetime;
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// Nonlinear fixed point function
static int FPFunction(N_Vector u, N_Vector f, void *user_data);

// Nonlinear function c(u)
static int c(N_Vector u, N_Vector z, void *user_data);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Create rhs = b - c(u) for FPFunction
static int SetupRHS(void *user_data);

// Set nonlinear function c(u)
static int SetC(UserData *udata);

// Setup the parallel decomposition
static int SetupDecomp(MPI_Comm comm_w, UserData *udata);

// Perform neighbor exchange
static int PostRecv(UserData *udata);
static int SendData(N_Vector y, UserData *udata);
static int WaitRecv(UserData *udata);

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

// Compute the true solution
static int Solution(N_Vector u, UserData *udata);

// Compute the solution error solution
static int SolutionError(N_Vector u,  N_Vector e, UserData *udata);

// Print the command line options
static void InputHelp();

// Print some UserData information
static int PrintUserData(UserData *udata);

// Print Fixed Point statistics
static int OutputStats(void *kinsol_mem, UserData *udata);

// Print integration timing
static int OutputTiming(UserData *udata);

// Write solution to a file
static int WriteSolution(N_Vector u, UserData *udata);

// Functions for outputting residual history to file
static int OpenOutput(UserData *udata);
static int WriteOutput(N_Vector u, N_Vector f, UserData *udata);
static int CloseOutput(UserData *udata);

// Check function return values
static int check_retval(void *flagvalue, const string funcname, int opt);

// -----------------------------------------------------------------------------
// Multiple nonlinear functions for testing
// -----------------------------------------------------------------------------

// c(u) = u
realtype c1(realtype u_val)
{
  return u_val;
}

// c(u) = u^3 - u
realtype c2(realtype u_val)
{
  return u_val * u_val * u_val - u_val;
}

// c(u) = u - u^2
realtype c3(realtype u_val)
{
  return u_val - u_val * u_val;
}

// c(u) = e^u
realtype c4(realtype u_val)
{
  return exp(u_val);
}

// c(u) = u^4
realtype c5(realtype u_val)
{
  return u_val * u_val * u_val * u_val;
}

// c(u) = cos^2(u) - sin^2(u)
realtype c6(realtype u_val)
{
  return (cos(u_val) * cos(u_val)) - (sin(u_val) * sin(u_val));
}

// c(u) = cos^2(u) - sin^2(u) - e^u
realtype c7(realtype u_val)
{
  return (cos(u_val) * cos(u_val)) - (sin(u_val) * sin(u_val)) - exp(u_val);
}

// c(u) = e^u * u^4 - u * e^{cos(u)}
realtype c8(realtype u_val)
{
  realtype u2 = u_val * u_val;
  return exp(u_val) * u2 * u2 - u_val * exp(cos(u_val));
}

// c(u) = e^(cos^2(u))
realtype c9(realtype u_val)
{
  realtype cos2u = cos(u_val) * cos(u_val);
  return exp(cos2u);
}

// c(u) = 10(u - u^2)
realtype c10(realtype u_val)
{
  realtype u2 = u_val * u_val;
  return 10.0 * (u_val - u2);
}

// c(u) = -13 + u + ((5-u)u - 2)u
realtype c11(realtype u_val)
{
  realtype temp = ((5.0 - u_val) * u_val) - 2.0;
  return -13.0 + u_val + temp * u_val;
}

// c(u) = sqrt(5) * (u - u^2)
realtype c12(realtype u_val)
{
  realtype temp = sqrt(5);
  realtype u2 = u_val * u_val;
  return temp * (u_val - u2);
}

// c(u) = (u - e^u)^2 + (u + u * sin(u) - cos(u))^2
realtype c13(realtype u_val)
{
  realtype eu   = u_val - exp(u_val);
  realtype usin = u_val * sin(u_val);
  realtype temp = (u_val + usin - cos(u_val));
  return eu * eu + temp * temp;
}

// c(u) = u + ue^u + ue^{-u}
realtype c14(realtype u_val)
{
  realtype ueu  = u_val * exp(u_val);
  realtype ue_u = u_val * exp(-u_val);
  return u_val + ueu + ue_u;
}

// c(u) = u + ue^u + ue^{-u} + (u - e^u)^2
realtype c15(realtype u_val)
{
  realtype ueu  = u_val * exp(u_val);
  realtype ue_u = u_val * exp(-u_val);
  realtype temp = u_val - exp(u_val);
  return u_val + ueu + ue_u + (temp * temp);
}

// c(u) = u + ue^u + ue^{-u} + (u - e^u)^2 + (u + usin(u) - cos(u))^2
realtype c16(realtype u_val)
{
  realtype ueu   = u_val * exp(u_val);
  realtype ue_u  = u_val * exp(-u_val);
  realtype temp  = u_val - exp(u_val);
  realtype temp2 = u_val + (u_val * sin(u_val)) - cos(u_val);
  return u_val + ueu + ue_u + (temp * temp) + (temp2 * temp2);
}

// c(u) = u + ue^{-u} + e^u*(u + sin(u) - cos(u))^3
realtype c17(realtype u_val)
{
  realtype ue_u = u_val * exp(-u_val);
  realtype eu   = exp(u_val);
  realtype temp = u_val + sin(u_val) - cos(u_val);
  return u_val + ue_u + eu * (temp * temp * temp);
}

#endif
