/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ UIUC/LLNL
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

#ifndef KIN_EM_MPICUDA_H
#define KIN_EM_MPICUDA_H

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <random>
#include <algorithm>

#include "kinsol/kinsol.h"             // access to KINSOL
#include <nvector/nvector_cuda.h>      // access to the cuda N_Vector
#include <nvector/nvector_mpiplusx.h>  // to be used with the MPI + X N_Vector
#include "mpi.h"                       // MPI header file

// Macros for problem constants
#define PI      RCONST(3.141592653589793238462643383279502884197169)
#define ZERO    RCONST(0.0)
#define ONE     RCONST(1.0)
#define HALF    RCONST(0.5)
#define PTTHREE RCONST(0.3)
#define PTFOUR  RCONST(0.4)
#define TWO     RCONST(2.0)
#define FIVE    RCONST(5.0)
#define TEN     RCONST(10.0)

// Maximum size of output directory string
#define MXSTR 2048

// Macro to access (x,mu) location in 1D NVector array
#define IDX(x,mu) ((3)*(x)+(mu))

using namespace std;

// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------

struct UserData
{
  // Sigmas
  realtype sigma;

  // Alphas - mixture proportions
  realtype alpha1;
  realtype alpha2;
  realtype alpha3;

  // Global total number of nodes
  sunindextype nodes;

  // Overall number of local nodes
  sunindextype nodes_loc;

  // MPI variables
  MPI_Comm comm; // Communicator in space
  int nprocs_w;  // total number of MPI processes in Comm world
  int myid;      // process ID in communicator

  // Fixed Point Solver settings
  realtype rtol;        // relative tolerance
  int      maa;         // m for Anderson Acceleration
  double   damping;     // daming for Anderson Acceleration
  int      orthaa;      // orthogonalization routine for AA
  int      maxits;      // max number of fixed point iterations

  // Vectors to help with FPFunction definition and execution
  N_Vector samples_local; // vector containing distribution samples
  N_Vector px;            // temporary vector
  N_Vector mu_bottom;     // temporary vector
  N_Vector mu_top;
  N_Vector mu_true;       // vector of true means

  int num_samples;

  // Ouput variables
  int      output; // output level
  N_Vector vtemp;  // error vector
  ofstream uout;   // output file stream
  ofstream rout;   // output residual file stream
  ofstream eout;   // error file stream

  // Timing variables
  bool   timing;     // print timings
  double totaltime;
  double fevaltime;

  bool debug;
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// Nonlinear fixed point function
static int FPFunction(N_Vector u, N_Vector f, void *user_data);

// Expectation Maximization Algorithm
static int EM(N_Vector u, N_Vector f, void *user_data);

// Setup up mean distribution samples
static int SetupSamples(UserData *udata);

// Random Vector
static int SetMus(UserData *udata);

// Starting Vector
static int SetStartGuess(N_Vector u, UserData *udata);

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

// Compute the solution error solution
static int SolutionError(N_Vector u_true, N_Vector u, N_Vector e,
                         UserData *udata);

// Print the command line options
static void InputHelp();

// Print some UserData information
static int PrintUserData(UserData *udata);

// Print Fixed Point statistics
static int OutputStats(void *kinsol_mem, UserData *udata);

// Print integration timing
static int OutputTiming(UserData *udata);

// Functions for outputting residual history to file
static int OpenOutput(UserData *udata);
static int WriteOutput(N_Vector u, N_Vector f, UserData *udata);
static int CloseOutput(UserData *udata);

// Check function return values
static int check_retval(void *flagvalue, const string funcname, int opt);

#endif
