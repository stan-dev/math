/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
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
 * Example problem:
 *
 * The following test simulates a simple anisotropic 2D heat equation,
 *
 *   u_t = kx u_xx + ky u_yy + b,
 *
 * for t in [0, 1] and (x,y) in [0, 1]^2, with initial conditions
 *
 *   u(0,x,y) = sin^2(pi x) sin^2(pi y),
 *
 * stationary boundary conditions
 *
 *   u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
 *
 * and the heat source
 *
 *   b(t,x,y) = -2 pi sin^2(pi x) sin^2(pi y) sin(pi t) cos(pi t)
 *              - kx 2 pi^2 (cos^2(pi x) - sin^2(pi x)) sin^2(pi y) cos^2(pi t)
 *              - ky 2 pi^2 (cos^2(pi y) - sin^2(pi y)) sin^2(pi x) cos^2(pi t).
 *
 * Under this setup, the problem has the analytical solution
 *
 *    u(t,x,y) = sin^2(pi x) sin^2(pi y) cos^2(pi t).
 *
 * The spatial derivatives are computed using second-order centered differences,
 * with the data distributed over nx * ny points on a uniform spatial grid. The
 * problem is advanced in time with BDF methods using an inexact Newton method
 * paired with the PCG or SPGMR linear solver. Several command line options are
 * available to change the problem parameters and CVODE settings. Use the flag
 * --help for more information.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <chrono>
#include <cmath>
#include <string>

#include "cvode/cvode.h"               // access to CVODE
#include "nvector/nvector_serial.h"    // access to the serial N_Vector
#include "sunlinsol/sunlinsol_pcg.h"   // access to PCG SUNLinearSolver
#include "sunlinsol/sunlinsol_spgmr.h" // access to SPGMR SUNLinearSolver


// Macros for problem constants
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define EIGHT RCONST(8.0)

// Macro to access (x,y) location in 1D NVector array
#define IDX(x,y,n) ((n)*(y)+(x))

using namespace std;

// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------

struct UserData
{
  // Diffusion coefficients in the x and y directions
  realtype kx;
  realtype ky;

  // Enable/disable forcing
  bool forcing;

  // Final time
  realtype tf;

  // Upper bounds in x and y directions
  realtype xu;
  realtype yu;

  // Number of nodes in the x and y directions
  sunindextype nx;
  sunindextype ny;

  // Total number of nodes
  sunindextype nodes;

  // Mesh spacing in the x and y directions
  realtype dx;
  realtype dy;

  // Integrator settings
  realtype rtol;        // relative tolerance
  realtype atol;        // absolute tolerance
  int      maxsteps;    // max number of steps between outputs

  // Linear solver and preconditioner settings
  bool     pcg;       // use PCG (true) or GMRES (false)
  bool     prec;      // preconditioner on/off
  bool     lsinfo;    // output residual history
  int      liniters;  // number of linear iterations
  int      msbp;      // max number of steps between preconditioner setups
  realtype epslin;    // linear solver tolerance factor

  // Inverse of Jacobian diagonal for preconditioner
  N_Vector d;

  // Ouput variables
  int      output; // output level
  int      nout;   // number of output times
  ofstream uout;   // output file stream
  ofstream eout;   // error file stream
  N_Vector e;      // error vector

  // Timing variables
  bool   timing;     // print timings
  double evolvetime;
  double rhstime;
  double psetuptime;
  double psolvetime;
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// ODE right hand side function
static int f(realtype t, N_Vector u, N_Vector f, void *user_data);

// Preconditioner setup and solve functions
static int PSetup(realtype t, N_Vector u, N_Vector f, booleantype jok,
                  booleantype *jcurPtr, realtype gamma, void *user_data);

static int PSolve(realtype t, N_Vector u, N_Vector f, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data);

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Set the default values in the UserData structure
static int InitUserData(UserData *udata);

// Free memory allocated within UserData
static int FreeUserData(UserData *udata);

// Read the command line inputs and set UserData values
static int ReadInputs(int *argc, char ***argv, UserData *udata);

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the true solution
static int Solution(realtype t, N_Vector u, UserData *udata);

// Compute the solution error solution
static int SolutionError(realtype t, N_Vector u,  N_Vector e, UserData *udata);

// Print the command line options
static void InputHelp();

// Print some UserData information
static int PrintUserData(UserData *udata);

// Output solution and error
static int OpenOutput(UserData *udata);
static int WriteOutput(realtype t, N_Vector u, UserData *udata);
static int CloseOutput(UserData *udata);

// Print integration statistics
static int OutputStats(void *cvode_mem, UserData *udata);

// Print integration timing
static int OutputTiming(UserData *udata);

// Check function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  sundials::Context sunctx;   // SUNDIALS context
  int flag;                   // reusable error-checking flag
  UserData *udata    = NULL;  // user data structure
  N_Vector u         = NULL;  // vector for storing solution
  SUNLinearSolver LS = NULL;  // linear solver memory structure
  void *cvode_mem    = NULL;  // CVODE memory structure
  FILE *diagfp       = NULL;  // diagnostics output file

  // Timing variables
  chrono::time_point<chrono::steady_clock> t1;
  chrono::time_point<chrono::steady_clock> t2;

  // ---------------
  // Setup UserData
  // ---------------

  // Allocate and initialize user data structure with default values. The
  // defaults may be overwritten by command line inputs in ReadInputs below.
  udata = new UserData;
  flag = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) return 1;

  // Parse command line inputs
  flag = ReadInputs(&argc, &argv, udata);
  if (flag != 0) return 1;

  // Output problem setup/options
  flag = PrintUserData(udata);
  if (check_flag(&flag, "PrintUserData", 1)) return 1;

  // Open diagnostics output file
  if (udata->lsinfo)
  {
    diagfp = fopen("diagnostics.txt", "w");
    if (check_flag((void *) diagfp, "fopen", 0)) return 1;
  }

  // ----------------------
  // Create serial vectors
  // ----------------------

  // Create vector for solution
  u = N_VNew_Serial(udata->nodes, sunctx);
  if (check_flag((void *) u, "N_VNew_Parallel", 0)) return 1;

  // Set initial condition
  flag = Solution(ZERO, u, udata);
  if (check_flag(&flag, "Solution", 1)) return 1;

  // Create vector for error
  udata->e = N_VClone(u);
  if (check_flag((void *) (udata->e), "N_VClone", 0)) return 1;

  // ---------------------
  // Create linear solver
  // ---------------------

  // Create linear solver
  int prectype = (udata->prec) ? SUN_PREC_RIGHT : SUN_PREC_NONE;

  if (udata->pcg)
  {
    LS = SUNLinSol_PCG(u, prectype, udata->liniters, sunctx);
    if (check_flag((void *) LS, "SUNLinSol_PCG", 0)) return 1;

    if (udata->lsinfo)
    {
      flag = SUNLinSolSetPrintLevel_PCG(LS, 1);
      if (check_flag(&flag, "SUNLinSolSetPrintLevel_PCG", 1)) return(1);

      flag = SUNLinSolSetInfoFile_PCG(LS, diagfp);
      if (check_flag(&flag, "SUNLinSolSetInfoFile_PCG", 1)) return(1);
    }
  }
  else
  {
    LS = SUNLinSol_SPGMR(u, prectype, udata->liniters, sunctx);
    if (check_flag((void *) LS, "SUNLinSol_SPGMR", 0)) return 1;

    if (udata->lsinfo)
    {
      flag = SUNLinSolSetPrintLevel_SPGMR(LS, 1);
      if (check_flag(&flag, "SUNLinSolSetPrintLevel_SPGMR", 1)) return(1);

      flag = SUNLinSolSetInfoFile_SPGMR(LS, diagfp);
      if (check_flag(&flag, "SUNLinSolSetInfoFile_SPGMR", 1)) return(1);
    }
  }

  // Allocate preconditioner workspace
  if (udata->prec)
  {
    udata->d = N_VClone(u);
    if (check_flag((void *) (udata->d), "N_VClone", 0)) return 1;
  }

  // --------------
  // Setup CVODE
  // --------------

  // Create integrator
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_flag((void *) cvode_mem, "CVodeCreate", 0)) return 1;

  // Initialize integrator
  flag = CVodeInit(cvode_mem, f, ZERO, u);
  if (check_flag(&flag, "CVodeInit", 1)) return 1;

  // Specify tolerances
  flag = CVodeSStolerances(cvode_mem, udata->rtol, udata->atol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return 1;

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, (void *) udata);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return 1;

  // Attach linear solver
  flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_flag(&flag, "CVodeSetLinearSolver", 1)) return 1;

  if (udata->prec)
  {
    // Attach preconditioner
    flag = CVodeSetPreconditioner(cvode_mem, PSetup, PSolve);
    if (check_flag(&flag, "CVodeSetPreconditioner", 1)) return 1;

    // Set linear solver setup frequency (update preconditioner)
    flag = CVodeSetLSetupFrequency(cvode_mem, udata->msbp);
    if (check_flag(&flag, "CVodeSetLSetupFrequency", 1)) return 1;
  }

  // Set linear solver tolerance factor
  flag = CVodeSetEpsLin(cvode_mem, udata->epslin);
  if (check_flag(&flag, "CVodeSetEpsLin", 1)) return 1;

  // Set max steps between outputs
  flag = CVodeSetMaxNumSteps(cvode_mem, udata->maxsteps);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return 1;

  // Set stopping time
  flag = CVodeSetStopTime(cvode_mem, udata->tf);
  if (check_flag(&flag, "CVodeSetStopTime", 1)) return 1;

  // -----------------------
  // Loop over output times
  // -----------------------

  realtype t     = ZERO;
  realtype dTout = udata->tf / udata->nout;
  realtype tout  = dTout;

  // Inital output
  flag = OpenOutput(udata);
  if (check_flag(&flag, "OpenOutput", 1)) return 1;

  flag = WriteOutput(t, u, udata);
  if (check_flag(&flag, "WriteOutput", 1)) return 1;

  for (int iout = 0; iout < udata->nout; iout++)
  {
    // Start timer
    t1 = chrono::steady_clock::now();

    // Evolve in time
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_flag(&flag, "CVode", 1)) break;

    // Stop timer
    t2 = chrono::steady_clock::now();

    // Update timer
    udata->evolvetime += chrono::duration<double>(t2 - t1).count();

    // Output solution and error
    flag = WriteOutput(t, u, udata);
    if (check_flag(&flag, "WriteOutput", 1)) return 1;

    // Update output time
    tout += dTout;
    tout = (tout > udata->tf) ? udata->tf : tout;
  }

  // Close output
  flag = CloseOutput(udata);
  if (check_flag(&flag, "CloseOutput", 1)) return 1;

  // --------------
  // Final outputs
  // --------------

  // Print final integrator stats
  if (udata->output > 0)
  {
    cout << "Final integrator statistics:" << endl;
    flag = OutputStats(cvode_mem, udata);
    if (check_flag(&flag, "OutputStats", 1)) return 1;
  }

  if (udata->forcing)
  {
    // Output final error
    flag = SolutionError(t, u, udata->e, udata);
    if (check_flag(&flag, "SolutionError", 1)) return 1;

    realtype maxerr = N_VMaxNorm(udata->e);

    cout << scientific;
    cout << setprecision(numeric_limits<realtype>::digits10);
    cout << "  Max error = " << maxerr << endl;
  }

  // Print timing
  if (udata->timing)
  {
    flag = OutputTiming(udata);
    if (check_flag(&flag, "OutputTiming", 1)) return 1;
  }

  // --------------------
  // Clean up and return
  // --------------------

  if (udata->lsinfo) fclose(diagfp);

  CVodeFree(&cvode_mem);     // Free integrator memory
  SUNLinSolFree(LS);         // Free linear solver
  N_VDestroy(u);             // Free vectors
  FreeUserData(udata);       // Free user data
  delete udata;
  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// f routine to compute the ODE RHS function f(t,y).
static int f(realtype t, N_Vector u, N_Vector f, void *user_data)
{
  // Timing variables
  chrono::time_point<chrono::steady_clock> t1;
  chrono::time_point<chrono::steady_clock> t2;

  // Start timer
  t1 = chrono::steady_clock::now();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to number of nodes
  sunindextype nx = udata->nx;
  sunindextype ny = udata->ny;

  // Constants for computing diffusion term
  realtype cx = udata->kx / (udata->dx * udata->dx);
  realtype cy = udata->ky / (udata->dy * udata->dy);
  realtype cc = -TWO * (cx + cy);

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  realtype *farray = N_VGetArrayPointer(f);
  if (check_flag((void *) farray, "N_VGetArrayPointer", 0)) return -1;

  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over domain interior and compute rhs forcing term
  if (udata->forcing)
  {
    realtype x, y;
    realtype sin_sqr_x, sin_sqr_y;
    realtype cos_sqr_x, cos_sqr_y;

    realtype bx = (udata->kx) * TWO * PI * PI;
    realtype by = (udata->ky) * TWO * PI * PI;

    realtype sin_t_cos_t = sin(PI * t) * cos(PI * t);
    realtype cos_sqr_t   = cos(PI * t) * cos(PI * t);

    for (sunindextype j = 1; j < ny - 1; j++)
    {
      for (sunindextype i = 1; i < nx - 1; i++)
      {
        x  = i * udata->dx;
        y  = j * udata->dy;

        sin_sqr_x = sin(PI * x) * sin(PI * x);
        sin_sqr_y = sin(PI * y) * sin(PI * y);

        cos_sqr_x = cos(PI * x) * cos(PI * x);
        cos_sqr_y = cos(PI * y) * cos(PI * y);

        farray[IDX(i,j,nx)] =
          -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t
          -bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t
          -by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
      }
    }
  }

  // Iterate over domain interior and add rhs diffusion term
  for (sunindextype j = 1; j < ny - 1; j++)
  {
    for (sunindextype i = 1; i < nx - 1; i++)
    {
      farray[IDX(i,j,nx)] +=
        cc * uarray[IDX(i,j,nx)]
        + cx * (uarray[IDX(i-1,j,nx)] + uarray[IDX(i+1,j,nx)])
        + cy * (uarray[IDX(i,j-1,nx)] + uarray[IDX(i,j+1,nx)]);
    }
  }

  // Stop timer
  t2 = chrono::steady_clock::now();

  // Update timer
  udata->rhstime += chrono::duration<double>(t2 - t1).count();

  // Return success
  return 0;
}

// Preconditioner setup routine
static int PSetup(realtype t, N_Vector u, N_Vector f, booleantype jok,
                  booleantype *jcurPtr, realtype gamma, void *user_data)
{
  // Timing variables
  chrono::time_point<chrono::steady_clock> t1;
  chrono::time_point<chrono::steady_clock> t2;

  // Start timer
  t1 = chrono::steady_clock::now();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Access data array
  realtype *diag = N_VGetArrayPointer(udata->d);
  if (check_flag((void *) diag, "N_VGetArrayPointer", 0)) return -1;

  // Constants for computing diffusion
  realtype cx = udata->kx / (udata->dx * udata->dx);
  realtype cy = udata->ky / (udata->dy * udata->dy);
  realtype cc = -TWO * (cx + cy);

  // Set all entries of d to the inverse diagonal values of interior
  // (since boundary RHS is 0, set boundary diagonals to the same)
  realtype c = ONE / (ONE - gamma * cc);
  N_VConst(c, udata->d);

  // Stop timer
  t2 = chrono::steady_clock::now();

  // Update timer
  udata->psetuptime += chrono::duration<double>(t2 - t1).count();

  // Return success
  return 0;
}

// Preconditioner solve routine for Pz = r
static int PSolve(realtype t, N_Vector u, N_Vector f, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data)
{
  // Timing variables
  chrono::time_point<chrono::steady_clock> t1;
  chrono::time_point<chrono::steady_clock> t2;

  // Start timer
  t1 = chrono::steady_clock::now();

  // Access user_data structure
  UserData *udata = (UserData *) user_data;

  // Perform Jacobi iteration
  N_VProd(udata->d, r, z);

  // Stop timer
  t2 = chrono::steady_clock::now();

  // Update timer
  udata->psolvetime += chrono::duration<double>(t2 - t1).count();

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
  // Diffusion coefficient
  udata->kx = ONE;
  udata->ky = ONE;

  // Enable forcing
  udata->forcing = true;

  // Final time
  udata->tf = ONE;

  // Upper bounds in x and y directions
  udata->xu = ONE;
  udata->yu = ONE;

  // Number of nodes in the x and y directions
  udata->nx    = 32;
  udata->ny    = 32;
  udata->nodes = udata->nx * udata->ny;

  // Mesh spacing in the x and y directions
  udata->dx = udata->xu / (udata->nx - 1);
  udata->dy = udata->yu / (udata->ny - 1);

  // Integrator settings
  udata->rtol     = RCONST(1.e-5);   // relative tolerance
  udata->atol     = RCONST(1.e-10);  // absolute tolerance
  udata->maxsteps = 0;               // use default

  // Linear solver and preconditioner options
  udata->pcg       = true;       // use PCG (true) or GMRES (false)
  udata->prec      = true;       // enable preconditioning
  udata->lsinfo    = false;      // output residual history
  udata->liniters  = 20;         // max linear iterations
  udata->msbp      = 0;          // use default (20 steps)
  udata->epslin    = ZERO;       // use default (0.05)

  // Inverse of Jacobian diagonal for preconditioner
  udata->d = NULL;

  // Output variables
  udata->output = 1;   // 0 = no output, 1 = stats output, 2 = output to disk
  udata->nout   = 20;  // Number of output times
  udata->e      = NULL;

  // Timing variables
  udata->timing     = false;
  udata->evolvetime = 0.0;
  udata->rhstime    = 0.0;
  udata->psetuptime = 0.0;
  udata->psolvetime = 0.0;

  // Return success
  return 0;
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{
  // Free preconditioner data
  if (udata->d)
  {
    N_VDestroy(udata->d);
    udata->d = NULL;
  }

  // Free error vector
  if (udata->e)
  {
    N_VDestroy(udata->e);
    udata->e = NULL;
  }

  // Return success
  return 0;
}

// Read command line inputs
static int ReadInputs(int *argc, char ***argv, UserData *udata)
{
  // Check for input args
  int arg_idx = 1;

  while (arg_idx < (*argc))
  {
    string arg = (*argv)[arg_idx++];

    // Mesh points
    if (arg == "--mesh")
    {
      udata->nx = stoi((*argv)[arg_idx++]);
      udata->ny = stoi((*argv)[arg_idx++]);
    }
    // Domain upper bounds
    else if (arg == "--domain")
    {
      udata->xu = stoi((*argv)[arg_idx++]);
      udata->yu = stoi((*argv)[arg_idx++]);
    }
    // Diffusion parameters
    else if (arg == "--k")
    {
      udata->kx = stod((*argv)[arg_idx++]);
      udata->ky = stod((*argv)[arg_idx++]);
    }
    // Disable forcing
    else if (arg == "--noforcing")
    {
      udata->forcing = false;
    }
    // Temporal domain settings
    else if (arg == "--tf")
    {
      udata->tf = stod((*argv)[arg_idx++]);
    }
    // Integrator settings
    else if (arg == "--rtol")
    {
      udata->rtol = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--atol")
    {
      udata->atol = stod((*argv)[arg_idx++]);
    }
    // Linear solver settings
    else if (arg == "--gmres")
    {
      udata->pcg = false;
    }
    else if (arg == "--lsinfo")
    {
      udata->lsinfo = true;
    }
    else if (arg == "--liniters")
    {
      udata->liniters = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--epslin")
    {
      udata->epslin = stod((*argv)[arg_idx++]);
    }
    // Preconditioner settings
    else if (arg == "--noprec")
    {
      udata->prec = false;
    }
    else if (arg == "--msbp")
    {
      udata->msbp = stoi((*argv)[arg_idx++]);
    }
    // Output settings
    else if (arg == "--output")
    {
      udata->output = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--nout")
    {
      udata->nout = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--maxsteps")
    {
      udata->maxsteps = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--timing")
    {
      udata->timing = true;
    }
    // Help
    else if (arg == "--help")
    {
      InputHelp();
      return -1;
    }
    // Unknown input
    else
    {
      cerr << "ERROR: Invalid input " << arg << endl;
      InputHelp();
      return -1;
    }
  }

  // Recompute total number of nodes
  udata->nodes = udata->nx * udata->ny;

  // Recompute x and y mesh spacing
  udata->dx = (udata->xu) / (udata->nx - 1);
  udata->dy = (udata->yu) / (udata->ny - 1);

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the exact solution
static int Solution(realtype t, N_Vector u, UserData *udata)
{
  realtype x, y;
  realtype cos_sqr_t;
  realtype sin_sqr_x, sin_sqr_y;

  // Constants for computing solution
  cos_sqr_t = cos(PI * t) * cos(PI * t);

  // Initialize u to zero (handles boundary conditions)
  N_VConst(ZERO, u);

  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  for (sunindextype j = 1; j < udata->ny - 1; j++)
  {
    for (sunindextype i = 1; i < udata->nx - 1; i++)
    {
      x = i * udata->dx;
      y = j * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uarray[IDX(i,j,udata->nx)] = sin_sqr_x * sin_sqr_y * cos_sqr_t;
    }
  }

  return 0;
}

// Compute the solution error
static int SolutionError(realtype t, N_Vector u, N_Vector e, UserData *udata)
{
  // Compute true solution
  int flag = Solution(t, e, udata);
  if (flag != 0) return -1;

  // Compute absolute error
  N_VLinearSum(ONE, u, -ONE, e, e);
  N_VAbs(e, e);

  return 0;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --mesh <nx> <ny>        : mesh points in the x and y directions" << endl;
  cout << "  --domain <xu> <yu>      : domain upper bound in the x and y direction" << endl;
  cout << "  --k <kx> <ky>           : diffusion coefficients" << endl;
  cout << "  --noforcing             : disable forcing term" << endl;
  cout << "  --tf <time>             : final time" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --atol <atol>           : absoltue tolerance" << endl;
  cout << "  --gmres                 : use GMRES linear solver" << endl;
  cout << "  --lsinfo                : output residual history" << endl;
  cout << "  --liniters <iters>      : max number of iterations" << endl;
  cout << "  --epslin <factor>       : linear tolerance factor" << endl;
  cout << "  --noprec                : disable preconditioner" << endl;
  cout << "  --msbp <steps>          : max steps between prec setups" << endl;
  cout << "  --output <level>        : output level" << endl;
  cout << "  --nout <nout>           : number of outputs" << endl;
  cout << "  --maxsteps <steps>      : max steps between outputs" << endl;
  cout << "  --timing                : print timing data" << endl;
  cout << "  --help                  : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData *udata)
{
  cout << endl;
  cout << "2D Heat PDE test problem:"                     << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  kx             = " << udata->kx              << endl;
  cout << "  ky             = " << udata->ky              << endl;
  cout << "  forcing        = " << udata->forcing         << endl;
  cout << "  tf             = " << udata->tf              << endl;
  cout << "  xu             = " << udata->xu              << endl;
  cout << "  yu             = " << udata->yu              << endl;
  cout << "  nx             = " << udata->nx              << endl;
  cout << "  ny             = " << udata->ny              << endl;
  cout << "  dx             = " << udata->dx              << endl;
  cout << "  dy             = " << udata->dy              << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  rtol           = " << udata->rtol            << endl;
  cout << "  atol           = " << udata->atol            << endl;
  cout << " --------------------------------- "           << endl;
  if (udata->pcg)
  {
    cout << "  linear solver  = PCG" << endl;
  }
  else
  {
    cout << "  linear solver  = GMRES" << endl;
  }
  cout << "  lin iters      = " << udata->liniters        << endl;
  cout << "  eps lin        = " << udata->epslin          << endl;
  cout << "  prec           = " << udata->prec            << endl;
  cout << "  msbp           = " << udata->msbp            << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  output         = " << udata->output          << endl;
  cout << " --------------------------------- "           << endl;
  cout << endl;

  return 0;
}

// Initialize output
static int OpenOutput(UserData *udata)
{
  // Header for status output
  if (udata->output > 0)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<realtype>::digits10);
    if (udata->forcing)
    {
      cout << "          t           ";
      cout << "          ||u||_rms      ";
      cout << "          max error      " << endl;
      cout << " ---------------------";
      cout << "-------------------------";
      cout << "-------------------------" << endl;
    }
    else
    {
      cout << "          t           ";
      cout << "          ||u||_rms      " << endl;
      cout << " ---------------------";
      cout << "-------------------------" << endl;
    }
  }

  // Output problem information and open output streams
  if (udata->output == 2)
  {
    // Each processor outputs subdomain information
    ofstream dout;
    dout.open("heat2d_info.txt");
    dout <<  "xu  " << udata->xu       << endl;
    dout <<  "yu  " << udata->yu       << endl;
    dout <<  "nx  " << udata->nx       << endl;
    dout <<  "ny  " << udata->ny       << endl;
    dout <<  "nt  " << udata->nout + 1 << endl;
    dout.close();

    // Open output streams for solution and error
    udata->uout.open("heat2d_solution.txt");
    udata->uout << scientific;
    udata->uout << setprecision(numeric_limits<realtype>::digits10);

    if (udata->forcing)
    {
      udata->eout.open("heat2d_error.txt");
      udata->eout << scientific;
      udata->eout << setprecision(numeric_limits<realtype>::digits10);
    }
  }

  return 0;
}

// Write output
static int WriteOutput(realtype t, N_Vector u, UserData *udata)
{
  int      flag;
  realtype max;

  if (udata->output > 0)
  {
    if (udata->forcing)
    {
      // Compute the error
      flag = SolutionError(t, u, udata->e, udata);
      if (check_flag(&flag, "SolutionError", 1)) return 1;

      // Compute max error
      max = N_VMaxNorm(udata->e);
    }

    // Compute rms norm of the state
    realtype urms = sqrt(N_VDotProd(u, u) / udata->nx / udata->ny);

    // Output current status
    if (udata->forcing)
    {
      cout << setw(22) << t << setw(25) << urms << setw(25) << max << endl;
    }
    else
    {
      cout << setw(22) << t << setw(25) << urms << endl;
    }

    // Write solution and error to disk
    if (udata->output == 2)
    {
      realtype *uarray = N_VGetArrayPointer(u);
      if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

      udata->uout << t << " ";
      for (sunindextype i = 0; i < udata->nodes; i++)
      {
        udata->uout << uarray[i] << " ";
      }
      udata->uout << endl;

      if (udata->forcing)
      {
        // Output error to disk
        realtype *earray = N_VGetArrayPointer(udata->e);
        if (check_flag((void *) earray, "N_VGetArrayPointer", 0)) return -1;

        udata->eout << t << " ";
        for (sunindextype i = 0; i < udata->nodes; i++)
        {
          udata->eout << earray[i] << " ";
        }
        udata->eout << endl;
      }
    }
  }

  return 0;
}

// Finalize output
static int CloseOutput(UserData *udata)
{
  // Footer for status output
  if (udata->output > 0)
  {
    if (udata->forcing)
    {
      cout << " ---------------------";
      cout << "-------------------------";
      cout << "-------------------------" << endl;
      cout << endl;
    }
    else
    {
      cout << " ---------------------";
      cout << "-------------------------" << endl;
      cout << endl;
    }
  }

  if (udata->output == 2)
  {
    // Close output streams
    udata->uout.close();
    if (udata->forcing) udata->eout.close();
  }

  return 0;
}

// Print integrator statistics
static int OutputStats(void *cvode_mem, UserData* udata)
{
  int flag;

  // Get integrator and solver stats
  long int nst, netf, nf, nni, ncfn, nli, nlcf, nsetups, nf_ls, nJv;
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  if (check_flag(&flag, "CVodeGetNumSteps", 1)) return -1;
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  if (check_flag(&flag, "CVodeGetNumErrTestFails", 1)) return -1;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nf);
  if (check_flag(&flag, "CVodeGetNumRhsEvals", 1)) return -1;
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  if (check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1)) return -1;
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  if (check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1)) return -1;
  flag = CVodeGetNumLinIters(cvode_mem, &nli);
  if (check_flag(&flag, "CVodeGetNumLinIters", 1)) return -1;
  flag = CVodeGetNumLinConvFails(cvode_mem, &nlcf);
  if (check_flag(&flag, "CVodeGetNumLinConvFails", 1)) return -1;
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  if (check_flag(&flag, "CVodeGetNumLinSolvSetups", 1)) return -1;
  flag = CVodeGetNumLinRhsEvals(cvode_mem, &nf_ls);
  if (check_flag(&flag, "CVodeGetNumLinRhsEvals", 1)) return -1;
  flag = CVodeGetNumJtimesEvals(cvode_mem, &nJv);
  if (check_flag(&flag, "CVodeGetNumJtimesEvals", 1)) return -1;

  cout << fixed;
  cout << setprecision(6);

  cout << "  Steps            = " << nst     << endl;
  cout << "  Error test fails = " << netf    << endl;
  cout << "  RHS evals        = " << nf      << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << "  NLS fails        = " << ncfn    << endl;
  cout << "  LS iters         = " << nli     << endl;
  cout << "  LS fails         = " << nlcf    << endl;
  cout << "  LS setups        = " << nsetups << endl;
  cout << "  LS RHS evals     = " << nf_ls   << endl;
  cout << "  Jv products      = " << nJv     << endl;
  cout << endl;

  // Compute average nls iters per step attempt and ls iters per nls iter
  realtype avgnli = (realtype) nni / (realtype) nst;
  realtype avgli  = (realtype) nli / (realtype) nni;
  cout << "  Avg NLS iters per step    = " << avgnli << endl;
  cout << "  Avg LS iters per NLS iter = " << avgli  << endl;
  cout << endl;

  // Get preconditioner stats
  if (udata->prec)
  {
    long int npe, nps;
    flag = CVodeGetNumPrecEvals(cvode_mem, &npe);
    if (check_flag(&flag, "CVodeGetNumPrecEvals", 1)) return -1;
    flag = CVodeGetNumPrecSolves(cvode_mem, &nps);
    if (check_flag(&flag, "CVodeGetNumPrecSolves", 1)) return -1;

    cout << "  Preconditioner setups = " << npe << endl;
    cout << "  Preconditioner solves = " << nps << endl;
    cout << endl;
  }

  return 0;
}

static int OutputTiming(UserData *udata)
{
  cout << scientific;
  cout << setprecision(6);

  cout << "  Evolve time = " << udata->evolvetime << " sec" << endl;
  cout << "  RHS time    = " << udata->rhstime    << " sec" << endl;
  cout << endl;

  if (udata->prec)
  {
    cout << "  PSetup time = " << udata->psetuptime << " sec" << endl;
    cout << "  PSolve time = " << udata->psolvetime << " sec" << endl;
    cout << endl;
  }

  return 0;
}

// Check function return value
static int check_flag(void *flagvalue, const string funcname, int opt)
{
  // Check if the function returned a NULL pointer
  if (opt == 0)
  {
    if (flagvalue == NULL)
    {
      cerr << endl << "ERROR: " << funcname << " returned NULL pointer" << endl
           << endl;
      return 1;
    }
  }
  // Check the function return flag value
  else if (opt == 1 || opt == 2)
  {
    int errflag = *((int *) flagvalue);
    if  ((opt == 1 && errflag < 0) || (opt == 2 && errflag != 0))
    {
      cerr << endl << "ERROR: " << funcname << " returned with flag = "
           << errflag << endl << endl;
      return 1;
    }
  }
  else
  {
    cerr << endl << "ERROR: check_flag called with an invalid option value"
         << endl;
    return 1;
  }

  return 0;
}

//---- end of file ----
