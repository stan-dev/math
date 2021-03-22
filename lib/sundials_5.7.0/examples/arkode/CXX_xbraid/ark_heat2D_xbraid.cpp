/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
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
 * problem is solved using the XBraid multigrid reduction in time library paired
 * with a diagonally implicit Runge-Kutta method from the ARKODE ARKStep module
 * using an inexact Newton method paired with the PCG or SPGMR linear solver.
 * Several command line options are available to change the problem parameters
 * and ARKStep settings. Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <chrono>
#include <cmath>

#include "arkode/arkode_arkstep.h"     // access to ARKStep
#include "nvector/nvector_serial.h"    // access to the serial N_Vector
#include "sunlinsol/sunlinsol_pcg.h"   // access to PCG SUNLinearSolver
#include "sunlinsol/sunlinsol_spgmr.h" // access to SPGMR SUNLinearSolver
#include "mpi.h"                       // MPI header file
#include "braid.h"                     // access to XBraid
#include "arkode/arkode_xbraid.h"      // access to ARKStep + XBraid interface


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

  // MPI variables
  MPI_Comm comm_w; // world communicator

  int nprocs_w; // total number of MPI processes in Comm world

  int myid_w; // process ID in space and time

  // Integrator settings
  realtype rtol;        // relative tolerance
  realtype atol;        // absolute tolerance
  int      order;       // ARKode method order
  bool     linear;      // enable/disable linearly implicit option
  bool     diagnostics; // output diagnostics

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
  double accesstime;

  // XBraid settings
  realtype x_tol;           // Xbraid stopping tolerance
  int      x_nt;            // number of fine grid time points
  int      x_skip;          // skip all work on first down cycle
  int      x_max_levels;    // max number of levels
  int      x_min_coarse;    // min possible coarse gird size
  int      x_nrelax;        // number of CF relaxation sweeps on all levels
  int      x_nrelax0;       // number of CF relaxation sweeps on level 0
  int      x_tnorm;         // temporal stopping norm
  int      x_cfactor;       // coarsening factor
  int      x_cfactor0;      // coarsening factor on level 0
  int      x_max_iter;      // max number of interations
  int      x_storage;       // Full storage on levels >= storage
  int      x_print_level;   // xbraid output level
  int      x_access_level;  // access level
  int      x_rfactor_limit; // refinement factor limit
  int      x_rfactor_fail;  // refinement factor on solver failure
  int      x_max_refine;    // max number of refinements
  bool     x_fmg;           // true = FMG cycle, false = V cycle
  bool     x_refine;        // enable refinement with XBraid
  bool     x_initseq;       // initialize with sequential solution
  bool     x_reltol;        // use relative tolerance
  bool     x_init_u0;       // initialize solution to initial condition
};

// -----------------------------------------------------------------------------
// Functions provided to XBraid
// -----------------------------------------------------------------------------

int MyInit(braid_App app, realtype t, braid_Vector *u_ptr);
int MyAccess(braid_App app, braid_Vector u, braid_AccessStatus astatus);

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
static int ReadInputs(int *argc, char ***argv, UserData *udata, bool outproc);

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

// Print integration statistics
static int OutputStats(void *arkode_mem, UserData *udata);

// Print integration timing
static int OutputTiming(UserData *udata);

// Check function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int flag;                   // reusable error-checking flag
  UserData *udata    = NULL;  // user data structure
  N_Vector u         = NULL;  // vector for storing solution
  SUNLinearSolver LS = NULL;  // linear solver memory structure
  void *arkode_mem   = NULL;  // ARKODE memory structure
  FILE *diagfp       = NULL;  // diagnostics output file
  braid_Core core    = NULL;  // XBraid memory structure
  braid_App app      = NULL;  // ARKode + XBraid interface structure

  // Timing variables
  chrono::time_point<chrono::steady_clock> t1;
  chrono::time_point<chrono::steady_clock> t2;

  // MPI variables
  MPI_Comm comm_w = MPI_COMM_WORLD; // MPI communicator
  int myid;                         // MPI process ID

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;

  flag = MPI_Comm_rank(comm_w, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

  // Set output process flag
  bool outproc = (myid == 0);

  // ---------------
  // Setup UserData
  // ---------------

  // Allocate and initialize user data structure with default values. The
  // defaults may be overwritten by command line inputs in ReadInputs below.
  udata = new UserData;
  flag = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) return 1;

  // Parse command line inputs
  flag = ReadInputs(&argc, &argv, udata, outproc);
  if (flag != 0) return 1;

  // Number of processes
  int nprocs_w;
  flag = MPI_Comm_size(comm_w, &nprocs_w);
  if (check_flag(&flag, "MPI_Comm_size", 1)) return 1;

  // Set communicator and number of processes in user data
  udata->comm_w   = comm_w;
  udata->nprocs_w = nprocs_w;
  udata->myid_w   = myid;

  // Output problem setup/options
  if (outproc)
  {
    flag = PrintUserData(udata);
    if (check_flag(&flag, "PrintUserData", 1)) return 1;
  }

  // Open diagnostics output file
  if (udata->diagnostics || udata->lsinfo)
  {
    stringstream fname;
    fname << "diagnostics." << setfill('0') << setw(5) << udata->myid_w
          << ".txt";

    const std::string tmp = fname.str();
    diagfp = fopen(tmp.c_str(), "w");
    if (check_flag((void *) diagfp, "fopen", 0)) return 1;
  }

  // ----------------------
  // Create serial vectors
  // ----------------------

  // Create vector for solution
  u = N_VNew_Serial(udata->nodes);
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
  int prectype = (udata->prec) ? PREC_RIGHT : PREC_NONE;

  if (udata->pcg)
  {
    LS = SUNLinSol_PCG(u, prectype, udata->liniters);
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
    LS = SUNLinSol_SPGMR(u, prectype, udata->liniters);
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
  // Setup ARKStep
  // --------------

  // Create integrator
  arkode_mem = ARKStepCreate(NULL, f, ZERO, u);
  if (check_flag((void *) arkode_mem, "ARKStepCreate", 0)) return 1;

  // Specify tolerances
  flag = ARKStepSStolerances(arkode_mem, udata->rtol, udata->atol);
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  // Attach user data
  flag = ARKStepSetUserData(arkode_mem, (void *) udata);
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

  // Attach linear solver
  flag = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

  if (udata->prec)
  {
    // Attach preconditioner
    flag = ARKStepSetPreconditioner(arkode_mem, PSetup, PSolve);
    if (check_flag(&flag, "ARKStepSetPreconditioner", 1)) return 1;

    // Set linear solver setup frequency (update preconditioner)
    flag = ARKStepSetLSetupFrequency(arkode_mem, udata->msbp);
    if (check_flag(&flag, "ARKStepSetLSetupFrequency", 1)) return 1;
  }

  // Set linear solver tolerance factor
  flag = ARKStepSetEpsLin(arkode_mem, udata->epslin);
  if (check_flag(&flag, "ARKStepSetEpsLin", 1)) return 1;

  // Select method order
  if (udata->order > 1)
  {
    // Use an ARKode provided table
    flag = ARKStepSetOrder(arkode_mem, udata->order);
    if (check_flag(&flag, "ARKStepSetOrder", 1)) return 1;
  }
  else
  {
    // Use implicit Euler (XBraid temporal refinement must be disabled)
    realtype c[1], A[1], b[1];
    ARKodeButcherTable B = NULL;

    // Create implicit Euler Butcher table
    c[0] = A[0] = b[0] = ONE;
    B = ARKodeButcherTable_Create(1, 1, 0, c, A, b, NULL);
    if (check_flag((void*) B, "ARKodeButcherTable_Create", 0)) return 1;

    // Attach the Butcher table
    flag = ARKStepSetTables(arkode_mem, 1, 0, B, NULL);
    if (check_flag(&flag, "ARKStepSetTables", 1)) return 1;

    // Free the Butcher table
    ARKodeButcherTable_Free(B);
  }

  // Specify linearly implicit non-time-dependent RHS
  if (udata->linear)
  {
    flag = ARKStepSetLinear(arkode_mem, 0);
    if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;
  }

  // Set adaptive stepping (XBraid with temporal refinement) options
  if (udata->x_refine)
  {
    // Use I controller
    flag = ARKStepSetAdaptivityMethod(arkode_mem, ARK_ADAPT_I, 1, 0, NULL);
    if (check_flag(&flag, "ARKStepSetAdaptivityMethod", 1)) return 1;

    // Set the step size reduction factor limit (1 / refinement factor limit)
    flag = ARKStepSetMinReduction(arkode_mem, ONE / udata->x_rfactor_limit);
    if (check_flag(&flag, "ARKStepSetMinReduction", 1)) return 1;

    // Set the failed solve step size reduction factor (1 / refinement factor)
    flag = ARKStepSetMaxCFailGrowth(arkode_mem, ONE / udata->x_rfactor_fail);
    if (check_flag(&flag, "ARKStepSetMaxCFailGrowth", 1)) return 1;
  }

  // Set diagnostics output file
  if (udata->diagnostics)
  {
    flag = ARKStepSetDiagnostics(arkode_mem, diagfp);
    if (check_flag(&flag, "ARKStepSetDiagnostics", 1)) return 1;
  }

  // ------------------------
  // Create XBraid interface
  // ------------------------

  // Create the ARKStep + XBraid interface
  flag = ARKBraid_Create(arkode_mem, &app);
  if (check_flag(&flag, "ARKBraid_Create", 1)) return 1;

  // Override the default initialization function
  flag = ARKBraid_SetInitFn(app, MyInit);
  if (check_flag(&flag, "ARKBraid_SetInitFn", 1)) return 1;

  // Override the default access function
  flag = ARKBraid_SetAccessFn(app, MyAccess);
  if (check_flag(&flag, "ARKBraid_SetAccesFn", 1)) return 1;

  // Initialize the ARKStep + XBraid interface
  flag = ARKBraid_BraidInit(comm_w, comm_w, ZERO, udata->tf,
                            udata->x_nt, app, &core);
  if (check_flag(&flag, "ARKBraid_BraidInit", 1)) return 1;

  // ----------------------
  // Set XBraid parameters
  // ----------------------

  flag = braid_SetTemporalNorm(core, udata->x_tnorm);
  if (check_flag(&flag, "braid_SetTemporalNorm", 1)) return 1;

  if (udata->x_reltol)
  {
    flag = braid_SetRelTol(core, udata->x_tol);
    if (check_flag(&flag, "braid_SetRelTol", 1)) return 1;
  }
  else
  {
    // Since we are using the Euclidean 2-norm in space, scale the tolerance so
    // it approximates to L2-norm.
    realtype tolfactor;
    if (udata->x_tnorm == 3)
    {
      // Infinity norm in time
      tolfactor = sqrt(udata->nx * udata->ny);
    }
    else
    {
      // 2-norm in time
      tolfactor = sqrt(udata->nx * udata->nx * udata->x_nt);
    }
    flag = braid_SetAbsTol(core, udata->x_tol * tolfactor);
    if (check_flag(&flag, "braid_SetAbsTol", 1)) return 1;
  }

  flag = braid_SetSkip(core, udata->x_skip);
  if (check_flag(&flag, "braid_SetSkip", 1)) return 1;

  flag = braid_SetMaxLevels( core, udata->x_max_levels );
  if (check_flag(&flag, "braid_SetMaxLevels", 1)) return 1;

  flag = braid_SetMinCoarse( core, udata->x_min_coarse );
  if (check_flag(&flag, "braid_SetMinCoarse", 1)) return 1;

  flag = braid_SetNRelax(core, -1, udata->x_nrelax);
  if (check_flag(&flag, "braid_SetNRelax", 1)) return 1;

  if (udata->x_nrelax0 > -1)
  {
    flag = braid_SetNRelax(core,  0, udata->x_nrelax0);
    if (check_flag(&flag, "braid_SetNRelax", 1)) return 1;
  }

  flag = braid_SetCFactor(core, -1, udata->x_cfactor);
  if (check_flag(&flag, "braid_SetCFactor", 1)) return 1;

  if (udata->x_cfactor0 > 0)
  {
    flag = braid_SetCFactor(core,  0, udata->x_cfactor0);
    if (check_flag(&flag, "braid_SetCFactor", 1)) return 1;
  }

  flag = braid_SetMaxIter(core, udata->x_max_iter);
  if (check_flag(&flag, "braid_SetMaxIter", 1)) return 1;

  if (udata->x_fmg)
  {
    // Use F-cycles
    flag = braid_SetFMG(core);
    if (check_flag(&flag, "braid_SetFMG", 1)) return 1;
  }

  flag = braid_SetPrintLevel(core, udata->x_print_level);
  if (check_flag(&flag, "braid_SetPrintLevel", 1)) return 1;

  flag = braid_SetAccessLevel(core, udata->x_access_level);
  if (check_flag(&flag, "braid_SetAccessLevel", 1)) return 1;

  if (udata->x_initseq) {
    flag =  braid_SetSeqSoln(core, 1);
    if (check_flag(&flag, "braid_SetSeqSoln", 1)) return 1;
  }

  // Temporal refinement
  if (udata->x_refine)
  {
    // Enable refinement
    flag = braid_SetRefine(core, 1);
    if (check_flag(&flag, "braid_SetRefine", 1)) return 1;

    // Set maximum number of refinements
    flag = braid_SetMaxRefinements(core, udata->x_max_refine);
    if (check_flag(&flag, "braid_SetMaxRefinements", 1)) return 1;

    // Use F-cycles
    flag = braid_SetFMG(core);
    if (check_flag(&flag, "braid_SetFMG", 1)) return 1;

    // Increase max levels after refinement
    flag = braid_SetIncrMaxLevels(core);
    if (check_flag(&flag, "braid_SetIncrMaxLevels", 1)) return 1;
  }

  // -----------------
  // "Loop" over time
  // -----------------

  // Start timer
  t1 = chrono::steady_clock::now();

  // Evolve in time
  flag = braid_Drive(core);
  if (check_flag(&flag, "braid_Drive", 1)) return 1;

  // Stop timer
  t2 = chrono::steady_clock::now();

  // Update timer
  udata->evolvetime += chrono::duration<double>(t2 - t1).count();

  // --------------
  // Final outputs
  // --------------

  // Print final integrator stats
  if (udata->output > 0)
  {
    if (outproc) cout << "Final max integrator statistics:" << endl;
    flag = OutputStats(arkode_mem, udata);
    if (check_flag(&flag, "OutputStats", 1)) return 1;
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

  if (udata->diagnostics || udata->lsinfo) fclose(diagfp);

  ARKStepFree(&arkode_mem);  // Free integrator memory
  SUNLinSolFree(LS);         // Free linear solver
  N_VDestroy(u);             // Free vectors
  FreeUserData(udata);       // Free user data
  delete udata;
  braid_Destroy(core);       // Free braid memory
  ARKBraid_Free(&app);       // Free interface memory
  flag = MPI_Finalize();     // Finalize MPI
  return 0;
}

// -----------------------------------------------------------------------------
// Functions provided to XBraid
// -----------------------------------------------------------------------------


// Create and initialize vectors
int MyInit(braid_App app, realtype t, braid_Vector *u_ptr)
{
  int      flag;
  void     *user_data;
  UserData *udata;

  // Get user data pointer
  ARKBraid_GetUserData(app, &user_data);
  udata = static_cast<UserData*>(user_data);

  // Create new vector
  N_Vector y = N_VNew_Serial(udata->nodes);
  flag = SUNBraidVector_New(y, u_ptr);
  if (flag != 0) return 1;

  // Set initial solution at all time points
  if (t == ZERO)
  {
    flag = Solution(t, y, udata);
    if (flag != 0) return 1;
  }
  else
  {
    N_VConst(ZERO, y);
  }

  return 0;
}

// Access XBraid and current vector
int MyAccess(braid_App app, braid_Vector u, braid_AccessStatus astatus)
{
  int       flag;    // return flag
  int       iter;    // current iteration number
  int       level;   // current level
  int       done;    // has XBraid finished
  realtype  t;       // current time
  void     *user_data;
  UserData *udata;

  // Timing variables
  chrono::time_point<chrono::steady_clock> t1;
  chrono::time_point<chrono::steady_clock> t2;

  // Start timer
  t1 = chrono::steady_clock::now();

  // Get user data pointer
  ARKBraid_GetUserData(app, &user_data);
  udata = static_cast<UserData*>(user_data);

  // Get current time, iteration, level, and status
  braid_AccessStatusGetTILD(astatus, &t, &iter, &level, &done);

  // Output on fine level when XBraid has finished
  if (level == 0 && done)
  {
    // Get current time index and number of fine grid points
    int index;
    int ntpts;
    braid_AccessStatusGetTIndex(astatus, &index);
    braid_AccessStatusGetNTPoints(astatus, &ntpts);

    // Extract NVector
    N_Vector y = NULL;
    flag = SUNBraidVector_GetNVector(u, &y);
    if (flag != 0) return 1;

    // Write visualization files
    if (udata->output == 2)
    {
      // Get output frequency (ensure the final time is output)
      int qout = ntpts / udata->nout;
      int rout = ntpts % udata->nout;
      int nout = (rout > 0) ? udata->nout + 2 : udata->nout + 1;

      // Output problem information
      if (index == 0)
      {
        ofstream dout;
        dout.open("heat2d_info.txt");
        dout <<  "xu  " << udata->xu << endl;
        dout <<  "yu  " << udata->yu << endl;
        dout <<  "nx  " << udata->nx << endl;
        dout <<  "ny  " << udata->ny << endl;
        dout <<  "nt  " << nout      << endl;
        dout.close();
      }

      // Output solution and error
      if (!(index % qout) || index == ntpts)
      {
        // Open output streams
        stringstream fname;
        fname << "heat2d_solution."
              << setfill('0') << setw(6) << index / qout << ".txt";

        udata->uout.open(fname.str());
        udata->uout << scientific;
        udata->uout << setprecision(numeric_limits<realtype>::digits10);

        fname.str("");
        fname.clear();
        fname << "heat2d_error."
              << setfill('0') << setw(6) << index / qout << ".txt";

        udata->eout.open(fname.str());
        udata->eout << scientific;
        udata->eout << setprecision(numeric_limits<realtype>::digits10);

        // Compute the error
        flag = SolutionError(t, y, udata->e, udata);
        if (check_flag(&flag, "SolutionError", 1)) return 1;

        // Output solution to disk
        realtype *yarray = N_VGetArrayPointer(y);
        if (check_flag((void *) yarray, "N_VGetArrayPointer", 0)) return -1;

        udata->uout << t << " ";
        for (sunindextype i = 0; i < udata->nodes; i++)
        {
          udata->uout << yarray[i] << " ";
        }
        udata->uout << endl;

        // Output error to disk
        realtype *earray = N_VGetArrayPointer(udata->e);
        if (check_flag((void *) earray, "N_VGetArrayPointer", 0)) return -1;

        udata->eout << t << " ";
        for (sunindextype i = 0; i < udata->nodes; i++)
        {
          udata->eout << earray[i] << " ";
        }
        udata->eout << endl;

        // Close output streams
        udata->uout.close();
        udata->eout.close();
      }
    }

    // Output final error
    if (index == ntpts)
    {
      // Compute the max error
      flag = SolutionError(t, y, udata->e, udata);
      if (check_flag(&flag, "SolutionError", 1)) return 1;

      realtype maxerr = N_VMaxNorm(udata->e);

      cout << scientific;
      cout << setprecision(numeric_limits<realtype>::digits10);
      cout << "  Max error = " << maxerr << endl << endl;
    }
  }

  // Stop timer
  t2 = chrono::steady_clock::now();

  // Update timing
  udata->accesstime += chrono::duration<double>(t2 - t1).count();

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

  // MPI variables
  udata->comm_w = MPI_COMM_NULL;

  udata->nprocs_w = 1;

  // Integrator settings
  udata->rtol        = RCONST(1.e-5);   // relative tolerance
  udata->atol        = RCONST(1.e-10);  // absolute tolerance
  udata->order       = 3;               // method order
  udata->linear      = true;            // linearly implicit problem
  udata->diagnostics = false;           // output diagnostics

  // Linear solver and preconditioner options
  udata->pcg       = true;       // use PCG (true) or GMRES (false)
  udata->prec      = true;       // enable preconditioning
  udata->lsinfo    = false;      // output residual history
  udata->liniters  = 100;        // max linear iterations
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
  udata->accesstime = 0.0;

  // Xbraid
  udata->x_tol           = 1.0e-6;
  udata->x_nt            = 300;
  udata->x_skip          = 1;
  udata->x_max_levels    = 15;
  udata->x_min_coarse    = 3;
  udata->x_nrelax        = 1;
  udata->x_nrelax0       = -1;
  udata->x_tnorm         = 2;
  udata->x_cfactor       = 2;
  udata->x_cfactor0      = -1;
  udata->x_max_iter      = 100;
  udata->x_storage       = -1;
  udata->x_print_level   = 1;
  udata->x_access_level  = 1;
  udata->x_rfactor_limit = 10;
  udata->x_rfactor_fail  = 4;
  udata->x_max_refine    = 8;
  udata->x_fmg           = false;
  udata->x_refine        = false;
  udata->x_initseq       = false;
  udata->x_reltol        = false;
  udata->x_init_u0       = false;

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
static int ReadInputs(int *argc, char ***argv, UserData *udata, bool outproc)
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
    else if (arg == "--order")
    {
      udata->order = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--nonlinear")
    {
      udata->linear = false;
    }
    else if (arg == "--diagnostics")
    {
      udata->diagnostics = true;
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
    // XBraid settings
    else if (arg == "--x_tol")
    {
      udata->x_tol = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--x_nt")
    {
      udata->x_nt = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_skip")
    {
      udata->x_skip = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_max_levels")
    {
      udata->x_max_levels = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_min_coarse")
    {
      udata->x_min_coarse = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_nrelax")
    {
      udata->x_nrelax = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_nrelax0")
    {
      udata->x_nrelax0 = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_tnorm")
    {
      udata->x_tnorm = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_cfactor")
    {
      udata->x_cfactor = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_cfactor0")
    {
      udata->x_cfactor0 = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_max_iter")
    {
      udata->x_max_iter = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_storage")
    {
      udata->x_storage = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_print_level")
    {
      udata->x_print_level = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_access_level")
    {
      udata->x_access_level = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_rfactor_limit")
    {
      udata->x_rfactor_limit = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_rfactor_fail")
    {
      udata->x_rfactor_fail = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_max_refine")
    {
      udata->x_max_refine = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--x_fmg")
    {
      udata->x_fmg = true;
    }
    else if (arg == "--x_refine")
    {
      udata->x_refine = true;
    }
    else if (arg == "--x_initseq")
    {
      udata->x_initseq = true;
    }
    else if (arg == "--x_reltol")
    {
      udata->x_reltol = true;
    }
    else if (arg == "--x_init_u0")
    {
      udata->x_init_u0 = true;
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
    else if (arg == "--timing")
    {
      udata->timing = true;
    }
    // Help
    else if (arg == "--help")
    {
      if (outproc) InputHelp();
      return -1;
    }
    // Unknown input
    else
    {
      if (outproc)
      {
        cerr << "ERROR: Invalid input " << arg << endl;
        InputHelp();
      }
      return -1;
    }
  }

  // Recompute total number of nodes
  udata->nodes = udata->nx * udata->ny;

  // Recompute x and y mesh spacing
  udata->dx = (udata->xu) / (udata->nx - 1);
  udata->dy = (udata->yu) / (udata->ny - 1);

  // If the method order is 1 the XBraid refinement must be disabled
  if (udata->order == 1 && !(udata->x_refine))
  {
    cerr << "ERROR: Method order 1 requires fixed time stepping" << endl;
    return -1;
  }

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
  cout << "  --nonlinear             : disable linearly implicit flag" << endl;
  cout << "  --order <ord>           : method order" << endl;
  cout << "  --diagnostics           : output diagnostics" << endl;
  cout << "  --gmres                 : use GMRES linear solver" << endl;
  cout << "  --lsinfo                : output residual history" << endl;
  cout << "  --liniters <iters>      : max number of iterations" << endl;
  cout << "  --epslin <factor>       : linear tolerance factor" << endl;
  cout << "  --noprec                : disable preconditioner" << endl;
  cout << "  --msbp <steps>          : max steps between prec setups" << endl;
  cout << "  --x_tol <tol>           : XBraid stopping tolerance" << endl;
  cout << "  --x_nt <nt>             : Initial number of time grid values" << endl;
  cout << "  --x_skip <0,1>          : Skip all work on first down cycle" << endl;
  cout << "  --x_max_levels <max>    : Max number of multigrid levels " << endl;
  cout << "  --x_min_coarse <size>   : Minimum coarse grid size" << endl;
  cout << "  --x_nrelax <num>        : Number of relaxation sweeps" << endl;
  cout << "  --x_nrelax0 <num>       : Number of relaxation sweeps on level 0" << endl;
  cout << "  --x_tnorm <1,2,3>       : Choice of temporal norm " << endl;
  cout << "  --x_cfactor <fac>       : Coarsening factor" << endl;
  cout << "  --x_cfactor0 <fac>      : Coarsening factor on level 0" << endl;
  cout << "  --x_max_iter <max>      : Max number of multigrid iterations" << endl;
  cout << "  --x_storage <lev>       : Full storage on levels >= <lev>" << endl;
  cout << "  --x_print_level <lev>   : Set print level" << endl;
  cout << "  --x_access_level <lev>  : Set access level" << endl;
  cout << "  --x_rfactor_limit <fac> : Max refinement factor" << endl;
  cout << "  --x_rfactor_fail <fac>  : Solver failure refinement factor" << endl;
  cout << "  --x_max_refine <max>    : Max number of grid refinements" << endl;
  cout << "  --x_fmg                 : Use FMG (F-cycles)" << endl;
  cout << "  --x_refine              : Enable temporal refinement" << endl;
  cout << "  --x_initseq             : Initialize with sequential solution (debug)" << endl;
  cout << "  --x_reltol              : Use relative stopping tolerance" << endl;
  cout << "  --x_init_u0             : Initialize all times with u0" << endl;
  cout << "  --output <level>        : output level" << endl;
  cout << "  --nout <nout>           : number of outputs" << endl;
  cout << "  --timing                : print timing data" << endl;
  cout << "  --help                  : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData *udata)
{
  cout << endl;
  cout << "2D Heat PDE test problem:"                     << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  nprocs         = " << udata->nprocs_w        << endl;
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
  cout << "  order          = " << udata->order           << endl;
  cout << "  linear         = " << udata->linear          << endl;
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
  cout << "  nt             = " << udata->x_nt            << endl;
  cout << "  xtol           = " << udata->x_tol           << endl;
  cout << "  refine         = " << udata->x_refine        << endl;
  cout << "  rfactor limit  = " << udata->x_rfactor_limit << endl;
  cout << "  rfactor fail   = " << udata->x_rfactor_fail  << endl;
  cout << "  init seq       = " << udata->x_initseq       << endl;
  cout << "  print level    = " << udata->x_print_level   << endl;
  cout << "  access level   = " << udata->x_access_level  << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  output         = " << udata->output          << endl;
  cout << " --------------------------------- "           << endl;
  cout << endl;

  return 0;
}

// Print integrator statistics
static int OutputStats(void *arkode_mem, UserData* udata)
{
  int flag;

  bool outproc = (udata->myid_w == 0);

  // Get integrator and solver stats
  long int nst, nst_a, netf, nfe, nfi, nni, ncfn, nli, nlcf, nsetups, nfi_ls, nJv;
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return -1;
  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(&flag, "ARKStepGetNumStepAttempts", 1)) return -1;
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(&flag, "ARKStepGetNumErrTestFails", 1)) return -1;
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return -1;
  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return -1;
  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1)) return -1;
  flag = ARKStepGetNumLinIters(arkode_mem, &nli);
  if (check_flag(&flag, "ARKStepGetNumLinIters", 1)) return -1;
  flag = ARKStepGetNumLinConvFails(arkode_mem, &nlcf);
  if (check_flag(&flag, "ARKStepGetNumLinConvFails", 1)) return -1;
  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1)) return -1;
  flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfi_ls);
  if (check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1)) return -1;
  flag = ARKStepGetNumJtimesEvals(arkode_mem, &nJv);
  if (check_flag(&flag, "ARKStepGetNumJtimesEvals", 1)) return -1;

  // Reduce stats across time
  MPI_Allreduce(MPI_IN_PLACE, &nst,     1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &nst_a,   1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &netf,    1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &nfi,     1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &nni,     1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &ncfn,    1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &nli,     1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &nlcf,    1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &nsetups, 1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &nfi_ls,  1, MPI_LONG, MPI_MAX, udata->comm_w);
  MPI_Allreduce(MPI_IN_PLACE, &nJv,     1, MPI_LONG, MPI_MAX, udata->comm_w);

  if (outproc)
  {
    cout << fixed;
    cout << setprecision(6);

    cout << "  Steps            = " << nst     << endl;
    cout << "  Step attempts    = " << nst_a   << endl;
    cout << "  Error test fails = " << netf    << endl;
    cout << "  RHS evals        = " << nfi     << endl;
    cout << "  NLS iters        = " << nni     << endl;
    cout << "  NLS fails        = " << ncfn    << endl;
    cout << "  LS iters         = " << nli     << endl;
    cout << "  LS fails         = " << nlcf    << endl;
    cout << "  LS setups        = " << nsetups << endl;
    cout << "  LS RHS evals     = " << nfi_ls  << endl;
    cout << "  Jv products      = " << nJv     << endl;
    cout << endl;

    // Compute average nls iters per step attempt and ls iters per nls iter
    realtype avgnli = (realtype) nni / (realtype) nst_a;
    realtype avgli  = (realtype) nli / (realtype) nni;
    cout << "  Avg NLS iters per step attempt = " << avgnli << endl;
    cout << "  Avg LS iters per NLS iter      = " << avgli  << endl;
    cout << endl;
  }

  // Get preconditioner stats
  if (udata->prec)
  {
    long int npe, nps;
    flag = ARKStepGetNumPrecEvals(arkode_mem, &npe);
    if (check_flag(&flag, "ARKStepGetNumPrecEvals", 1)) return -1;
    flag = ARKStepGetNumPrecSolves(arkode_mem, &nps);
    if (check_flag(&flag, "ARKStepGetNumPrecSolves", 1)) return -1;

    MPI_Allreduce(MPI_IN_PLACE, &npe, 1, MPI_LONG, MPI_MAX, udata->comm_w);
    MPI_Allreduce(MPI_IN_PLACE, &nps, 1, MPI_LONG, MPI_MAX, udata->comm_w);

    if (outproc)
    {
      cout << "  Preconditioner setups = " << npe << endl;
      cout << "  Preconditioner solves = " << nps << endl;
      cout << endl;
    }
  }

  return 0;
}

static int OutputTiming(UserData *udata)
{
  bool outproc = (udata->myid_w == 0);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double maxtime = 0.0;

  MPI_Reduce(&(udata->evolvetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_w);
  if (outproc)
  {
    cout << "  Evolve time   = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->rhstime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_w);
  if (outproc)
  {
    cout << "  RHS time      = " << maxtime << " sec" << endl;
  }

  if (udata->prec)
  {
    MPI_Reduce(&(udata->psetuptime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
               udata->comm_w);
    if (outproc)
    {
      cout << "  PSetup time   = " << maxtime << " sec" << endl;
    }

    MPI_Reduce(&(udata->psolvetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
               udata->comm_w);
    if (outproc)
    {
      cout << "  PSolve time   = " << maxtime << " sec" << endl;
      cout << endl;
    }
  }

  MPI_Reduce(&(udata->accesstime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_w);
  if (outproc)
  {
    cout << "  Access time   = " << maxtime << " sec" << endl;
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
