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
 * This example simulates the 2D diffusion-reaction (Brusselator) equation,
 *
 *   u_t = Dux u_xx + Duy u_yy + A + u * u * v - (B + 1) * u
 *   v_t = Dvx u_xx + Dvy u_yy + B * u - u * u * v
 *
 * where u and v represent the concentrations of the two chemcial species, the
 * diffusion rates are Dux = Duy = Dvx = Dvy = 1e-3, and the species with
 * constant concentration over time are A = 1 and B = 3.
 *
 * The system is evolved from t = 0 to t = 10 on a square domain centered at
 * the origin with sizes of length 1. The initial condition is
 *
 *   u(x,y) = A + 0.5 * sin(2 pi (x - xl) / wx) * sin(2 pi (y - yl) / wy)
 *   v(x,y) = B / A
 *
 * where xl and yl are the lower bounds of the domain in the x and y directions
 * respectively, wx is the width of the domain, and wy is the height of the
 * domain.
 *
 * The system is evolved in time using one of the following approaches:
 *
 *   1. A single rate IMEX method (ARKStep) with implicit diffusion and
 *      explicit reactions.
 *
 *   2. A solve-decoupled implicit MRI-GARK method (MRIStep) paired with one of
 *      the following fast time scale integrators:
 *
 *      a. An explicit method (ARKStep) integrating all the reaction systems
 *         simultaneously.
 *
 *      b. A user-defined custom inner stepper wrapping CVODE and integrating
 *         all the reaction systems simultaneously (default).
 *
 *      c. A user-defined custom inner stepper wrapping CVODE and integrating
 *         the MPI task-local reaction systems independently.
 *
 * When CVODE is used as the fast time scale integrator variable order implicit
 * Adams methods are used and the nonlinear implicit systems are solved with the
 * Anderson accelerated fixed point solver.
 *
 * Several command line options are available to change the problem parameters
 * and ARKStep/MRIStep/CVODE settings. Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <random>

// Include MPI
#include "mpi.h"

// Include desired integrators, vectors, linear solvers, and nonlinear sovlers
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_mristep.h"
#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"
#include "nvector/nvector_mpiplusx.h"
#include "sunlinsol/sunlinsol_pcg.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

// Macros for problem constants
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define NSPECIES 2

#define WIDTH (10 + numeric_limits<realtype>::digits10)

// Macro to access each species at an (x,y) location in a 1D array
#define UIDX(x,y,nx) (NSPECIES * ((nx) * (y) + (x)))
#define VIDX(x,y,nx) (NSPECIES * ((nx) * (y) + (x)) + 1)

using namespace std;


// -----------------------------------------------------------------------------
// Simple timer class
// -----------------------------------------------------------------------------


class Timer
{
public:
  Timer() : total_(0.0), start_(0.0), end_(0.0) {}
  void start()
  {
    start_ = MPI_Wtime();
  }
  void stop()
  {
    end_ = MPI_Wtime();
    total_ += (end_-start_);
  }
  realtype total() const
  {
    return total_;
  }
  realtype max(MPI_Comm comm)
  {
    double maxtime = 0.0;
    MPI_Reduce(&total_, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    return maxtime;
  }
private:
  realtype total_;
  realtype start_;
  realtype end_;
};


// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------


struct UserData
{
  // ------------------
  // Problem parameters
  // ------------------

  // Diffusion coefficients for u and v
  realtype Dux = RCONST(1.0e-3);
  realtype Duy = RCONST(1.0e-3);
  realtype Dvx = RCONST(1.0e-3);
  realtype Dvy = RCONST(1.0e-3);

  // Feed and reaction rates
  realtype A = RCONST(1.0);
  realtype B = RCONST(3.0);

  // Final simulation time
  realtype tf = RCONST(10.0);

  // Domain boundaries in x and y directions
  realtype xl = RCONST(-0.5);
  realtype xu = RCONST(0.5);
  realtype yl = RCONST(-0.5);
  realtype yu = RCONST(0.5);

  // Enable/disable RHS terms
  bool diffusion = true;
  bool reaction  = true;

  // --------------------------
  // Discretization parameteres
  // --------------------------

  // Global and local number of nodes in the x and y directions
  sunindextype nx     = 128;
  sunindextype ny     = 128;
  sunindextype nx_loc = 0;
  sunindextype ny_loc = 0;

  // Mesh spacing in the x and y directions
  realtype dx = (xu - xl) / nx;
  realtype dy = (yu - yl) / ny;

  // Global and local number of equations
  sunindextype neq     = NSPECIES * nx * ny;
  sunindextype neq_loc = 0;

  // Subdomain global starting and ending x and y indices
  sunindextype is = 0;
  sunindextype ie = 0;
  sunindextype js = 0;
  sunindextype je = 0;

  // -------------
  // MPI variables
  // -------------

  // Cartesian communicator
  MPI_Comm comm = MPI_COMM_NULL;

  // MPI processes total, in the x and y directions, and process ID
  int nprocs = 1;
  int npx    = 0;
  int npy    = 0;
  int myid   = 0;

  // Output from this process
  bool outproc = false;

  // ------------------
  // Exchange variables
  // ------------------

  // Neighbor IDs
  int ipW  = -1;
  int ipE  = -1;
  int ipS  = -1;
  int ipN  = -1;
  int ipSW = -1;
  int ipNE = -1;

  // Number of elements in buffers
  int xbufcount = 0;
  int ybufcount = 0;

  // Receive and send buffers
  realtype *Wrecv = NULL;
  realtype *Erecv = NULL;
  realtype *Srecv = NULL;
  realtype *Nrecv = NULL;

  realtype *Wsend = NULL;
  realtype *Esend = NULL;
  realtype *Ssend = NULL;
  realtype *Nsend = NULL;

  realtype *SWsend = NULL;
  realtype *NErecv = NULL;

  // Recieve and send requests
  MPI_Request reqRW, reqRE, reqRS, reqRN;
  MPI_Request reqSW, reqSE, reqSS, reqSN;
  MPI_Request reqRC, reqSC;

  // ------------------
  // Integrator options
  // ------------------

  // Flag to change integration method
  //   0 = ARKStep IMEX
  //   1 = MRIStep with ARKStep global inner integrator
  //   2 = MRIStep with CVODE global inner integrator
  //   3 = MRIStep with CVODE local inner integrator
  int integrator = 2;

  // -------------
  // IMEX settings
  // -------------

  // Relative and absolute tolerances
  realtype rtol_imex = RCONST(1.e-4);
  realtype atol_imex = RCONST(1.e-8);

  // Step size selection (ZERO = adaptive steps)
  realtype h_imex = ZERO;

  // Method order
  int order_imex = 3;

  // ------------
  // MRI settings
  // ------------

  // Relative and absolute tolerances (slow and fast)
  realtype rtol_slow = RCONST(1.e-4);
  realtype atol_slow = RCONST(1.e-8);
  realtype rtol_fast = RCONST(1.e-5);
  realtype atol_fast = RCONST(1.e-9);

  // Fixed step size (slow and fast)
  realtype h_slow = RCONST(-1.0);  // use multiple of CFL
  realtype h_fast = ZERO;          // use adaptive stepping

  // Inner ARKODE method order
  int order_fast = 3;

  // Inner stepper memory
  MRIStepInnerStepper stepper = NULL;

  // ----------------------------
  // Shared IMEX and MRI settings
  // ----------------------------

  int  controller  = 0;      // step size adaptivity method (0 = use default)
  int  maxsteps    = 0;      // max steps between outputs (0 = use default)
  bool linear      = true;   // enable/disable linearly implicit option
  bool diagnostics = false;  // output diagnostics
  FILE *diagfp     = NULL;   // diagnostics output file

  // -----------------------------------------
  // Nonlinear solver settings
  // -----------------------------------------

  int fp_iters = 10;  // max number of fixed-point iterations with CVODE
  int fp_aa    = 3;   // Anderson acceleration depth with fixed-point

  // -----------------------------------------
  // Linear solver and preconditioner settings
  // -----------------------------------------

  bool     pcg      = true;   // use PCG (true) or GMRES (false)
  bool     prec     = true;   // preconditioner on/off
  bool     lsinfo   = false;  // output residual history
  int      liniters = 5;      // number of linear iterations
  int      msbp     = 0;      // preconditioner setup frequency (0 = default)
  realtype epslin   = ZERO;   // linear solver tolerance factor (ZERO = default)
  N_Vector diag     = NULL;   // inverse of Jacobian diagonal

  // ---------------
  // Ouput variables
  // ---------------

  int output = 1;   // 0 = no output, 1 = output stats, 2 = write to disk
  int nout   = 20;  // number of output times
  ofstream uout;    // output file stream

  // ----------------
  // Timing variables
  // ----------------

  bool  timing = false;  // print timings
  Timer evolve;          // evolve time (excluding output)
  Timer rhsD;            // diffustion rhs time (including exchange)
  Timer rhsR;            // reaction rhs time
  Timer psolve;          // preconditioner solve time
  Timer exchange;        // MPI exchange time

  // ---------
  // Debugging
  // ---------

  // Run in one step mode for fixed number of steps (0 = normal mode)
  int onestep = 0;
};


// -----------------------------------------------------------------------------
// Custom inner stepper content and functions
// -----------------------------------------------------------------------------


struct InnerStepperContent
{
  void* cvode_mem = NULL;   // CVODE memory structure
  void* user_data = NULL;   // user data pointer
  bool  local     = false;  // global or task-local inner integrator

  // saved integrator stats
  long int nst  = 0;  // time steps
  long int netf = 0;  // error test fails
  long int nfe  = 0;  // rhs evals
  long int nni  = 0;  // nonlinear iterations
  long int nncf = 0;  // nonlinear convergence failures
};

static int CVodeInnerStepper_Evolve(MRIStepInnerStepper stepper,
                                    realtype t0, realtype tout, N_Vector y);
static int CVodeInnerStepper_FullRhs(MRIStepInnerStepper stepper,
                                     realtype t, N_Vector y, N_Vector f,
                                     int mode);
static int CVodeInnerStepper_Reset(MRIStepInnerStepper stepper,
                                   realtype tR, N_Vector yR);


// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------


// ODE right hand side functions
static int diffusion(realtype t, N_Vector u, N_Vector f, void *user_data);
static int reaction(realtype t, N_Vector u, N_Vector f, void *user_data);

// Preconditioner solve function
static int PSolve(realtype t, N_Vector u, N_Vector f, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data);


// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------


// Setup the parallel decomposition
static int SetupDecomp(UserData *udata);

// Integrator setup functions
static int SetupARK(SUNContext ctx, UserData *udata, N_Vector u,
                    SUNLinearSolver LS, void **arkode_mem);
static int SetupMRI(SUNContext ctx, UserData *udata, N_Vector u,
                    SUNLinearSolver LS, void **arkode_mem,
                    MRIStepInnerStepper *stepper);
static int SetupMRICVODE(SUNContext ctx, UserData *udata, N_Vector u,
                         SUNLinearSolver LS, SUNNonlinearSolver *NLS,
                         void **arkode_mem, MRIStepInnerStepper *stepper);

// Perform neighbor exchange
static int StartExchange(N_Vector y, UserData *udata);
static int EndExchange(UserData *udata);

// Exchange boundary data for output
static int ExchangeBC(N_Vector y, UserData *udata);


// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------


// Free memory allocated within UserData
static int FreeUserData(UserData *udata);

// Read the command line inputs and set UserData values
static int ReadInputs(int *argc, char ***argv, UserData *udata);


// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------


// Compute the initial condition
static int SetIC(N_Vector u, UserData *udata);

// Print the command line options
static void InputHelp();

// Print some UserData information
static int PrintUserData(UserData *udata);

// Output solution
static int OpenOutput(UserData *udata);
static int WriteOutput(realtype t, N_Vector u, UserData *udata);
static int CloseOutput(UserData *udata);

// Print integration statistics
static int OutputStatsIMEX(void *arkode_mem, UserData *udata);
static int OutputStatsMRI(void *arkode_mem, MRIStepInnerStepper stepper,
                          UserData *udata);
static int OutputStatsMRICVODE(void *arkode_mem,
                               MRIStepInnerStepper stepper,
                               UserData* udata);

// Print integration timing
static int OutputTiming(UserData *udata);

// Check function return values
static int check_flag(void *flagvalue, const string funcname, int opt);


// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------


int main(int argc, char* argv[])
{
  // Reusable error-checking flag
  int flag = 0;

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;

  // Create the SUNDIALS context object for this simulation.
  SUNContext ctx = NULL;
  MPI_Comm comm  = MPI_COMM_WORLD;
  SUNContext_Create((void*) &comm, &ctx);

  // MPI process ID
  int myid;
  flag = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

  // Set output process flag
  bool outproc = (myid == 0);

  // ---------------
  // Setup user data
  // ---------------

  UserData udata;

  udata.outproc = outproc;

  flag = ReadInputs(&argc, &argv, &udata);
  if (flag != 0) return 1;

  // ----------------------------
  // Setup parallel decomposition
  // ----------------------------

  flag = SetupDecomp(&udata);
  if (check_flag(&flag, "SetupDecomp", 1)) return 1;

  // Output problem setup/options
  if (outproc)
  {
    flag = PrintUserData(&udata);
    if (check_flag(&flag, "PrintUserData", 1)) return 1;
  }

  // --------------
  // Create vectors
  // --------------

  N_Vector u = N_VMake_MPIPlusX(udata.comm,
                                N_VNew_Serial(udata.neq_loc, ctx), ctx);
  if (check_flag((void *) u, "N_VNew_MPIPlusX", 0)) return 1;

  // --------------------
  // Create linear solver
  // --------------------

  // Open diagnostics file
  if (outproc && (udata.diagnostics || udata.lsinfo))
  {
    udata.diagfp = fopen("diagnostics.txt", "w");
    if (check_flag((void *) (udata.diagfp), "fopen", 0)) return 1;
  }

  // Preconditioning type
  int prectype = (udata.prec) ? SUN_PREC_RIGHT : SUN_PREC_NONE;

  // Linear solver memory structure
  SUNLinearSolver LS = NULL;

  if (udata.pcg)
  {
    LS = SUNLinSol_PCG(u, prectype, udata.liniters, ctx);
    if (check_flag((void *) LS, "SUNLinSol_PCG", 0)) return 1;

    if (udata.lsinfo && outproc)
    {
      flag = SUNLinSolSetPrintLevel_PCG(LS, 1);
      if (check_flag(&flag, "SUNLinSolSetPrintLevel_PCG", 1)) return(1);

      flag = SUNLinSolSetInfoFile_PCG(LS, udata.diagfp);
      if (check_flag(&flag, "SUNLinSolSetInfoFile_PCG", 1)) return(1);
    }
  }
  else
  {
    LS = SUNLinSol_SPGMR(u, prectype, udata.liniters, ctx);
    if (check_flag((void *) LS, "SUNLinSol_SPGMR", 0)) return 1;

    if (udata.lsinfo && outproc)
    {
      flag = SUNLinSolSetPrintLevel_SPGMR(LS, 1);
      if (check_flag(&flag, "SUNLinSolSetPrintLevel_SPGMR", 1)) return(1);

      flag = SUNLinSolSetInfoFile_SPGMR(LS, udata.diagfp);
      if (check_flag(&flag, "SUNLinSolSetInfoFile_SPGMR", 1)) return(1);
    }
  }

  // Allocate preconditioner workspace
  if (udata.prec)
  {
    udata.diag = N_VClone(u);
    if (check_flag((void *) (udata.diag), "N_VClone", 0)) return 1;
  }

  // ---------------------
  // Set initial condition
  // ---------------------

  flag = SetIC(u, &udata);
  if (check_flag(&flag, "SetIC", 1)) return 1;

  // ----------------
  // Setup Integrator
  // ----------------

  // ARKODE memory structure
  void *arkode_mem = NULL;

  // Custom inner stepper memory (CVODE)
  MRIStepInnerStepper stepper = NULL;

  // Inner stepper nonlinear solver (CVODE)
  SUNNonlinearSolver NLS = NULL;

  // Create integrator
  switch(udata.integrator)
  {
  case(0):
    flag = SetupARK(ctx, &udata, u, LS, &arkode_mem);
    if (check_flag((void *) arkode_mem, "SetupARK", 0)) return 1;
    break;
  case(1):
    flag = SetupMRI(ctx, &udata, u, LS, &arkode_mem, &stepper);
    if (check_flag((void *) arkode_mem, "SetupMRI", 0)) return 1;
    break;
  case(2):
  case(3):
    flag = SetupMRICVODE(ctx, &udata, u, LS, &NLS, &arkode_mem, &stepper);
    if (check_flag((void *) arkode_mem, "SetupMRICVODE", 0)) return 1;
    break;
  default:
    cerr << "Invalid integrator option" << endl;
    break;
  }

  // ----------------------
  // Evolve problem in time
  // ----------------------

  // Set the step mode
  int stepmode = ARK_NORMAL;

  if (udata.onestep)
  {
    udata.nout = udata.onestep;
    stepmode   = ARK_ONE_STEP;
  }

  // Initial time, time between outputs, output time
  realtype t     = ZERO;
  realtype dTout = udata.tf / udata.nout;
  realtype tout  = dTout;

  // Inital output
  flag = OpenOutput(&udata);
  if (check_flag(&flag, "OpenOutput", 1)) return 1;

  flag = WriteOutput(t, u, &udata);
  if (check_flag(&flag, "WriteOutput", 1)) return 1;

  // Loop over output times
  for (int iout = 0; iout < udata.nout; iout++)
  {
    // Start timer
    udata.evolve.start();

    // Evolve
    if (udata.integrator)
    {
      flag = MRIStepEvolve(arkode_mem, tout, u, &t, stepmode);
      if (check_flag(&flag, "MRIStepEvolve", 1)) break;
    }
    else
    {
      flag = ARKStepEvolve(arkode_mem, tout, u, &t, stepmode);
      if (check_flag(&flag, "ARKStepEvolve", 1)) break;
    }

    // Stop timer
    udata.evolve.stop();

    // Output solution
    flag = WriteOutput(t, u, &udata);
    if (check_flag(&flag, "WriteOutput", 1)) return 1;

    // Update output time
    tout += dTout;
    tout = (tout > udata.tf) ? udata.tf : tout;
  }

  // Close output
  flag = CloseOutput(&udata);
  if (check_flag(&flag, "CloseOutput", 1)) return 1;

  // -------------
  // Final outputs
  // -------------

  // Print final integrator stats
  if (udata.output > 0 && outproc)
  {
    cout << "Final integrator statistics:" << endl;
    switch(udata.integrator)
    {
    case(0):
      flag = OutputStatsIMEX(arkode_mem, &udata);
      if (check_flag(&flag, "OutputStatsIMEX", 1)) return 1;
      break;
    case(1):
      flag = OutputStatsMRI(arkode_mem, stepper, &udata);
      if (check_flag(&flag, "OutputStatsMRI", 1)) return 1;
      break;
    case(2):
    case(3):
      flag = OutputStatsMRICVODE(arkode_mem, stepper, &udata);
      if (check_flag(&flag, "OutputStatsMRICVODE", 1)) return 1;
      break;
    default:
      cerr << "Invalid integrator option" << endl;
      break;
    }
  }

  // Print timing
  if (udata.timing)
  {
    flag = OutputTiming(&udata);
    if (check_flag(&flag, "OutputTiming", 1)) return 1;
  }

  // --------------------
  // Clean up and return
  // --------------------

  if (outproc && (udata.diagnostics || udata.lsinfo)) fclose(udata.diagfp);

  switch(udata.integrator)
  {
  case(0):
    ARKStepFree(&arkode_mem);
    break;
  case(1):
    {
      void* inner_arkode_mem = NULL;
      MRIStepInnerStepper_GetContent(stepper, &inner_arkode_mem);
      ARKStepFree(&inner_arkode_mem);
      MRIStepInnerStepper_Free(&stepper);
      MRIStepFree(&arkode_mem);
      break;
    }
  case(2):
  case(3):
    {
      void* inner_content = NULL;
      MRIStepInnerStepper_GetContent(stepper, &inner_content);
      InnerStepperContent* content = (InnerStepperContent *) inner_content;
      CVodeFree(&(content->cvode_mem));
      delete content;
      MRIStepInnerStepper_Free(&stepper);
      SUNNonlinSolFree(NLS);
      MRIStepFree(&arkode_mem);
      break;
    }
  default:
    cerr << "Invalid integrator option" << endl;
    break;
  }

  SUNLinSolFree(LS);
  N_VDestroy(N_VGetLocalVector_MPIPlusX(u));
  N_VDestroy(u);
  FreeUserData(&udata);
  SUNContext_Free(&ctx);
  flag = MPI_Finalize();
  return 0;
}


// -----------------------------------------------------------------------------
// Setup the parallel decomposition
// -----------------------------------------------------------------------------


static int SetupDecomp(UserData *udata)
{
  int flag;

  // Check that this has not been called before
  if (udata->Erecv != NULL || udata->Wrecv != NULL ||
      udata->Srecv != NULL || udata->Nrecv != NULL)
  {
    cerr << "SetupDecomp error: parallel decomposition already set up" << endl;
    return -1;
  }

  // Get the number of processes
  flag = MPI_Comm_size(MPI_COMM_WORLD, &(udata->nprocs));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << flag << endl;
    return -1;
  }

  // Set up 2D Cartesian communicator
  int dims[2];
  dims[0] = (udata->npx > 0) ? udata->npx : 0;
  dims[1] = (udata->npy > 0) ? udata->npy : 0;

  int periods[2];
  periods[0] = 1;
  periods[1] = 1;

  flag = MPI_Dims_create(udata->nprocs, 2, dims);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Dims_create = " << flag << endl;
    return -1;
  }

  udata->npx = dims[0];
  udata->npy = dims[1];

  flag = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &(udata->comm));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_create = " << flag << endl;
    return -1;
  }

  // Get my rank in the new Cartesian communicator
  flag = MPI_Comm_rank(udata->comm, &(udata->myid));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << flag << endl;
    return -1;
  }

  if (udata->myid == 0) udata->outproc = true;

  // Get dimension of the Cartesian communicator and my coordinates
  int coords[2];
  flag = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_get = " << flag << endl;
    return -1;
  }

  // Determine local extents in x-direction
  int idx         = coords[0];
  sunindextype qx = udata->nx / dims[0];
  sunindextype rx = udata->nx % dims[0];

  udata->is = qx * idx + (idx < rx ? idx : rx);
  udata->ie = udata->is + qx - 1 + (idx < rx ? 1 : 0);

  // Sanity check
  if (udata->ie > (udata->nx - 1))
  {
    cerr << "Error ie > nx - 1" << endl;
    return -1;
  }

  // Determine local extents in y-direction
  int idy         = coords[1];
  sunindextype qy = udata->ny / dims[1];
  sunindextype ry = udata->ny % dims[1];

  udata->js = qy * idy + (idy < ry ? idy : ry);
  udata->je = udata->js + qy - 1 + (idy < ry ? 1 : 0);

  // Sanity check
  if (udata->je > (udata->ny - 1))
  {
    cerr << "Error je > ny - 1" << endl;
    return -1;
  }

  // Number of local nodes
  udata->nx_loc = (udata->ie) - (udata->is) + 1;
  udata->ny_loc = (udata->je) - (udata->js) + 1;

  // Initialize global and local vector lengths
  udata->neq     = NSPECIES * udata->nx * udata->ny;
  udata->neq_loc = NSPECIES * udata->nx_loc * udata->ny_loc;

  // Allocate exchange buffers if necessary
  udata->ybufcount = NSPECIES * udata->ny_loc;
  udata->Wrecv     = new realtype[udata->ybufcount];
  udata->Wsend     = new realtype[udata->ybufcount];
  udata->Erecv     = new realtype[udata->ybufcount];
  udata->Esend     = new realtype[udata->ybufcount];

  udata->xbufcount = NSPECIES * udata->nx_loc;
  udata->Srecv     = new realtype[udata->xbufcount];
  udata->Ssend     = new realtype[udata->xbufcount];
  udata->Nrecv     = new realtype[udata->xbufcount];
  udata->Nsend     = new realtype[udata->xbufcount];

  udata->SWsend = new realtype[NSPECIES];
  udata->NErecv = new realtype[NSPECIES];

  // MPI neighborhood information
  int nbcoords[2];

  // West neighbor
  nbcoords[0] = coords[0]-1;
  nbcoords[1] = coords[1];
  flag = MPI_Cart_rank(udata->comm, nbcoords, &(udata->ipW));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_rank = " << flag << endl;
    return -1;
  }

  // East neighbor
  nbcoords[0] = coords[0]+1;
  nbcoords[1] = coords[1];
  flag = MPI_Cart_rank(udata->comm, nbcoords, &(udata->ipE));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_rank = " << flag << endl;
    return -1;
  }

  // South neighbor
  nbcoords[0] = coords[0];
  nbcoords[1] = coords[1]-1;
  flag = MPI_Cart_rank(udata->comm, nbcoords, &(udata->ipS));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_rank = " << flag << endl;
    return -1;
  }

  // North neighbor
  nbcoords[0] = coords[0];
  nbcoords[1] = coords[1]+1;
  flag = MPI_Cart_rank(udata->comm, nbcoords, &(udata->ipN));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_rank = " << flag << endl;
    return -1;
  }

  // Opposite Corners for periodic BC output
  if (udata->is == 0 && udata->js == 0)
  {
    nbcoords[0] = coords[0] - 1;
    nbcoords[1] = coords[1] - 1;
    flag = MPI_Cart_rank(udata->comm, nbcoords, &(udata->ipSW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  if (udata->ie == udata->nx - 1 && udata->je == udata->ny - 1)
  {
    nbcoords[0] = coords[0] + 1;
    nbcoords[1] = coords[1] + 1;
    flag = MPI_Cart_rank(udata->comm, nbcoords, &(udata->ipNE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // Return success
  return 0;
}


// -----------------------------------------------------------------------------
// Setup the integrator
// -----------------------------------------------------------------------------


static int SetupARK(SUNContext ctx, UserData* udata, N_Vector u,
                    SUNLinearSolver LS, void** arkode_mem)
{
  int flag;

  // Optionally enable/disable diffusion or reactions (helpful for debugging)
  ARKRhsFn fe = udata->reaction  ? reaction  : NULL;
  ARKRhsFn fi = udata->diffusion ? diffusion : NULL;

  // Create ARKStep memory with explicit reactions and implicit diffusion
  *arkode_mem = ARKStepCreate(fe, fi, ZERO, u, ctx);
  if (check_flag((void *) *arkode_mem, "ARKStepCreate", 0)) return 1;

  // Specify tolerances
  flag = ARKStepSStolerances(*arkode_mem, udata->rtol_imex, udata->atol_imex);
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  // Attach user data
  flag = ARKStepSetUserData(*arkode_mem, (void *) udata);
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

  if (udata->diffusion)
  {
    // Attach linear solver
    flag = ARKStepSetLinearSolver(*arkode_mem, LS, NULL);
    if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

    if (udata->prec)
    {
      // Attach preconditioner
      flag = ARKStepSetPreconditioner(*arkode_mem, NULL, PSolve);
      if (check_flag(&flag, "ARKStepSetPreconditioner", 1)) return 1;

      // Set linear solver setup frequency (update preconditioner)
      flag = ARKStepSetLSetupFrequency(*arkode_mem, udata->msbp);
      if (check_flag(&flag, "ARKStepSetLSetupFrequency", 1)) return 1;
    }

    // Set linear solver tolerance factor
    flag = ARKStepSetEpsLin(*arkode_mem, udata->epslin);
    if (check_flag(&flag, "ARKStepSetEpsLin", 1)) return 1;

    // Specify linearly implicit non-time-dependent RHS
    if (udata->linear)
    {
      flag = ARKStepSetLinear(*arkode_mem, 0);
      if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;
    }
  }

  // Select method order
  flag = ARKStepSetOrder(*arkode_mem, udata->order_imex);
  if (check_flag(&flag, "ARKStepSetOrder", 1)) return 1;

  // Set fixed step size or adaptivity method
  if (udata->h_imex > ZERO)
  {
    flag = ARKStepSetFixedStep(*arkode_mem, udata->h_imex);
    if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;
  }
  else
  {
    flag = ARKStepSetAdaptivityMethod(*arkode_mem, udata->controller, SUNTRUE,
                                      SUNFALSE, NULL);
    if (check_flag(&flag, "ARKStepSetAdaptivityMethod", 1)) return 1;
  }

  // Set max steps between outputs
  flag = ARKStepSetMaxNumSteps(*arkode_mem, udata->maxsteps);
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;

  // Set stopping time
  flag = ARKStepSetStopTime(*arkode_mem, udata->tf);
  if (check_flag(&flag, "ARKStepSetStopTime", 1)) return 1;

  // Set diagnostics output file
  if (udata->diagnostics && udata->outproc)
  {
    flag = ARKStepSetDiagnostics(*arkode_mem, udata->diagfp);
    if (check_flag(&flag, "ARKStepSetDiagnostics", 1)) return 1;
  }

  return 0;
}


static int SetupMRI(SUNContext ctx, UserData* udata, N_Vector y,
                    SUNLinearSolver LS, void** arkode_mem,
                    MRIStepInnerStepper* stepper)
{
  int flag;

  // -------------------------
  // Setup the fast integrator
  // -------------------------

  // Create fast explicit integrator for reactions
  void *inner_arkode_mem = ARKStepCreate(reaction, NULL, ZERO, y, ctx);
  if (check_flag((void *) inner_arkode_mem, "ARKStepCreate", 0)) return 1;

  // Specify tolerances
  flag = ARKStepSStolerances(inner_arkode_mem, udata->rtol_fast,
                             udata->atol_fast);
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  // Attach user data
  flag = ARKStepSetUserData(inner_arkode_mem, (void *) udata);
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

  // Select method order
  flag = ARKStepSetOrder(inner_arkode_mem, udata->order_fast);
  if (check_flag(&flag, "ARKStepSetOrder", 1)) return 1;

  // Set fixed step size or adaptivity method
  if (udata->h_fast > ZERO)
  {
    flag = ARKStepSetFixedStep(inner_arkode_mem, udata->h_fast);
    if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;
  }
  else
  {
    flag = ARKStepSetAdaptivityMethod(inner_arkode_mem, udata->controller,
                                      SUNTRUE, SUNFALSE, NULL);
    if (check_flag(&flag, "ARKStepSetAdaptivityMethod", 1)) return 1;
  }

  // Set max steps between outputs
  flag = ARKStepSetMaxNumSteps(inner_arkode_mem, udata->maxsteps);
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;

  // Set diagnostics output file
  if (udata->diagnostics && udata->outproc)
  {
    flag = ARKStepSetDiagnostics(inner_arkode_mem, udata->diagfp);
    if (check_flag(&flag, "ARKStepSetDiagnostics", 1)) return 1;
  }

  // Wrap ARKODE as an MRIStepInnerStepper
  flag = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem,
                                          stepper);
  if (check_flag(&flag, "ARKStepCreateMRIStepInnerStepper", 1)) return 1;

  // -------------------------
  // Setup the slow integrator
  // -------------------------

  // Create slow integrator for diffusion and attach fast integrator
  *arkode_mem = MRIStepCreate(NULL, diffusion, ZERO, y, *stepper, ctx);
  if (check_flag((void *)*arkode_mem, "MRIStepCreate", 0)) return 1;

  // Set method coupling table (solve-decoupled implicit method)
  MRIStepCoupling C = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ESDIRK34a);
  if (check_flag((void*)C, "MRIStepCoupling_LoadTable", 1)) return 1;

  flag = MRIStepSetCoupling(*arkode_mem, C);
  if (check_flag(&flag, "MRIStepSetCoupling", 1)) return 1;

  MRIStepCoupling_Free(C);

  // Set the slow step size
  flag = MRIStepSetFixedStep(*arkode_mem, udata->h_slow);
  if (check_flag(&flag, "MRIStepSetFixedStep", 1)) return 1;

  // Specify tolerances
  flag = MRIStepSStolerances(*arkode_mem, udata->rtol_slow, udata->atol_slow);
  if (check_flag(&flag, "MRIStepSStolerances", 1)) return 1;

  // Attach user data
  flag = MRIStepSetUserData(*arkode_mem, (void *) udata);
  if (check_flag(&flag, "MRIStepSetUserData", 1)) return 1;

  // Attach linear solver
  flag = MRIStepSetLinearSolver(*arkode_mem, LS, NULL);
  if (check_flag(&flag, "MRIStepSetLinearSolver", 1)) return 1;

  if (udata->prec)
  {
    // Attach preconditioner
    flag = MRIStepSetPreconditioner(*arkode_mem, NULL, PSolve);
    if (check_flag(&flag, "MRIStepSetPreconditioner", 1)) return 1;

    // Set linear solver setup frequency (update preconditioner)
    flag = MRIStepSetLSetupFrequency(*arkode_mem, udata->msbp);
    if (check_flag(&flag, "MRIStepSetLSetupFrequency", 1)) return 1;
  }

  // Set linear solver tolerance factor
  flag = MRIStepSetEpsLin(*arkode_mem, udata->epslin);
  if (check_flag(&flag, "MRIStepSetEpsLin", 1)) return 1;

  // Specify linearly implicit non-time-dependent RHS
  if (udata->linear)
  {
    flag = MRIStepSetLinear(*arkode_mem, 0);
    if (check_flag(&flag, "MRIStepSetLinear", 1)) return 1;
  }

  // Set max steps between outputs
  flag = MRIStepSetMaxNumSteps(*arkode_mem, udata->maxsteps);
  if (check_flag(&flag, "MRIStepSetMaxNumSteps", 1)) return 1;

  // Set stopping time
  flag = MRIStepSetStopTime(*arkode_mem, udata->tf);
  if (check_flag(&flag, "MRIStepSetStopTime", 1)) return 1;

  // Set diagnostics output file
  if (udata->diagnostics && udata->outproc)
  {
    flag = MRIStepSetDiagnostics(*arkode_mem, udata->diagfp);
    if (check_flag(&flag, "MRIStepSetDiagnostics", 1)) return 1;
  }

  return 0;
}


static int SetupMRICVODE(SUNContext ctx, UserData *udata, N_Vector y,
                         SUNLinearSolver LS, SUNNonlinearSolver *NLS,
                         void **arkode_mem, MRIStepInnerStepper* stepper)
{
  int flag;

  // -------------------------
  // Setup the fast integrator
  // -------------------------

  // Use the global or local state vector to create the inner integrator
  N_Vector y_vec;
  if (udata->integrator == 2)
  {
    y_vec = y;
  }
  else if (udata->integrator == 3)
  {
    y_vec = N_VGetLocalVector_MPIPlusX(y);
  }
  else
  {
    cerr << "ERROR: Invalid MRIStep + CVODE option" << endl;
    return -1;
  }

  // Create the solver memory and specify the Adams methods
  void *cvode_mem = CVodeCreate(CV_ADAMS, ctx);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  // Initialize the integrator memory
  flag = CVodeInit(cvode_mem, reaction, ZERO, y_vec);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  // Specify tolerances
  flag = CVodeSStolerances(cvode_mem, udata->rtol_fast, udata->atol_fast);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, (void *) udata);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return 1;

  // Create and attach fixed-point nonlinear solver
  *NLS = SUNNonlinSol_FixedPoint(y_vec, udata->fp_aa, ctx);
  if(check_flag((void *)*NLS, "SUNNonlinSol_FixedPoint", 0)) return(1);

  flag = CVodeSetNonlinearSolver(cvode_mem, *NLS);
  if(check_flag(&flag, "CVodeSetNonlinearSolver", 1)) return(1);

  // Set max number of fixed-point iterations
  flag = CVodeSetMaxNonlinIters(cvode_mem, udata->fp_iters);
  if(check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

  // Set max steps between outputs
  flag = CVodeSetMaxNumSteps(cvode_mem, udata->maxsteps);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return 1;

  // Create the inner stepper wrapper
  flag = MRIStepInnerStepper_Create(ctx, stepper);
  if (check_flag(&flag, "MRIStepInnerStepper_Create", 1)) return 1;

  // Attach memory and operations
  InnerStepperContent *inner_content = new InnerStepperContent;

  inner_content->cvode_mem = cvode_mem;
  inner_content->user_data = udata;
  inner_content->local = (udata->integrator == 2) ? false : true;

  flag = MRIStepInnerStepper_SetContent(*stepper, inner_content);
  if (check_flag(&flag, "MRIStepInnerStepper_SetContent", 1)) return 1;

  flag = MRIStepInnerStepper_SetEvolveFn(*stepper,
                                         CVodeInnerStepper_Evolve);
  if (check_flag(&flag, "MRIStepInnerStepper_SetEvolve", 1)) return 1;

  flag = MRIStepInnerStepper_SetFullRhsFn(*stepper,
                                          CVodeInnerStepper_FullRhs);
  if (check_flag(&flag, "MRIStepInnerStepper_SetFullRhsFn", 1)) return 1;

  flag = MRIStepInnerStepper_SetResetFn(*stepper,
                                        CVodeInnerStepper_Reset);
  if (check_flag(&flag, "MRIStepInnerStepper_SetResetFn", 1)) return 1;

  // Attach inner stepper memory to user data
  udata->stepper = *stepper;

  // -------------------------
  // Setup the slow integrator
  // -------------------------

  // Create slow integrator for diffusion and attach fast integrator
  *arkode_mem = MRIStepCreate(NULL, diffusion, ZERO, y, *stepper, ctx);
  if (check_flag((void *)*arkode_mem, "MRIStepCreate", 0)) return 1;

  // Set method coupling table (solve-decoupled implicit method)
  MRIStepCoupling C = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ESDIRK34a);
  if (check_flag((void*)C, "MRIStepCoupling_LoadTable", 1)) return 1;

  flag = MRIStepSetCoupling(*arkode_mem, C);
  if (check_flag(&flag, "MRIStepSetCoupling", 1)) return 1;

  MRIStepCoupling_Free(C);

  // Set the slow step size
  flag = MRIStepSetFixedStep(*arkode_mem, udata->h_slow);
  if (check_flag(&flag, "MRIStepSetFixedStep", 1)) return 1;

  // Specify tolerances
  flag = MRIStepSStolerances(*arkode_mem, udata->rtol_slow, udata->atol_slow);
  if (check_flag(&flag, "MRIStepSStolerances", 1)) return 1;

  // Attach user data
  flag = MRIStepSetUserData(*arkode_mem, (void *) udata);
  if (check_flag(&flag, "MRIStepSetUserData", 1)) return 1;

  // Attach linear solver
  flag = MRIStepSetLinearSolver(*arkode_mem, LS, NULL);
  if (check_flag(&flag, "MRIStepSetLinearSolver", 1)) return 1;

  if (udata->prec)
  {
    // Attach preconditioner
    flag = MRIStepSetPreconditioner(*arkode_mem, NULL, PSolve);
    if (check_flag(&flag, "MRIStepSetPreconditioner", 1)) return 1;

    // Set linear solver setup frequency (update preconditioner)
    flag = MRIStepSetLSetupFrequency(*arkode_mem, udata->msbp);
    if (check_flag(&flag, "MRIStepSetLSetupFrequency", 1)) return 1;
  }

  // Set linear solver tolerance factor
  flag = MRIStepSetEpsLin(*arkode_mem, udata->epslin);
  if (check_flag(&flag, "MRIStepSetEpsLin", 1)) return 1;

  // Specify linearly implicit non-time-dependent RHS
  if (udata->linear)
  {
    flag = MRIStepSetLinear(*arkode_mem, 0);
    if (check_flag(&flag, "MRIStepSetLinear", 1)) return 1;
  }

  // Set max steps between outputs
  flag = MRIStepSetMaxNumSteps(*arkode_mem, udata->maxsteps);
  if (check_flag(&flag, "MRIStepSetMaxNumSteps", 1)) return 1;

  // Set stopping time
  flag = MRIStepSetStopTime(*arkode_mem, udata->tf);
  if (check_flag(&flag, "MRIStepSetStopTime", 1)) return 1;

  // Set diagnostics output file
  if (udata->diagnostics && udata->outproc)
  {
    flag = MRIStepSetDiagnostics(*arkode_mem, udata->diagfp);
    if (check_flag(&flag, "MRIStepSetDiagnostics", 1)) return 1;
  }

  return 0;
}


// -----------------------------------------------------------------------------
// Custom inner stepper functions
// -----------------------------------------------------------------------------


static int CVodeInnerStepper_Evolve(MRIStepInnerStepper stepper,
                                    realtype t0, realtype tout, N_Vector y)
{
  int      flag;
  realtype tret;
  void*    inner_content = NULL;

  flag = MRIStepInnerStepper_GetContent(stepper, &inner_content);
  if (check_flag(&flag, "MRIStepInnerStepper_GetContent", 1)) return -1;

  InnerStepperContent *content = (InnerStepperContent*) inner_content;

  N_Vector y_vec;
  if (content->local)
  {
    // Using local inner stepper, extract the local serial vector
    y_vec = N_VGetLocalVector_MPIPlusX(y);
  }
  else
  {
    // Using global inner stepper, use the MPIPlusX vector
    y_vec = y;
  }

  flag = CVodeSetStopTime(content->cvode_mem, tout);
  if (check_flag(&flag, "CVodeSetStopTime", 1)) return -1;

  flag = CVode(content->cvode_mem, tout, y_vec, &tret, CV_NORMAL);
  if (flag < 0) return -1;

  return 0;
}


static int CVodeInnerStepper_FullRhs(MRIStepInnerStepper stepper,
                                     realtype t, N_Vector y, N_Vector f,
                                     int mode)
{
  int   flag;
  void* inner_content = NULL;

  flag = MRIStepInnerStepper_GetContent(stepper, &inner_content);
  if (check_flag(&flag, "MRIStepInnerStepper_GetContent", 1)) return -1;

  InnerStepperContent *content = (InnerStepperContent*) inner_content;
  UserData *udata = (UserData *) content->user_data;

  // disable forcing
  int integrator = udata->integrator;
  udata->integrator = 0;

  flag = reaction(t, y, f, content->user_data);
  if (flag != 0) return -1;

  // enable forcing
  udata->integrator = integrator;

  return 0;
}


static int CVodeInnerStepper_Reset(MRIStepInnerStepper stepper,
                                   realtype tR, N_Vector yR)
{
  int   flag;
  void* inner_content = NULL;

  flag = MRIStepInnerStepper_GetContent(stepper, &inner_content);
  if (check_flag(&flag, "MRIStepInnerStepper_GetContent", 1)) return -1;

  InnerStepperContent *content = (InnerStepperContent*) inner_content;

  N_Vector yR_vec;
  if (content->local)
  {
    yR_vec = N_VGetLocalVector_MPIPlusX(yR);
  }
  else
  {
    yR_vec = yR;
  }

  // Save current stats before reinitializing
  long int nst, netf, nfe, nni, nncf;

  flag = CVodeGetNumSteps(content->cvode_mem, &nst);
  if (check_flag(&flag, "CVodeGetNumSteps", 1)) return -1;
  content->nst += nst;

  flag = CVodeGetNumErrTestFails(content->cvode_mem, &netf);
  if (check_flag(&flag, "CVodeGetNumErrTestFails", 1)) return -1;
  content->netf += netf;

  flag = CVodeGetNumRhsEvals(content->cvode_mem, &nfe);
  if (check_flag(&flag, "CVodeGetNumRhsEvals", 1)) return -1;
  content->nfe += nfe;

  flag = CVodeGetNumNonlinSolvIters(content->cvode_mem, &nni);
  if (check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1)) return -1;
  content->nni += nni;

  flag = CVodeGetNumNonlinSolvConvFails(content->cvode_mem, &nncf);
  if (check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1)) return -1;
  content->nncf += nncf;

  // Reinitialize CVODE with new state
  flag = CVodeReInit(content->cvode_mem, tR, yR_vec);
  if (check_flag(&flag, "CVodeReInit", 1)) return -1;

  return 0;
}


// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------


// Routine to compute the ODE diffusion RHS function
static int diffusion(realtype t, N_Vector y, N_Vector f, void *user_data)
{
  int flag;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Start timer
  udata->rhsD.start();

  // Open exchange receives and exchange data
  flag = StartExchange(y, udata);
  if (check_flag(&flag, "StartExchange", 1)) return -1;

  // Constants for computing diffusion term
  realtype cxu = udata->Dux / (udata->dx * udata->dx);
  realtype cyu = udata->Duy / (udata->dy * udata->dy);
  realtype ccu = -TWO * (cxu + cyu);

  realtype cxv = udata->Dvx / (udata->dx * udata->dx);
  realtype cyv = udata->Dvy / (udata->dy * udata->dy);
  realtype ccv = -TWO * (cxv + cyv);

  // Access data arrays
  realtype *ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) ydata, "N_VGetArrayPointer", 0)) return -1;

  realtype *fdata = N_VGetArrayPointer(f);
  if (check_flag((void *) fdata, "N_VGetArrayPointer", 0)) return -1;

  // Shortcuts to array indices (center, west, east, south, north)
  sunindextype uc, uw, ue, us, un;
  sunindextype vc, vw, ve, vs, vn;

  // Shortcuts to local number of nodes
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Compute diffusion term on subdomain
  for (sunindextype j = 1; j < ny_loc - 1; j++)
  {
    for (sunindextype i = 1; i < nx_loc - 1; i++)
    {
      uc = UIDX(i,   j,   nx_loc);
      uw = UIDX(i-1, j,   nx_loc);
      ue = UIDX(i+1, j,   nx_loc);
      us = UIDX(i,   j-1, nx_loc);
      un = UIDX(i,   j+1, nx_loc);

      vc = VIDX(i,   j,   nx_loc);
      vw = VIDX(i-1, j,   nx_loc);
      ve = VIDX(i+1, j,   nx_loc);
      vs = VIDX(i,   j-1, nx_loc);
      vn = VIDX(i,   j+1, nx_loc);

      fdata[uc] = ccu * ydata[uc]
        + cxu * (ydata[uw] + ydata[ue])
        + cyu * (ydata[us] + ydata[un]);

      fdata[vc] = ccv * ydata[vc]
        + cxv * (ydata[vw] + ydata[ve])
        + cyv * (ydata[vs] + ydata[vn]);
    }
  }

  // Wait for exchange receives and compute diffusion term on subdomain boundary
  flag = EndExchange(udata);
  if (check_flag(&flag, "EndExchange", 1)) return -1;

  realtype *Wdata = udata->Wrecv;
  realtype *Edata = udata->Erecv;
  realtype *Sdata = udata->Srecv;
  realtype *Ndata = udata->Nrecv;

  sunindextype i, j;

  // -----------------------------------------------------
  // West face (updates south-west and north-west corners)
  // -----------------------------------------------------
  i = 0;

  // South-West corner
  j = 0;

  uc = UIDX(i,   j,   nx_loc);
  ue = UIDX(i+1, j,   nx_loc);
  un = UIDX(i,   j+1, nx_loc);

  vc = VIDX(i,   j,   nx_loc);
  ve = VIDX(i+1, j,   nx_loc);
  vn = VIDX(i,   j+1, nx_loc);

  fdata[uc] = ccu * ydata[uc]
    + cxu * (Wdata[NSPECIES * j] + ydata[ue])
    + cyu * (Sdata[NSPECIES * i] + ydata[un]);

  fdata[vc] = ccv * ydata[vc]
    + cxv * (Wdata[NSPECIES * j + 1] + ydata[ve])
    + cyv * (Sdata[NSPECIES * i + 1] + ydata[vn]);

  // West face interior
  for (j = 1; j < ny_loc - 1; j++)
  {
    uc = UIDX(i,   j,   nx_loc);
    ue = UIDX(i+1, j,   nx_loc);
    us = UIDX(i,   j-1, nx_loc);
    un = UIDX(i,   j+1, nx_loc);

    vc = VIDX(i,   j,   nx_loc);
    ve = VIDX(i+1, j,   nx_loc);
    vs = VIDX(i,   j-1, nx_loc);
    vn = VIDX(i,   j+1, nx_loc);

    fdata[uc] = ccu * ydata[uc]
      + cxu * (Wdata[NSPECIES * j] + ydata[ue])
      + cyu * (ydata[us] + ydata[un]);

    fdata[vc] = ccv * ydata[vc]
      + cxv * (Wdata[NSPECIES * j + 1] + ydata[ve])
      + cyv * (ydata[vs] + ydata[vn]);
  }

  // North-West corner
  j = ny_loc - 1;

  uc = UIDX(i,   j,   nx_loc);
  ue = UIDX(i+1, j,   nx_loc);
  us = UIDX(i,   j-1, nx_loc);

  vc = VIDX(i,   j,   nx_loc);
  ve = VIDX(i+1, j,   nx_loc);
  vs = VIDX(i,   j-1, nx_loc);

  fdata[uc] = ccu * ydata[uc]
    + cxu * (Wdata[NSPECIES * j] + ydata[ue])
    + cyu * (ydata[us] + Ndata[NSPECIES * i]);

  fdata[vc] = ccv * ydata[vc]
    + cxv * (Wdata[NSPECIES * j + 1] + ydata[ve])
    + cyv * (ydata[vs] + Ndata[NSPECIES * i + 1]);

  // -----------------------------------------------------
  // East face (updates south-east and north-east corners)
  // -----------------------------------------------------
  i = nx_loc - 1;

  // South-East corner
  j = 0;

  uc = UIDX(i,   j,   nx_loc);
  uw = UIDX(i-1, j,   nx_loc);
  un = UIDX(i,   j+1, nx_loc);

  vc = VIDX(i,   j,   nx_loc);
  vw = VIDX(i-1, j,   nx_loc);
  vn = VIDX(i,   j+1, nx_loc);

  fdata[uc] = ccu * ydata[uc]
    + cxu * (ydata[uw] + Edata[NSPECIES * j])
    + cyu * (Sdata[NSPECIES * i] + ydata[un]);

  fdata[vc] = ccv * ydata[vc]
    + cxv * (ydata[vw] + Edata[NSPECIES * j + 1])
    + cyv * (Sdata[NSPECIES * i + 1] + ydata[vn]);

  // East face interior
  for (j = 1; j < ny_loc - 1; j++)
  {
    uc = UIDX(i,   j,   nx_loc);
    uw = UIDX(i-1, j,   nx_loc);
    us = UIDX(i,   j-1, nx_loc);
    un = UIDX(i,   j+1, nx_loc);

    vc = VIDX(i,   j,   nx_loc);
    vw = VIDX(i-1, j,   nx_loc);
    vs = VIDX(i,   j-1, nx_loc);
    vn = VIDX(i,   j+1, nx_loc);

    fdata[uc] = ccu * ydata[uc]
      + cxu * (ydata[uw] + Edata[NSPECIES * j])
      + cyu * (ydata[us] + ydata[un]);

    fdata[vc] = ccv * ydata[vc]
      + cxv * (ydata[vw] + Edata[NSPECIES * j + 1])
      + cyv * (ydata[vs] + ydata[vn]);
  }

  // North-East corner
  j = ny_loc - 1;

  uc = UIDX(i,   j,   nx_loc);
  uw = UIDX(i-1, j,   nx_loc);
  us = UIDX(i,   j-1, nx_loc);

  vc = VIDX(i,   j,   nx_loc);
  vw = VIDX(i-1, j,   nx_loc);
  vs = VIDX(i,   j-1, nx_loc);

  fdata[uc] = ccu * ydata[uc]
    + cxu * (ydata[uw] + Edata[NSPECIES * j])
    + cyu * (ydata[us] + Ndata[NSPECIES * i]);

  fdata[vc] = ccv * ydata[vc]
    + cxv * (ydata[vw] + Edata[NSPECIES * j + 1])
    + cyv * (ydata[vs] + Ndata[NSPECIES * i + 1]);

  // -----------------------------
  // South face (excludes corners)
  // -----------------------------
  j = 0;

  for (i = 1; i < nx_loc - 1; i++)
  {
    uc = UIDX(i,   j,   nx_loc);
    uw = UIDX(i-1, j,   nx_loc);
    ue = UIDX(i+1, j,   nx_loc);
    un = UIDX(i,   j+1, nx_loc);

    vc = VIDX(i,   j,   nx_loc);
    vw = VIDX(i-1, j,   nx_loc);
    ve = VIDX(i+1, j,   nx_loc);
    vn = VIDX(i,   j+1, nx_loc);

    fdata[uc] = ccu * ydata[uc]
      + cxu * (ydata[uw] + ydata[ue])
      + cyu * (Sdata[NSPECIES * i] + ydata[un]);

    fdata[vc] = ccv * ydata[vc]
      + cxv * (ydata[vw] + ydata[ve])
      + cyv * (Sdata[NSPECIES * i + 1] + ydata[vn]);
  }

  // -----------------------------
  // North face (excludes corners)
  // -----------------------------
  j = ny_loc - 1;

  for (i = 1; i < nx_loc - 1; i++)
  {
    uc = UIDX(i,   j,   nx_loc);
    uw = UIDX(i-1, j,   nx_loc);
    ue = UIDX(i+1, j,   nx_loc);
    us = UIDX(i,   j-1, nx_loc);

    vc = VIDX(i,   j,   nx_loc);
    vw = VIDX(i-1, j,   nx_loc);
    ve = VIDX(i+1, j,   nx_loc);
    vs = VIDX(i,   j-1, nx_loc);

    fdata[uc] = ccu * ydata[uc]
      + cxu * (ydata[uw] + ydata[ue])
      + cyu * (ydata[us] + Ndata[NSPECIES * i]);

    fdata[vc] = ccv * ydata[vc]
      + cxv * (ydata[vw] + ydata[ve])
      + cyv * (ydata[vs] + Ndata[NSPECIES * i + 1]);
  }

  // Stop timer
  udata->rhsD.stop();

  // Return success
  return 0;
}


// Routine to compute the ODE reaction RHS function
static int reaction(realtype t, N_Vector y, N_Vector f, void *user_data)
{
  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Start timer
  udata->rhsR.start();

  // Access data arrays
  realtype *ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) ydata, "N_VGetArrayPointer", 0)) return -1;

  realtype *fdata = N_VGetArrayPointer(f);
  if (check_flag((void *) fdata, "N_VGetArrayPointer", 0)) return -1;

  // Shortcuts to local number of nodes
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Compute reaction term on the subdomain
  realtype u, v;

  for (sunindextype j = 0; j < ny_loc; j++)
  {
    for (sunindextype i = 0; i < nx_loc; i++)
    {
      u = ydata[UIDX(i,j,nx_loc)];
      v = ydata[VIDX(i,j,nx_loc)];

      fdata[UIDX(i,j,nx_loc)] = udata->A + u * u * v - (udata->B + 1) * u;
      fdata[VIDX(i,j,nx_loc)] = udata->B * u - u * u * v;
    }
  }

  // Apply inner forcing for MRI + CVODE
  if (udata->integrator > 1)
  {
    int flag;

    if (udata->integrator == 2)
    {
      // With a global inner stepper the RHS vector f and the forcing vectors
      // from the outer integrator are both MPIPlusX vectors as such we can use
      // a utility function to add the forcing to the RHS vector
      MRIStepInnerStepper_AddForcing(udata->stepper, t, f);
    }
    else if (udata->integrator == 3)
    {
      int      nforcing;
      realtype tshift, tscale;
      N_Vector *forcing;

      // With a local inner stepper the RHS vector f is a serial vector and the
      // forcing vectors from the outer integrator are MPIPlusX vectors as such
      // we need to extract the local serial vectors and apply the forcing
      flag = MRIStepInnerStepper_GetForcingData(udata->stepper,
                                                &tshift, &tscale,
                                                &forcing, &nforcing);
      if (flag != 0) return flag;

      N_Vector forcing_loc;
      realtype tau  = (t - tshift) / tscale;
      realtype taui = ONE;

      for (int i = 0; i < nforcing; i++)
      {
        forcing_loc = N_VGetLocalVector_MPIPlusX(forcing[i]);
        N_VLinearSum(ONE, f, taui, forcing_loc, f);
        taui *= tau;
      }
    }
    else
    {
      cerr << "ERROR: Invalid MRIStep + CVODE option" << endl;
      return -1;
    }
  }

  // Stop timer
  udata->rhsR.stop();

  // Return success
  return 0;
}


// Preconditioner solve routine for Pz = r
static int PSolve(realtype t, N_Vector u, N_Vector f, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data)
{
  // Access user_data structure
  UserData *udata = (UserData *) user_data;

  // Start timer
  udata->psolve.start();

  // Constants for computing diffusion
  realtype cxu = udata->Dux / (udata->dx * udata->dx);
  realtype cyu = udata->Duy / (udata->dy * udata->dy);
  realtype ccu = -TWO * (cxu + cyu);

  realtype cxv = udata->Dvx / (udata->dx * udata->dx);
  realtype cyv = udata->Dvy / (udata->dy * udata->dy);
  realtype ccv = -TWO * (cxv + cyv);

  // Access data arrays
  realtype *rdata = N_VGetArrayPointer(r);
  if (check_flag((void *) rdata, "N_VGetArrayPointer", 0)) return -1;

  realtype *zdata = N_VGetArrayPointer(z);
  if (check_flag((void *) zdata, "N_VGetArrayPointer", 0)) return -1;

  // Shortcuts to local number of nodes
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Set all entries of diag to the inverse diagonal values
  realtype du = ONE / (ONE - gamma * ccu);
  realtype dv = ONE / (ONE - gamma * ccv);

  for (sunindextype j = 0; j < ny_loc; j++)
  {
    for (sunindextype i = 0; i < nx_loc; i++)
    {
      zdata[UIDX(i, j, nx_loc)] = du * rdata[UIDX(i, j, nx_loc)];
      zdata[VIDX(i, j, nx_loc)] = dv * rdata[VIDX(i, j, nx_loc)];
    }
  }

  // Stop timer
  udata->psolve.stop();

  // Return success
  return 0;
}


// -----------------------------------------------------------------------------
// RHS helper functions
// -----------------------------------------------------------------------------


// Open exchange receives
static int StartExchange(N_Vector y, UserData *udata)
{
  int flag;

  // Start timer
  udata->exchange.start();

  // East face (from neighbor's West face)
  flag = MPI_Irecv(udata->Erecv, udata->ybufcount, MPI_SUNREALTYPE,
                   udata->ipE, 0, udata->comm, &(udata->reqRE));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Irecv = " << flag << endl;
    return -1;
  }

  // West face (from neighbor's East face)
  flag = MPI_Irecv(udata->Wrecv, udata->ybufcount, MPI_SUNREALTYPE,
                   udata->ipW, 1, udata->comm, &(udata->reqRW));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Irecv = " << flag << endl;
    return -1;
  }

  // North face (from neighbor's South face)
  flag = MPI_Irecv(udata->Nrecv, udata->xbufcount, MPI_SUNREALTYPE,
                   udata->ipN, 2, udata->comm, &(udata->reqRN));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Irecv = " << flag << endl;
    return -1;
  }

  // South face (from neighbor's North face)
  flag = MPI_Irecv(udata->Srecv, udata->xbufcount, MPI_SUNREALTYPE,
                   udata->ipS, 3, udata->comm, &(udata->reqRS));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Irecv = " << flag << endl;
    return -1;
  }

  // Shortcuts to local number of nodes
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Access data array
  realtype *ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) ydata, "N_VGetArrayPointer", 0)) return -1;

  // Send West face data to neighbor's East face
  for (sunindextype i = 0; i < ny_loc; i++)
  {
    udata->Wsend[NSPECIES * i]     = ydata[UIDX(0,i,nx_loc)];
    udata->Wsend[NSPECIES * i + 1] = ydata[VIDX(0,i,nx_loc)];
  }
  flag = MPI_Isend(udata->Wsend, udata->ybufcount, MPI_SUNREALTYPE,
                   udata->ipW, 0, udata->comm, &(udata->reqSW));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Isend = " << flag << endl;
    return -1;
  }

  // Send East face data to neighbor's West face
  for (sunindextype i = 0; i < ny_loc; i++)
  {
    udata->Esend[NSPECIES * i]     = ydata[UIDX(nx_loc-1,i,nx_loc)];
    udata->Esend[NSPECIES * i + 1] = ydata[VIDX(nx_loc-1,i,nx_loc)];
  }
  flag = MPI_Isend(udata->Esend, udata->ybufcount, MPI_SUNREALTYPE,
                   udata->ipE, 1, udata->comm, &(udata->reqSE));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Isend = " << flag << endl;
    return -1;
  }

  // Send South face data to neighbor's North face
  for (sunindextype i = 0; i < nx_loc; i++)
  {
    udata->Ssend[NSPECIES * i]     = ydata[UIDX(i,0,nx_loc)];
    udata->Ssend[NSPECIES * i + 1] = ydata[VIDX(i,0,nx_loc)];
  }
  flag = MPI_Isend(udata->Ssend, udata->xbufcount, MPI_SUNREALTYPE,
                   udata->ipS, 2, udata->comm, &(udata->reqSS));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Isend = " << flag << endl;
    return -1;
  }

  // Send North face data to neighbor's South face
  for (sunindextype i = 0; i < nx_loc; i++)
  {
    udata->Nsend[NSPECIES * i]     = ydata[UIDX(i,ny_loc-1,nx_loc)];
    udata->Nsend[NSPECIES * i + 1] = ydata[VIDX(i,ny_loc-1,nx_loc)];
  }
  flag = MPI_Isend(udata->Nsend, udata->xbufcount, MPI_SUNREALTYPE,
                   udata->ipN, 3, udata->comm, &(udata->reqSN));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Isend = " << flag << endl;
    return -1;
  }

  // Stop timer
  udata->exchange.stop();

  // Return success
  return 0;
}


// Wait for exchange data
static int EndExchange(UserData *udata)
{
  // Local variables
  int flag;
  MPI_Status stat;

  // Start timer
  udata->exchange.start();

  // Wait for messages to finish
  flag = MPI_Wait(&(udata->reqRW), &stat);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Wait = " << flag << endl;
    return -1;
  }
  flag = MPI_Wait(&(udata->reqSW), &stat);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Wait = " << flag << endl;
    return -1;
  }

  flag = MPI_Wait(&(udata->reqRE), &stat);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Wait = " << flag << endl;
    return -1;
  }
  flag = MPI_Wait(&(udata->reqSE), &stat);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Wait = " << flag << endl;
    return -1;
  }

  flag = MPI_Wait(&(udata->reqRS), &stat);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Wait = " << flag << endl;
    return -1;
  }
  flag = MPI_Wait(&(udata->reqSS), &stat);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Wait = " << flag << endl;
    return -1;
  }

  flag = MPI_Wait(&(udata->reqRN), &stat);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Wait = " << flag << endl;
    return -1;
  }
  flag = MPI_Wait(&(udata->reqSN), &stat);
  if (flag != MPI_SUCCESS)
  {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
  }

  // Stop timer
  udata->exchange.stop();

  // Return success
  return 0;
}


// Send exchange data
static int ExchangeBC(N_Vector y, UserData *udata)
{
  int flag;
  MPI_Status stat;

  // Shortcuts to local number of nodes
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Post East face exchange receives
  if (udata->ie == udata->nx - 1)
  {
    flag = MPI_Irecv(udata->Erecv, udata->ybufcount, MPI_SUNREALTYPE,
                     udata->ipE, MPI_ANY_TAG, udata->comm, &(udata->reqRE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  // Post North face exchange receives
  if (udata->je == udata->ny - 1)
  {
    flag = MPI_Irecv(udata->Nrecv, udata->xbufcount, MPI_SUNREALTYPE,
                     udata->ipN, MPI_ANY_TAG, udata->comm, &(udata->reqRN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  // Post North-East corner exchange receive
  if (udata->ie == udata->nx - 1 && udata->je == udata->ny - 1)
  {
    flag = MPI_Irecv(udata->NErecv, NSPECIES, MPI_SUNREALTYPE,
                     udata->ipNE, MPI_ANY_TAG, udata->comm, &(udata->reqRC));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  realtype *ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) ydata, "N_VGetArrayPointer", 0)) return -1;

  // Send West face data
  if (udata->is == 0)
  {
    for (sunindextype i = 0; i < ny_loc; i++)
    {
      udata->Wsend[NSPECIES * i]     = ydata[UIDX(0,i,nx_loc)];
      udata->Wsend[NSPECIES * i + 1] = ydata[VIDX(0,i,nx_loc)];
    }
    flag = MPI_Isend(udata->Wsend, udata->ybufcount, MPI_SUNREALTYPE,
                     udata->ipW, 0, udata->comm, &(udata->reqSW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  // Send South face data
  if (udata->js == 0)
  {
    for (sunindextype i = 0; i < nx_loc; i++)
    {
      udata->Ssend[NSPECIES * i]     = ydata[UIDX(i,0,nx_loc)];
      udata->Ssend[NSPECIES * i + 1] = ydata[VIDX(i,0,nx_loc)];
    }
    flag = MPI_Isend(udata->Ssend, udata->xbufcount, MPI_SUNREALTYPE,
                     udata->ipS, 2, udata->comm, &(udata->reqSS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  // Send South-West corner data
  if (udata->is == 0 && udata->js == 0)
  {
    udata->SWsend[0] = ydata[UIDX(0,0,nx_loc)];
    udata->SWsend[1] = ydata[VIDX(0,0,nx_loc)];
    flag = MPI_Isend(udata->SWsend, NSPECIES, MPI_SUNREALTYPE,
                     udata->ipSW, 2, udata->comm, &(udata->reqSC));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  // Wait for messages to finish
  if (udata->ie == udata->nx - 1)
  {
    flag = MPI_Wait(&(udata->reqRE), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->je == udata->ny - 1)
  {
    flag = MPI_Wait(&(udata->reqRN), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->ie == udata->nx - 1 && udata->je == udata->ny - 1)
  {
    flag = MPI_Wait(&(udata->reqRC), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->is == 0)
  {
    flag = MPI_Wait(&(udata->reqSW), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->js == 0)
  {
    flag = MPI_Wait(&(udata->reqSS), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->is == 0 && udata->js == 0)
  {
    flag = MPI_Wait(&(udata->reqSC), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  return(0);
}


// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------


// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{
  // Free exchange buffers
  if (udata->Wrecv  != NULL) delete[] udata->Wrecv;
  if (udata->Wsend  != NULL) delete[] udata->Wsend;
  if (udata->Erecv  != NULL) delete[] udata->Erecv;
  if (udata->Esend  != NULL) delete[] udata->Esend;
  if (udata->Srecv  != NULL) delete[] udata->Srecv;
  if (udata->Ssend  != NULL) delete[] udata->Ssend;
  if (udata->Nrecv  != NULL) delete[] udata->Nrecv;
  if (udata->Nsend  != NULL) delete[] udata->Nsend;
  if (udata->NErecv != NULL) delete[] udata->NErecv;
  if (udata->SWsend != NULL) delete[] udata->SWsend;

  // Free preconditioner data
  if (udata->diag)
  {
    N_VDestroy(udata->diag);
    udata->diag = NULL;
  }

  // Free MPI Cartesian communicator
  if (udata->comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&(udata->comm));
    udata->comm = MPI_COMM_NULL;
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
    // MPI processes
    else if (arg == "--np")
    {
      udata->npx = stoi((*argv)[arg_idx++]);
      udata->npy = stoi((*argv)[arg_idx++]);
    }
    // Domain bounds
    else if (arg == "--domain")
    {
      udata->xl = stoi((*argv)[arg_idx++]);
      udata->xu = stoi((*argv)[arg_idx++]);
      udata->yl = stoi((*argv)[arg_idx++]);
      udata->yu = stoi((*argv)[arg_idx++]);
    }
    // Diffusion parameters
    else if (arg == "--D")
    {
      udata->Dux = stod((*argv)[arg_idx++]);
      udata->Duy = stod((*argv)[arg_idx++]);
      udata->Dvx = stod((*argv)[arg_idx++]);
      udata->Dvy = stod((*argv)[arg_idx++]);
    }
    // Reaction parameters
    else if (arg == "--A")
    {
      udata->A = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--B")
    {
      udata->B = stod((*argv)[arg_idx++]);
    }
    // Temporal domain settings
    else if (arg == "--tf")
    {
      udata->tf = stod((*argv)[arg_idx++]);
    }
    // Integrator options
    else if (arg == "--imex")
    {
      udata->integrator = 0;
    }
    else if (arg == "--mri-arkstep")
    {
      udata->integrator = 1;
    }
    else if (arg == "--mri-cvode-global")
    {
      udata->integrator = 2;
    }
    else if (arg == "--mri-cvode-local")
    {
      udata->integrator = 3;
    }
    // IMEX integrator settings
    else if (arg == "--rtol_imex")
    {
      udata->rtol_imex = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--atol_imex")
    {
      udata->atol_imex = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--h_imex")
    {
      udata->h_imex = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--order_imex")
    {
      udata->order_imex = stoi((*argv)[arg_idx++]);
    }
    // MRI integrator settings
    else if (arg == "--rtol_slow")
    {
      udata->rtol_fast = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--atol_slow")
    {
      udata->atol_fast = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--rtol_fast")
    {
      udata->rtol_fast = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--atol_fast")
    {
      udata->atol_fast = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--h_slow")
    {
      udata->h_slow = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--h_fast")
    {
      udata->h_fast = stod((*argv)[arg_idx++]);
    }
    // Shared IMEX and MRI settings
    else if (arg == "--controller")
    {
      udata->controller = stoi((*argv)[arg_idx++]);
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
    // Debugging
    else if (arg == "--onestep")
    {
      udata->onestep = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--no_diffusion")
    {
      udata->diffusion = false;
    }
    else if (arg == "--no_reaction")
    {
      udata->reaction = false;
    }
    // Help
    else if (arg == "--help")
    {
      if (udata->outproc) InputHelp();
      return -1;
    }
    // Unknown input
    else
    {
      if (udata->outproc)
      {
        cerr << "ERROR: Invalid input " << arg << endl;
        InputHelp();
      }
      return -1;
    }
  }

  // Recompute total number of equations
  udata->neq = NSPECIES * udata->nx * udata->ny;

  // Recompute x and y mesh spacing with periodic boundary conditions
  udata->dx = (udata->xu - udata->xl) / udata->nx;
  udata->dy = (udata->yu - udata->yl) / udata->ny;

  // Compute slow step size based on CFL if not set by input
  if (udata->h_slow < ZERO)
  {
    realtype cfl_u = RCONST(0.5) / ((udata->Dux / (udata->dx * udata->dx)) +
                                    (udata->Duy / (udata->dy * udata->dy)));
    realtype cfl_v = RCONST(0.5) / ((udata->Dvx / (udata->dx * udata->dx)) +
                                    (udata->Dvy / (udata->dy * udata->dy)));
    udata->h_slow = RCONST(5.0) * min(cfl_u, cfl_v);
  }

  // Return success
  return 0;
}


// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------


// Compute the initial condition
static int SetIC(N_Vector u, UserData *udata)
{
  realtype x, y, a, b;

  // Shortcuts to local number of nodes
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Gaussian random number generator
  default_random_engine generator;
  normal_distribution<double> dist(RCONST(0.0), RCONST(0.001));

  realtype *data = N_VGetArrayPointer(u);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return -1;

  for (sunindextype j = 0; j < ny_loc; j++)
  {
    for (sunindextype i = 0; i < nx_loc; i++)
    {
      x = udata->xl + (udata->is + i) * udata->dx;
      y = udata->yl + (udata->js + j) * udata->dy;

      a = TWO * PI * (x - udata->xl) / (udata->xu - udata->xl);
      b = TWO * PI * (y - udata->yl) / (udata->yu - udata->yl);

      data[UIDX(i,j,nx_loc)] = udata->A + RCONST(0.5) * sin(a) * sin(b);
      data[VIDX(i,j,nx_loc)] = udata->B / udata->A;
    }
  }

  return 0;
}


// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --mesh <nx> <ny>             : number of mesh points"       << endl;
  cout << "  --np <npx> <npy>             : number of MPI processes"     << endl;
  cout << "  --domain <xl> <xu> <yl> <yu> : domain boundaries"           << endl;
  cout << "  --D <Dux> <Duy> <Dvx> <Dvy>  : diffusion coefficients"      << endl;
  cout << "  --A <A>                      : species A concentration"     << endl;
  cout << "  --B <A>                      : species B concentration"     << endl;
  cout << "  --tf <time>                  : final time"                  << endl;
  cout << "  --imex                       : use an IMEX method"          << endl;
  cout << "  --mri-arkstep                : use MRI with ARKStep"        << endl;
  cout << "  --mri-cvode-global           : use MRI with CVODE global stepper"     << endl;
  cout << "  --mri-cvode-local            : use MRI with CVODE task-local stepper" << endl;
  cout << "  --rtol_imex <rtol>           : IMEX relative tolerance"     << endl;
  cout << "  --atol_imex <atol>           : IMEX absoltue tolerance"     << endl;
  cout << "  --h_imex <h>                 : IMEX fixed step size"        << endl;
  cout << "  --order_imex <ord>           : IMEX method order"           << endl;
  cout << "  --rtol_slow <rtol>           : MRI slow relative tolerance" << endl;
  cout << "  --atol_slow <atol>           : MRI slow absoltue tolerance" << endl;
  cout << "  --h_slow <h>                 : MRI slow step size"          << endl;
  cout << "  --rtol_fast <rtol>           : MRI fast relative tolerance" << endl;
  cout << "  --atol_fast <atol>           : MRI fast absoltue tolerance" << endl;
  cout << "  --h_fast <h>                 : MRI fast step size"          << endl;
  cout << "  --controller <ctr>           : time step adaptivity"        << endl;
  cout << "  --nonlinear                  : nonlinearly implicit"        << endl;
  cout << "  --diagnostics                : output diagnostics"          << endl;
  cout << "  --gmres                      : use GMRES linear solver"     << endl;
  cout << "  --lsinfo                     : output residual history"     << endl;
  cout << "  --liniters <iters>           : max number of iterations"    << endl;
  cout << "  --epslin <factor>            : linear tolerance factor"     << endl;
  cout << "  --noprec                     : disable preconditioner"      << endl;
  cout << "  --msbp <steps>               : prec setup frequency"        << endl;
  cout << "  --output <level>             : output level"                << endl;
  cout << "  --nout <nout>                : number of outputs"           << endl;
  cout << "  --maxsteps <steps>           : max steps between outputs"   << endl;
  cout << "  --timing                     : print timing data"           << endl;
  cout << "  --onestep <steps>            : fixed number of steps"       << endl;
  cout << "  --nodiffusion                : no diffusion (IMEX only)"    << endl;
  cout << "  --noreaction                 : no reactions (IMEX only)"    << endl;
  cout << "  --help                       : print options and exit"      << endl;
}


// Print user data
static int PrintUserData(UserData *udata)
{
  cout << endl;
  cout << "2D Heat PDE test problem:"                << endl;
  cout << " --------------------------------- "      << endl;
  cout << "  nprocs         = " << udata->nprocs     << endl;
  cout << "  npx            = " << udata->npx        << endl;
  cout << "  npy            = " << udata->npy        << endl;
  cout << " --------------------------------- "      << endl;
  cout << "  Dux            = " << udata->Dux        << endl;
  cout << "  Duy            = " << udata->Duy        << endl;
  cout << "  Dvx            = " << udata->Dvx        << endl;
  cout << "  Dvy            = " << udata->Dvy        << endl;
  cout << "  A              = " << udata->A          << endl;
  cout << "  B              = " << udata->B          << endl;
  cout << " --------------------------------- "      << endl;
  cout << "  tf             = " << udata->tf         << endl;
  cout << "  xl             = " << udata->xl         << endl;
  cout << "  xu             = " << udata->xu         << endl;
  cout << "  yl             = " << udata->yl         << endl;
  cout << "  yu             = " << udata->yu         << endl;
  cout << " --------------------------------- "      << endl;
  cout << "  nx             = " << udata->nx         << endl;
  cout << "  ny             = " << udata->ny         << endl;
  cout << "  dx             = " << udata->dx         << endl;
  cout << "  dy             = " << udata->dy         << endl;
  cout << "  nxl (proc 0)   = " << udata->nx_loc     << endl;
  cout << "  nyl (proc 0)   = " << udata->ny_loc     << endl;
  cout << "  is  (proc 0)   = " << udata->is         << endl;
  cout << "  ie  (proc 0)   = " << udata->ie         << endl;
  cout << "  je  (proc 0)   = " << udata->js         << endl;
  cout << "  je  (proc 0)   = " << udata->je         << endl;
  cout << " --------------------------------- "      << endl;
  if (udata->integrator)
  {
    cout << "  rtol_slow      = " << udata->rtol_slow  << endl;
    cout << "  atol_slow      = " << udata->atol_slow  << endl;
    cout << "  rtol_fast      = " << udata->rtol_fast  << endl;
    cout << "  atol_fast      = " << udata->atol_fast  << endl;
    cout << "  order_fast     = " << udata->order_fast << endl;
    cout << "  fixed h slow   = " << udata->h_slow     << endl;
    cout << "  fixed h fast   = " << udata->h_fast     << endl;
  }
  else
  {
    cout << "  rtol           = " << udata->rtol_imex  << endl;
    cout << "  atol           = " << udata->atol_imex  << endl;
    cout << "  order          = " << udata->order_imex << endl;
    cout << "  fixed h        = " << udata->h_imex     << endl;
  }
  cout << "  controller     = " << udata->controller << endl;
  cout << "  linear         = " << udata->linear     << endl;
  cout << " --------------------------------- "      << endl;
  if (udata->pcg)
  {
    cout << "  linear solver  = PCG" << endl;
  }
  else
  {
    cout << "  linear solver  = GMRES" << endl;
  }
  cout << "  lin iters      = " << udata->liniters << endl;
  cout << "  eps lin        = " << udata->epslin   << endl;
  cout << "  prec           = " << udata->prec     << endl;
  cout << "  msbp           = " << udata->msbp     << endl;
  cout << " --------------------------------- "    << endl;
  cout << "  output         = " << udata->output   << endl;
  cout << " --------------------------------- "    << endl;
  cout << endl;

  return 0;
}


// Initialize output
static int OpenOutput(UserData *udata)
{
  // Header for status output
  if (udata->output > 0 && udata->outproc)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<realtype>::digits10);
    cout << "          t           ";
    cout << "          ||u||_rms      " << endl;
    cout << " ---------------------";
    cout << "-------------------------" << endl;
  }

  // Open output stream and output problem information
  if (udata->output == 2)
  {
    // Open output stream
    stringstream fname;
    fname << "diffusion_reaction." << setfill('0') << setw(5)
          << udata->myid << ".out";
    udata->uout.open(fname.str());

    udata->uout << scientific;
    udata->uout << setprecision(numeric_limits<realtype>::digits10);

    // Add 1 to the total number of nodes in the x and y directions and to the
    // end indices in the x and y direction at the North and East boundary to
    // account for additional output at the periodic boundary
    udata->uout <<  "# title Diffusion-Reaction (Brusselator)" << endl;
    udata->uout <<  "# nprocs "  << udata->nprocs << endl;
    udata->uout <<  "# npx "     << udata->npx << endl;
    udata->uout <<  "# npy "     << udata->npy << endl;
    udata->uout <<  "# nvar 2"   << endl;
    udata->uout <<  "# vars u v" << endl;
    udata->uout <<  "# nt "      << udata->nout + 1 << endl;
    udata->uout <<  "# nx "      << udata->nx + 1 << endl;
    udata->uout <<  "# xl "      << udata->xl << endl;
    udata->uout <<  "# xu "      << udata->xu << endl;
    udata->uout <<  "# is "      << udata->is << endl;
    if (udata->ie == udata->nx - 1)
      udata->uout <<  "# ie " << udata->ie + 1 << endl;
    else
      udata->uout <<  "# ie " << udata->ie << endl;
    udata->uout <<  "# ny "      << udata->ny + 1 << endl;
    udata->uout <<  "# yl "      << udata->yl << endl;
    udata->uout <<  "# yu "      << udata->yu << endl;
    udata->uout <<  "# js " << udata->js << endl;
    if (udata->je == udata->ny - 1)
      udata->uout <<  "# je " << udata->je + 1 << endl;
    else
      udata->uout <<  "# je " << udata->je << endl;
  }

  return 0;
}


// Write output
static int WriteOutput(realtype t, N_Vector y, UserData *udata)
{
  int flag;

  if (udata->output > 0)
  {
    // Compute rms norm of the state
    realtype urms = sqrt(N_VDotProd(y, y) / udata->nx / udata->ny);

    // Output current status
    if (udata->outproc) cout << setw(22) << t << setw(25) << urms << endl;

    // Write solution to disk
    if (udata->output == 2)
    {
      // Shortcuts to local number of nodes
      sunindextype ny_loc = udata->ny_loc;
      sunindextype nx_loc = udata->nx_loc;

      flag = ExchangeBC(y, udata);
      if (check_flag(&flag, "ExchangeBC", 1)) return -1;

      realtype *ydata = N_VGetArrayPointer(y);
      if (check_flag((void *) ydata, "N_VGetArrayPointer", 0)) return -1;

      udata->uout << t;
      for (sunindextype j = 0; j < ny_loc; j++)
      {
        for (sunindextype i = 0; i < nx_loc; i++)
        {
          udata->uout << setw(WIDTH) << ydata[UIDX(i,j,nx_loc)];
          udata->uout << setw(WIDTH) << ydata[VIDX(i,j,nx_loc)];
        }
        // East boundary (same as West face)
        if (udata->ie == udata->nx - 1)
        {
          udata->uout << setw(WIDTH) << udata->Erecv[NSPECIES * j];
          udata->uout << setw(WIDTH) << udata->Erecv[NSPECIES * j + 1];
        }
      }
      // North boundary (same as South face)
      if (udata->je == udata->ny - 1)
      {
        for (sunindextype i = 0; i < udata->nx_loc; i++)
        {
          udata->uout << setw(WIDTH) << udata->Nrecv[NSPECIES * i];
          udata->uout << setw(WIDTH) << udata->Nrecv[NSPECIES * i + 1];
        }
        // North-East corner (same as South-West corner)
        if (udata->ie == udata->nx - 1)
        {
          udata->uout << setw(WIDTH) << udata->NErecv[0];
          udata->uout << setw(WIDTH) << udata->NErecv[1];
        }
      }
      udata->uout << endl;
    }
  }

  return 0;
}


// Finalize output
static int CloseOutput(UserData *udata)
{
  // Footer for status output
  if (udata->outproc && (udata->output > 0))
  {
    cout << " ---------------------";
    cout << "-------------------------" << endl;
    cout << endl;
  }

  if (udata->output == 2)
  {
    // Close output streams
    udata->uout.close();
  }

  return 0;
}

// Print integrator statistics
static int OutputStatsIMEX(void *arkode_mem, UserData* udata)
{
  int flag;

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

  if (udata->diffusion)
  {
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
  }

  cout << fixed;
  cout << setprecision(6);

  cout << "  Steps            = " << nst     << endl;
  cout << "  Step attempts    = " << nst_a   << endl;
  cout << "  Error test fails = " << netf    << endl;
  if (udata->reaction)
  {
    cout << "  RHS reaction     = " << nfe     << endl;
  }
  if (udata->diffusion)
  {
    cout << "  RHS diffusion    = " << nfi     << endl;
    cout << "  NLS iters        = " << nni     << endl;
    cout << "  NLS fails        = " << ncfn    << endl;
    cout << "  LS iters         = " << nli     << endl;
    cout << "  LS fails         = " << nlcf    << endl;
    cout << "  LS setups        = " << nsetups << endl;
    cout << "  LS RHS evals     = " << nfi_ls  << endl;
    cout << "  Jv products      = " << nJv     << endl;
  }
  cout << endl;

  if (udata->diffusion)
  {
    // Compute average nls iters per step attempt and ls iters per nls iter
    realtype avgnli = (realtype) nni / (realtype) nst_a;
    realtype avgli  = (realtype) nli / (realtype) nni;
    cout << "  Avg NLS iters per step attempt = " << avgnli << endl;
    cout << "  Avg LS iters per NLS iter      = " << avgli  << endl;
    cout << endl;

    // Get preconditioner stats
    if (udata->prec)
    {
      long int npe, nps;
      flag = ARKStepGetNumPrecEvals(arkode_mem, &npe);
      if (check_flag(&flag, "ARKStepGetNumPrecEvals", 1)) return -1;
      flag = ARKStepGetNumPrecSolves(arkode_mem, &nps);
      if (check_flag(&flag, "ARKStepGetNumPrecSolves", 1)) return -1;

      cout << "  Preconditioner setups = " << npe << endl;
      cout << "  Preconditioner solves = " << nps << endl;
      cout << endl;
    }
  }

  return 0;
}


// Print integrator statistics
static int OutputStatsMRI(void *arkode_mem, MRIStepInnerStepper stepper,
                          UserData* udata)
{
  int flag;

  // Get slow integrator and solver stats
  long int nsts, nfse, nfsi, nni, ncfn, nli, nlcf, nsetups, nfi_ls, nJv;
  flag = MRIStepGetNumSteps(arkode_mem, &nsts);
  if (check_flag(&flag, "MRIStepGetNumSteps", 1)) return -1;
  flag = MRIStepGetNumRhsEvals(arkode_mem, &nfse, &nfsi);
  if (check_flag(&flag, "MRIStepGetNumRhsEvals", 1)) return -1;
  flag = MRIStepGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(&flag, "MRIStepGetNumNonlinSolvIters", 1)) return -1;
  flag = MRIStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(&flag, "MRIStepGetNumNonlinSolvConvFails", 1)) return -1;
  flag = MRIStepGetNumLinIters(arkode_mem, &nli);
  if (check_flag(&flag, "MRIStepGetNumLinIters", 1)) return -1;
  flag = MRIStepGetNumLinConvFails(arkode_mem, &nlcf);
  if (check_flag(&flag, "MRIStepGetNumLinConvFails", 1)) return -1;
  flag = MRIStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(&flag, "MRIStepGetNumLinSolvSetups", 1)) return -1;
  flag = MRIStepGetNumLinRhsEvals(arkode_mem, &nfi_ls);
  if (check_flag(&flag, "MRIStepGetNumLinRhsEvals", 1)) return -1;
  flag = MRIStepGetNumJtimesEvals(arkode_mem, &nJv);
  if (check_flag(&flag, "MRIStepGetNumJtimesEvals", 1)) return -1;

  cout << fixed;
  cout << setprecision(6);

  cout << endl << "Slow Integrator:" << endl;

  cout << "  Steps            = " << nsts    << endl;
  cout << "  RHS diffusion    = " << nfsi    << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << "  NLS fails        = " << ncfn    << endl;
  cout << "  LS iters         = " << nli     << endl;
  cout << "  LS fails         = " << nlcf    << endl;
  cout << "  LS setups        = " << nsetups << endl;
  cout << "  LS RHS evals     = " << nfi_ls  << endl;
  cout << "  Jv products      = " << nJv     << endl;
  cout << endl;

  // Compute average nls iters per step and ls iters per nls iter
  realtype avgnli = (realtype) nni / (realtype) nsts;
  realtype avgli  = (realtype) nli / (realtype) nni;
  cout << "  Avg NLS iters per step attempt = " << avgnli << endl;
  cout << "  Avg LS iters per NLS iter      = " << avgli  << endl;
  cout << endl;

  // Get preconditioner stats
  if (udata->prec)
  {
    long int npe, nps;
    flag = MRIStepGetNumPrecEvals(arkode_mem, &npe);
    if (check_flag(&flag, "MRIStepGetNumPrecEvals", 1)) return -1;
    flag = MRIStepGetNumPrecSolves(arkode_mem, &nps);
    if (check_flag(&flag, "MRIStepGetNumPrecSolves", 1)) return -1;

    cout << "  Preconditioner setups = " << npe << endl;
    cout << "  Preconditioner solves = " << nps << endl;
    cout << endl;
  }

  // Get fast integrator stats
  void* inner_arkode_mem;
  MRIStepInnerStepper_GetContent(stepper, &inner_arkode_mem);

  long int nstf, nstf_a, netff, nffe, nffi;

  flag = ARKStepGetNumSteps(inner_arkode_mem, &nstf);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return -1;
  flag = ARKStepGetNumStepAttempts(inner_arkode_mem, &nstf_a);
  if (check_flag(&flag, "ARKStepGetNumStepAttempts", 1)) return -1;
  flag = ARKStepGetNumErrTestFails(inner_arkode_mem, &netff);
  if (check_flag(&flag, "ARKStepGetNumErrTestFails", 1)) return -1;
  flag = ARKStepGetNumRhsEvals(inner_arkode_mem, &nffe, &nffi);
  if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return -1;

  cout << "Fast Integrator:" << endl;
  cout << "  Steps            = " << nstf   << endl;
  cout << "  Step attempts    = " << nstf_a << endl;
  cout << "  Error test fails = " << netff  << endl;
  cout << "  RHS reaction     = " << nffe   << endl;

  return 0;
}


// Print integrator statistics
static int OutputStatsMRICVODE(void *arkode_mem,
                               MRIStepInnerStepper stepper,
                               UserData* udata)
{
  int flag;

  // Get slow integrator and solver stats
  long int nsts, nfse, nfsi, nni, ncfn, nli, nlcf, nsetups, nfi_ls, nJv;
  flag = MRIStepGetNumSteps(arkode_mem, &nsts);
  if (check_flag(&flag, "MRIStepGetNumSteps", 1)) return -1;
  flag = MRIStepGetNumRhsEvals(arkode_mem, &nfse, &nfsi);
  if (check_flag(&flag, "MRIStepGetNumRhsEvals", 1)) return -1;
  flag = MRIStepGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(&flag, "MRIStepGetNumNonlinSolvIters", 1)) return -1;
  flag = MRIStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(&flag, "MRIStepGetNumNonlinSolvConvFails", 1)) return -1;
  flag = MRIStepGetNumLinIters(arkode_mem, &nli);
  if (check_flag(&flag, "MRIStepGetNumLinIters", 1)) return -1;
  flag = MRIStepGetNumLinConvFails(arkode_mem, &nlcf);
  if (check_flag(&flag, "MRIStepGetNumLinConvFails", 1)) return -1;
  flag = MRIStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(&flag, "MRIStepGetNumLinSolvSetups", 1)) return -1;
  flag = MRIStepGetNumLinRhsEvals(arkode_mem, &nfi_ls);
  if (check_flag(&flag, "MRIStepGetNumLinRhsEvals", 1)) return -1;
  flag = MRIStepGetNumJtimesEvals(arkode_mem, &nJv);
  if (check_flag(&flag, "MRIStepGetNumJtimesEvals", 1)) return -1;

  cout << fixed;
  cout << setprecision(6);

  cout << endl << "Slow Integrator:" << endl;

  cout << "  Steps            = " << nsts    << endl;
  cout << "  RHS diffusion    = " << nfsi    << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << "  NLS fails        = " << ncfn    << endl;
  cout << "  LS iters         = " << nli     << endl;
  cout << "  LS fails         = " << nlcf    << endl;
  cout << "  LS setups        = " << nsetups << endl;
  cout << "  LS RHS evals     = " << nfi_ls  << endl;
  cout << "  Jv products      = " << nJv     << endl;
  cout << endl;

  // Compute average nls iters per step and ls iters per nls iter
  realtype avgnli = (realtype) nni / (realtype) nsts;
  realtype avgli  = (realtype) nli / (realtype) nni;
  cout << "  Avg NLS iters per step attempt = " << avgnli << endl;
  cout << "  Avg LS iters per NLS iter      = " << avgli  << endl;
  cout << endl;

  // Get preconditioner stats
  if (udata->prec)
  {
    long int npe, nps;
    flag = MRIStepGetNumPrecEvals(arkode_mem, &npe);
    if (check_flag(&flag, "MRIStepGetNumPrecEvals", 1)) return -1;
    flag = MRIStepGetNumPrecSolves(arkode_mem, &nps);
    if (check_flag(&flag, "MRIStepGetNumPrecSolves", 1)) return -1;

    cout << "  Preconditioner setups = " << npe << endl;
    cout << "  Preconditioner solves = " << nps << endl;
    cout << endl;
  }

  // Get fast integrator stats and solver stats
  void* inner_content;
  MRIStepInnerStepper_GetContent(stepper, &inner_content);
  InnerStepperContent* content = (InnerStepperContent *) inner_content;

  cout << fixed;
  cout << setprecision(6);

  cout << "Fast Integrator:" << endl;
  cout << "  Steps            = " << content->nst  << endl;
  cout << "  Error test fails = " << content->netf << endl;
  cout << "  RHS reaction     = " << content->nfe  << endl;
  cout << "  NLS iters        = " << content->nni  << endl;
  cout << "  NLS fails        = " << content->nncf << endl;
  cout << endl;

  return 0;
}


// Ouput timing stats
static int OutputTiming(UserData *udata)
{
  if (udata->outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double max = 0.0;

  MPI_Reduce(&(udata->evolve), &max, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm);
  if (udata->outproc)
  {
    cout << "  Evolve time   = " << max << " sec" << endl;
  }

  MPI_Reduce(&(udata->rhsD), &max, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm);
  if (udata->outproc)
  {
    cout << "  Diffusion RHS time = " << max << " sec" << endl;
  }

  MPI_Reduce(&(udata->rhsD), &max, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm);
  if (udata->outproc)
  {
    cout << "  Reaction RHS time = " << max << " sec" << endl;
  }

  MPI_Reduce(&(udata->exchange), &max, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm);
  if (udata->outproc)
  {
    cout << "  Exchange time = " << max << " sec" << endl;
    cout << endl;
  }

  if (udata->prec)
  {
    MPI_Reduce(&(udata->psolve), &max, 1, MPI_DOUBLE, MPI_MAX, 0,
               udata->comm);
    if (udata->outproc)
    {
      cout << "  PSolve time   = " << max << " sec" << endl;
      cout << endl;
    }
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
