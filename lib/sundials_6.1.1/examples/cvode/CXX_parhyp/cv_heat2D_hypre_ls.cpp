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
 * paired with the hypre's PCG or GMRES linear solver and PFMG preconditioner.
 * Several command line options are available to change the problem parameters
 * and CVODE settings. Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>

#include "cvode/cvode.h"                     // access to CVODE
#include "nvector/nvector_parallel.h"        // access to the MPI N_Vector
#include "sundials/sundials_linearsolver.h"  // definition SUNLinearSolver
#include "sundials/sundials_matrix.h"        // definition SUNMatrix
#include "HYPRE_struct_ls.h"                 // HYPRE structured grid solver interface
#include "mpi.h"                             // MPI header file


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
  UserData(sundials::Context&);
  ~UserData();

  // SUNDIALS simulation context
  sundials::Context& sunctx;

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

  // Integrator settings
  realtype rtol;        // relative tolerance
  realtype atol;        // absolute tolerance
  int      maxsteps;    // max number of steps between outputs

  // Linear solver and preconditioner settings
  bool     pcg;       // use PCG (true) or GMRES (false)
  bool     prec;      // preconditioner on/off
  int      liniters;  // number of linear iterations
  int      msbp;      // max number of steps between preconditioner setups
  realtype epslin;    // linear solver tolerance factor

  // hypre PFMG settings (hypre defaults)
  HYPRE_Int pfmg_relax;  // type of relaxation:
                         //   0 - Jacobi
                         //   1 - Weighted Jacobi
                         //   2 - symmetric R/B Gauss-Seidel (*)
                         //   3 - nonsymmetric R/B Gauss-Seidel
  HYPRE_Int pfmg_nrelax; // number of pre and post relaxation sweeps (2)

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
  double matfilltime;
  double setuptime;
  double solvetime;
  double exchangetime;
};

// -----------------------------------------------------------------------------
// Custom hypre 5-point structured grid matrix definition
// -----------------------------------------------------------------------------

struct Hypre5ptMatrixContent
{
  // hypre objects
  HYPRE_StructGrid    grid;
  HYPRE_StructStencil stencil;
  HYPRE_StructMatrix  matrix;

  // hypre grid extents
  HYPRE_Int ilower[2];
  HYPRE_Int iupper[2];

  // hypre workspace
  HYPRE_Int   nwork;
  HYPRE_Real *work;

  // User data
  UserData *udata;
};

// Accessor macros
#define H5PM_CONTENT(A)  ( (Hypre5ptMatrixContent*)(A->content) )
#define H5PM_ILOWER(A)   ( H5PM_CONTENT(A)->ilower )
#define H5PM_IUPPER(A)   ( H5PM_CONTENT(A)->iupper )
#define H5PM_GRID(A)     ( H5PM_CONTENT(A)->grid )
#define H5PM_STENCIL(A)  ( H5PM_CONTENT(A)->stencil )
#define H5PM_MATRIX(A)   ( H5PM_CONTENT(A)->matrix )
#define H5PM_WORK(A)     ( H5PM_CONTENT(A)->work )
#define H5PM_NWORK(A)    ( H5PM_CONTENT(A)->nwork )
#define H5PM_UDATA(A)    ( H5PM_CONTENT(A)->udata )

// Matrix function prototypes
SUNMatrix Hypre5ptMatrix(UserData *udata);
SUNMatrix_ID Hypre5ptMatrix_GetID(SUNMatrix A);
SUNMatrix Hypre5ptMatrix_Clone(SUNMatrix A);
void Hypre5ptMatrix_Destroy(SUNMatrix A);
int Hypre5ptMatrix_Copy(SUNMatrix A, SUNMatrix B);
int Hypre5ptMatrix_ScaleAddI(realtype c, SUNMatrix A);

// -----------------------------------------------------------------------------
// Custom hypre linear solver definition
// -----------------------------------------------------------------------------

struct HypreLSContent
{
  // hypre objects
  HYPRE_StructVector bvec;
  HYPRE_StructVector xvec;
  HYPRE_StructSolver precond;
  HYPRE_StructSolver solver;

  // Counters
  HYPRE_Int iters;

  // PCG (true) or GMRES (false)
  bool pcg;
};

// Accessor macros
#define HLS_CONTENT(S)  ( (HypreLSContent*)(S->content) )
#define HLS_X(S)        ( HLS_CONTENT(S)->xvec )
#define HLS_B(S)        ( HLS_CONTENT(S)->bvec )
#define HLS_PRECOND(S)  ( HLS_CONTENT(S)->precond )
#define HLS_SOLVER(S)   ( HLS_CONTENT(S)->solver )
#define HLS_ITERS(S)    ( HLS_CONTENT(S)->iters )
#define HLS_PCG(S)      ( HLS_CONTENT(S)->pcg )

// Solver function prototypes
SUNLinearSolver HypreLS(SUNMatrix A, UserData *udata);
SUNLinearSolver_Type HypreLS_GetType(SUNLinearSolver S);
int HypreLS_Initialize(SUNLinearSolver S);
int HypreLS_Setup(SUNLinearSolver S, SUNMatrix A);
int HypreLS_Solve(SUNLinearSolver S, SUNMatrix A,
                   N_Vector x, N_Vector b, realtype tol);
int HypreLS_NumIters(SUNLinearSolver S);
int HypreLS_Free(SUNLinearSolver S);

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// ODE right hand side function
static int f(realtype t, N_Vector u, N_Vector f, void *user_data);

// Jacobian evaluation function
static int Jac(realtype t, N_Vector u, N_Vector f, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Setup the parallel decomposition
static int SetupDecomp(MPI_Comm comm_w, UserData *udata);

// Perform neighbor exchange
static int PostRecv(UserData *udata);
static int SendData(N_Vector y, UserData *udata);
static int WaitRecv(UserData *udata);

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

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
  int flag; // reusable error-checking flag

  // Timing variables
  double t1 = 0.0;
  double t2 = 0.0;

  // MPI variables
  MPI_Comm comm_w = MPI_COMM_WORLD; // MPI communicator
  int myid;                         // MPI process ID

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;

  flag = MPI_Comm_rank(comm_w, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

  // Create a new scope so that sundials::Context is deleted
  // prior to the MPI_Finalize() call.
  {
    UserData *udata    = NULL;  // user data structure
    N_Vector u         = NULL;  // vector for storing solution
    SUNMatrix A        = NULL;  // matrix for Jacobian
    SUNLinearSolver LS = NULL;  // linear solver memory structure
    void *cvode_mem    = NULL;  // CVODE memory structure

    // SUNDIALS context
    sundials::Context sunctx(&comm_w);

    // Set output process flag
    bool outproc = (myid == 0);

    // ------------------------------------------
    // Setup UserData and parallel decomposition
    // ------------------------------------------

    // Allocate and initialize user data structure with default values. The
    // defaults may be overwritten by command line inputs in ReadInputs below.
    udata = new UserData(sunctx);

    // Parse command line inputs
    flag = ReadInputs(&argc, &argv, udata, outproc);
    if (flag != 0) return 1;

    // Setup parallel decomposition
    flag = SetupDecomp(comm_w, udata);
    if (check_flag(&flag, "SetupDecomp", 1)) return 1;

    // Output problem setup/options
    if (outproc)
    {
      flag = PrintUserData(udata);
      if (check_flag(&flag, "PrintUserData", 1)) return 1;
    }

    // ------------------------
    // Create parallel vectors
    // ------------------------

    // Create vector for solution
    u = N_VNew_Parallel(udata->comm_c, udata->nodes_loc, udata->nodes, sunctx);
    if (check_flag((void *) u, "N_VNew_Parallel", 0)) return 1;

    // Set initial condition
    flag = Solution(ZERO, u, udata);
    if (check_flag(&flag, "Solution", 1)) return 1;

    // Create vector for error
    udata->e = N_VClone(u);
    if (check_flag((void *) (udata->e), "N_VClone", 0)) return 1;

    // --------------------------------
    // Create matrix and linear solver
    // --------------------------------

    // Create custom matrix
    A = Hypre5ptMatrix(udata);
    if (check_flag((void *) A, "Hypre5ptMatrix", 0)) return 1;

    // Create linear solver
    LS = HypreLS(A, udata);
    if (check_flag((void *) LS, "HypreLS", 0)) return 1;

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
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1)) return 1;

    // Specify the Jacobian evaluation function
    flag = CVodeSetJacFn(cvode_mem, Jac);
    if (check_flag(&flag, "CVodeSetJacFn", 1)) return 1;

    // Set linear solver setup frequency (update linear system matrix)
    flag = CVodeSetLSetupFrequency(cvode_mem, udata->msbp);
    if (check_flag(&flag, "CVodeSetLSetupFrequency", 1)) return 1;

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
      t1 = MPI_Wtime();

      // Evolve in time
      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      if (check_flag(&flag, "CVode", 1)) break;

      // Stop timer
      t2 = MPI_Wtime();

      // Update timer
      udata->evolvetime += t2 - t1;

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
    if (udata->output > 0 && outproc)
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

      if (outproc)
      {
        cout << scientific;
        cout << setprecision(numeric_limits<realtype>::digits10);
        cout << "  Max error = " << maxerr << endl;
      }
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

    CVodeFree(&cvode_mem);     // Free integrator memory
    SUNLinSolFree(LS);         // Free linear solver
    SUNMatDestroy(A);          // Free matrix
    N_VDestroy(u);             // Free vectors
    delete udata;
  }

  flag = MPI_Finalize();     // Finalize MPI
  return 0;
}

// -----------------------------------------------------------------------------
// Setup the parallel decomposition
// -----------------------------------------------------------------------------

static int SetupDecomp(MPI_Comm comm_w, UserData *udata)
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
  flag = MPI_Comm_size(comm_w, &(udata->nprocs_w));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << flag << endl;
    return -1;
  }

  // Check the processor grid
  if ((udata->npx * udata->npy) != udata->nprocs_w)
  {
    cerr << "Error: npx * npy != nproc" << endl;
    return -1;
  }

  // Set up 2D Cartesian communicator
  int dims[2];
  dims[0] = udata->npx;
  dims[1] = udata->npy;

  int periods[2];
  periods[0] = 0;
  periods[1] = 0;

  flag = MPI_Cart_create(comm_w, 2, dims, periods, 0, &(udata->comm_c));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_create = " << flag << endl;
    return -1;
  }

  // Get my rank in the new Cartesian communicator
  flag = MPI_Comm_rank(udata->comm_c, &(udata->myid_c));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << flag << endl;
    return -1;
  }

  // Get dimension of the Cartesian communicator and my coordinates
  int coords[2];
  flag = MPI_Cart_get(udata->comm_c, 2, dims, periods, coords);
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
  udata->nodes     = udata->nx * udata->ny;
  udata->nodes_loc = udata->nx_loc * udata->ny_loc;

  // Determine if this proc has neighbors
  udata->HaveNbrW = (udata->is != 0);
  udata->HaveNbrE = (udata->ie != udata->nx-1);
  udata->HaveNbrS = (udata->js != 0);
  udata->HaveNbrN = (udata->je != udata->ny-1);

  // Allocate exchange buffers if necessary
  if (udata->HaveNbrW)
  {
    udata->Wrecv = new realtype[udata->ny_loc];
    udata->Wsend = new realtype[udata->ny_loc];
  }
  if (udata->HaveNbrE)
  {
    udata->Erecv = new realtype[udata->ny_loc];
    udata->Esend = new realtype[udata->ny_loc];
  }
  if (udata->HaveNbrS)
  {
    udata->Srecv = new realtype[udata->nx_loc];
    udata->Ssend = new realtype[udata->nx_loc];
  }
  if (udata->HaveNbrN)
  {
    udata->Nrecv = new realtype[udata->nx_loc];
    udata->Nsend = new realtype[udata->nx_loc];
  }

  // MPI neighborhood information
  int nbcoords[2];

  // West neighbor
  if (udata->HaveNbrW)
  {
    nbcoords[0] = coords[0]-1;
    nbcoords[1] = coords[1];
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // East neighbor
  if (udata->HaveNbrE)
  {
    nbcoords[0] = coords[0]+1;
    nbcoords[1] = coords[1];
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // South neighbor
  if (udata->HaveNbrS)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]-1;
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // North neighbor
  if (udata->HaveNbrN)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]+1;
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipN));
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
// Functions called by the integrator
// -----------------------------------------------------------------------------

// f routine to compute the ODE RHS function f(t,y).
static int f(realtype t, N_Vector u, N_Vector f, void *user_data)
{
  int          flag;
  sunindextype i, j;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Open exchange receives
  flag = PostRecv(udata);
  if (check_flag(&flag, "PostRecv", 1)) return -1;

  // Send exchange data
  flag = SendData(u, udata);
  if (check_flag(&flag, "SendData", 1)) return -1;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

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

  // Iterate over subdomain and compute rhs forcing term
  if (udata->forcing)
  {
    realtype x, y;
    realtype sin_sqr_x, sin_sqr_y;
    realtype cos_sqr_x, cos_sqr_y;

    realtype bx = (udata->kx) * TWO * PI * PI;
    realtype by = (udata->ky) * TWO * PI * PI;

    realtype sin_t_cos_t = sin(PI * t) * cos(PI * t);
    realtype cos_sqr_t   = cos(PI * t) * cos(PI * t);

    for (j = jstart; j < jend; j++)
    {
      for (i = istart; i < iend; i++)
      {
        x = (udata->is + i) * udata->dx;
        y = (udata->js + j) * udata->dy;

        sin_sqr_x = sin(PI * x) * sin(PI * x);
        sin_sqr_y = sin(PI * y) * sin(PI * y);

        cos_sqr_x = cos(PI * x) * cos(PI * x);
        cos_sqr_y = cos(PI * y) * cos(PI * y);

        farray[IDX(i,j,nx_loc)] =
          -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t
          -bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t
          -by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
      }
    }
  }

  // Iterate over subdomain interior and add rhs diffusion term
  for (j = 1; j < ny_loc - 1; j++)
  {
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }
  }

  // Wait for exchange receives
  flag = WaitRecv(udata);
  if (check_flag(&flag, "WaitRecv", 1)) return -1;

  // Iterate over subdomain boundaries and add rhs diffusion term
  realtype *Warray = udata->Wrecv;
  realtype *Earray = udata->Erecv;
  realtype *Sarray = udata->Srecv;
  realtype *Narray = udata->Nrecv;

  // West face (updates south-west and north-west corners if necessary)
  if (udata->HaveNbrW)
  {
    i = 0;
    if (udata->HaveNbrS)  // South-West corner
    {
      j = 0;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }

    if (udata->HaveNbrN)  // North-West corner
    {
      j = ny_loc - 1;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // East face (updates south-east and north-east corners if necessary)
  if (udata->HaveNbrE)
  {
    i = nx_loc - 1;
    if (udata->HaveNbrS)  // South-East corner
    {
      j = 0;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }

    if (udata->HaveNbrN)  // North-East corner
    {
      j = ny_loc - 1;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // South face (excludes corners)
  if (udata->HaveNbrS)
  {
    j = 0;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }
  }

  // North face (excludes corners)
  if (udata->HaveNbrN)
  {
    j = udata->ny_loc - 1;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->rhstime += t2 - t1;

  // Return success
  return 0;
}

// Jac function to compute the ODE RHS function Jacobian, (df/dy)(t,y).
static int Jac(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Shortcuts to hypre matrix and grid extents, work array, etc.
  HYPRE_StructMatrix Jmatrix = H5PM_MATRIX(J);

  UserData *udata = H5PM_UDATA(J);

  HYPRE_Int ilower[2];
  HYPRE_Int iupper[2];

  ilower[0] = H5PM_ILOWER(J)[0];
  ilower[1] = H5PM_ILOWER(J)[1];

  iupper[0] = H5PM_IUPPER(J)[0];
  iupper[1] = H5PM_IUPPER(J)[1];

  HYPRE_Int   nwork = H5PM_NWORK(J);
  HYPRE_Real *work  = H5PM_WORK(J);

  sunindextype nx_loc = iupper[0] - ilower[0] + 1;
  sunindextype ny_loc = iupper[1] - ilower[1] + 1;

  // Matrix stencil: center, left, right, bottom, top
  HYPRE_Int entries[5] = {0, 1, 2, 3, 4};
  HYPRE_Int entry[1];

  // Grid extents for setting boundary entries
  HYPRE_Int bc_ilower[2];
  HYPRE_Int bc_iupper[2];

  // Loop counters
  HYPRE_Int idx, ix, iy;

  // hypre return flag
  int flag;

  // ----------
  // Compute J
  // ----------

  // Start timer
  double t1 = MPI_Wtime();

  // Only do work if the box is non-zero in size
  if ((ilower[0] <= iupper[0]) &&
      (ilower[1] <= iupper[1]))
  {
    // Jacobian values
    realtype cx = udata->kx / (udata->dx * udata->dx);
    realtype cy = udata->ky / (udata->dy * udata->dy);
    realtype cc = -TWO * (cx + cy);

    // --------------------------------
    // Set matrix values for all nodes
    // --------------------------------

    // Set the matrix interior entries (center, left, right, bottom, top)
    idx = 0;
    for (iy = 0; iy < ny_loc; iy++)
    {
      for (ix = 0; ix < nx_loc; ix++)
      {
        work[idx]     = cc;
        work[idx + 1] = cx;
        work[idx + 2] = cx;
        work[idx + 3] = cy;
        work[idx + 4] = cy;
        idx += 5;
      }
    }

    // Modify the matrix
    flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                          ilower, iupper,
                                          5, entries, work);
    if (flag != 0) return -1;

    // ----------------------------------------
    // Correct matrix values at boundary nodes
    // ----------------------------------------

    // Set the matrix boundary entries (center, left, right, bottom, top)
    if (ilower[1] == 0 ||
        iupper[1] == (udata->ny - 1) ||
        ilower[0] == 0 ||
        iupper[0] == (udata->nx - 1))
    {
      idx = 0;
      for (iy = 0; iy < ny_loc; iy++)
      {
        for (ix = 0; ix < nx_loc; ix++)
        {
          work[idx]     = ONE;
          work[idx + 1] = ZERO;
          work[idx + 2] = ZERO;
          work[idx + 3] = ZERO;
          work[idx + 4] = ZERO;
          idx += 5;
        }
      }
    }

    // Set cells on western boundary
    if (ilower[0] == 0)
    {
      // Grid cell on south-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on north-west corner
      bc_iupper[0] = ilower[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                              bc_ilower, bc_iupper,
                                              5, entries, work);
        if (flag != 0) return -1;
      }
    }

    // Set cells on eastern boundary
    if (iupper[0] == (udata->nx - 1))
    {
      // Grid cell on south-east corner
      bc_ilower[0] = iupper[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                              bc_ilower, bc_iupper,
                                              5, entries, work);
        if (flag != 0) return -1;
      }
    }

    // Correct cells on southern boundary
    if (ilower[1] == 0)
    {
      // Grid cell on south-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on south-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = ilower[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                              bc_ilower, bc_iupper,
                                              5, entries, work);
        if (flag != 0) return -1;
      }
    }

    // Set cells on northern boundary
    if (iupper[1] == (udata->ny - 1))
    {
      // Grid cell on north-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = iupper[1];

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                              bc_ilower, bc_iupper,
                                              5, entries, work);
        if (flag != 0) return -1;
      }
    }

    // -----------------------------------------------------------
    // Remove connections between the interior and boundary nodes
    // -----------------------------------------------------------

    // Zero out work array
    for (ix = 0; ix < nwork; ix++)
    {
      work[ix] = ZERO;
    }

    // Second column of nodes (depends on western boundary)
    if ((ilower[0] <= 1) && (iupper[0] >= 1))
    {
      // Remove western dependency
      entry[0] = 1;

      // Grid cell on south-west corner
      bc_ilower[0] = 1;
      bc_ilower[1] = ilower[1];

      // Grid cell on north-west corner
      bc_iupper[0] = 1;
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                              bc_ilower, bc_iupper,
                                              1, entry, work);
        if (flag != 0) return -1;
      }
    }

    // Next to last column (depends on eastern boundary)
    if ((ilower[0] <= (udata->nx - 2)) &&
        (iupper[0] >= (udata->nx - 2)))
    {
      // Remove eastern dependency
      entry[0] = 2;

      // Grid cell on south-east corner
      bc_ilower[0] = udata->nx - 2;
      bc_ilower[1] = ilower[1];

      // Grid cell on north-east corner
      bc_iupper[0] = udata->nx - 2;
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                              bc_ilower, bc_iupper,
                                              1, entry, work);
        if (flag != 0) return -1;
      }
    }

    // Second row of nodes (depends on southern boundary)
    if ((ilower[1] <= 1) && (iupper[1] >= 1))
    {
      // Remove southern dependency
      entry[0] = 3;

      // Grid cell on south-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = 1;

      // Grid cell on south-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = 1;

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                              bc_ilower, bc_iupper,
                                              1, entry, work);
        if (flag != 0) return -1;
      }
    }

    // Next to last row of nodes (depends on northern boundary)
    if ((ilower[1] <= (udata->ny - 2)) &&
        (iupper[1] >= (udata->ny - 2)))
    {
      // Remove northern dependency
      entry[0] = 4;

      // Grid cell on north-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = udata->ny - 2;

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = udata->ny - 2;

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        flag = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                              bc_ilower, bc_iupper,
                                              1, entry, work);
        if (flag != 0) return -1;
      }
    }
  }

  // The matrix is assembled matrix in linear solver setup

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->matfilltime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// RHS helper functions
// -----------------------------------------------------------------------------

// Post exchange receives
static int PostRecv(UserData *udata)
{
  int flag;

  // Start timer
  double t1 = MPI_Wtime();

  // Open Irecv buffers
  if (udata->HaveNbrW)
  {
    flag = MPI_Irecv(udata->Wrecv, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipW, MPI_ANY_TAG, udata->comm_c, &(udata->reqRW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    flag = MPI_Irecv(udata->Erecv, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipE, MPI_ANY_TAG, udata->comm_c, &(udata->reqRE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    flag = MPI_Irecv(udata->Srecv, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipS, MPI_ANY_TAG, udata->comm_c, &(udata->reqRS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    flag = MPI_Irecv(udata->Nrecv, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipN, MPI_ANY_TAG, udata->comm_c, &(udata->reqRN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Send exchange data
static int SendData(N_Vector y, UserData *udata)
{
  int flag, i;
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Start timer
  double t1 = MPI_Wtime();

  // Access data array
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;

  // Send data
  if (udata->HaveNbrW)
  {
    for (i = 0; i < ny_loc; i++) udata->Wsend[i] = Y[IDX(0,i,nx_loc)];
    flag = MPI_Isend(udata->Wsend, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipW, 0, udata->comm_c, &(udata->reqSW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    for (i = 0; i < ny_loc; i++) udata->Esend[i] = Y[IDX(nx_loc-1,i,nx_loc)];
    flag = MPI_Isend(udata->Esend, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipE, 1, udata->comm_c, &(udata->reqSE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    for (i = 0; i < nx_loc; i++) udata->Ssend[i] = Y[IDX(i,0,nx_loc)];
    flag = MPI_Isend(udata->Ssend, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipS, 2, udata->comm_c, &(udata->reqSS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    for (i = 0; i < nx_loc; i++) udata->Nsend[i] = Y[IDX(i,ny_loc-1,nx_loc)];
    flag = MPI_Isend(udata->Nsend, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipN, 3, udata->comm_c, &(udata->reqSN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Wait for exchange data
static int WaitRecv(UserData *udata)
{
  // Local variables
  int flag;
  MPI_Status stat;

  // Start timer
  double t1 = MPI_Wtime();

  // Wait for messages to finish
  if (udata->HaveNbrW)
  {
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
  }

  if (udata->HaveNbrE)
  {
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
  }

  if (udata->HaveNbrS)
  {
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
  }

  if (udata->HaveNbrN)
  {
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
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Initialize memory allocated within Userdata
UserData::UserData(sundials::Context& sunctx)
  : sunctx(sunctx)
{
  // Diffusion coefficient
  kx = ONE;
  ky = ONE;

  // Enable forcing
  forcing = true;

  // Final time
  tf = ONE;

  // Upper bounds in x and y directions
  xu = ONE;
  yu = ONE;

  // Global number of nodes in the x and y directions
  nx    = 64;
  ny    = 64;
  nodes = nx * ny;

  // Mesh spacing in the x and y directions
  dx = xu / (nx - 1);
  dy = yu / (ny - 1);

  // Locals number of nodes in the x and y directions (set in SetupDecomp)
  nx_loc    = 0;
  ny_loc    = 0;
  nodes_loc = 0;

  // Global indices of this subdomain (set in SetupDecomp)
  is = 0;
  ie = 0;
  js = 0;
  je = 0;

  // MPI variables (set in SetupDecomp)
  comm_c = MPI_COMM_NULL;

  nprocs_w = 1;
  npx      = 1;
  npy      = 1;

  myid_c = 0;

  // Flags denoting neighbors (set in SetupDecomp)
  HaveNbrW = true;
  HaveNbrE = true;
  HaveNbrS = true;
  HaveNbrN = true;

  // Exchange receive buffers (allocated in SetupDecomp)
  Erecv = NULL;
  Wrecv = NULL;
  Nrecv = NULL;
  Srecv = NULL;

  // Exchange send buffers (allocated in SetupDecomp)
  Esend = NULL;
  Wsend = NULL;
  Nsend = NULL;
  Ssend = NULL;

  // Neighbors IDs (set in SetupDecomp)
  ipW = -1;
  ipE = -1;
  ipS = -1;
  ipN = -1;

  // Integrator settings
  rtol     = RCONST(1.e-5);   // relative tolerance
  atol     = RCONST(1.e-10);  // absolute tolerance
  maxsteps = 0;               // use default

  // Linear solver and preconditioner options
  pcg       = true;       // use PCG (true) or GMRES (false)
  prec      = true;       // enable preconditioning
  liniters  = 10;         // max linear iterations
  msbp      = 0;          // use default (20 steps)
  epslin    = ZERO;       // use default (0.05)

  // hypre PFMG settings
  pfmg_relax  = 2;
  pfmg_nrelax = 2;

  // Output variables
  output = 1;   // 0 = no output, 1 = stats output, 2 = output to disk
  nout   = 20;  // Number of output times
  e      = NULL;

  // Timing variables
  timing       = false;
  evolvetime   = 0.0;
  rhstime      = 0.0;
  matfilltime  = 0.0;
  setuptime    = 0.0;
  solvetime    = 0.0;
  exchangetime = 0.0;
}

// Free memory allocated within Userdata
UserData::~UserData()
{
  // Free exchange buffers
  if (Wrecv != NULL)  delete[] Wrecv;
  if (Wsend != NULL)  delete[] Wsend;
  if (Erecv != NULL)  delete[] Erecv;
  if (Esend != NULL)  delete[] Esend;
  if (Srecv != NULL)  delete[] Srecv;
  if (Ssend != NULL)  delete[] Ssend;
  if (Nrecv != NULL)  delete[] Nrecv;
  if (Nsend != NULL)  delete[] Nsend;

  // Free MPI Cartesian communicator
  if (comm_c != MPI_COMM_NULL)
    MPI_Comm_free(&(comm_c));

  // Free error vector
  if (e)
  {
    N_VDestroy(e);
    e = NULL;
  }
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
    // MPI processes
    else if (arg == "--np")
    {
      udata->npx = stoi((*argv)[arg_idx++]);
      udata->npy = stoi((*argv)[arg_idx++]);
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
    // PFMG settings
    else if (arg == "--pfmg_relax")
    {
      udata->pfmg_relax = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--pfmg_nrelax")
    {
      udata->pfmg_nrelax = stoi((*argv)[arg_idx++]);
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

  // Iterative over domain interior
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? udata->nx_loc : udata->nx_loc - 1;

  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? udata->ny_loc : udata->ny_loc - 1;

  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  for (sunindextype j = jstart; j < jend; j++)
  {
    for (sunindextype i = istart; i < iend; i++)
    {
      x  = (udata->is + i) * udata->dx;
      y  = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uarray[IDX(i,j,udata->nx_loc)] = sin_sqr_x * sin_sqr_y * cos_sqr_t;
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
  cout << "  --np <npx> <npy>        : number of MPI processes in the x and y directions" << endl;
  cout << "  --domain <xu> <yu>      : domain upper bound in the x and y direction" << endl;
  cout << "  --k <kx> <ky>           : diffusion coefficients" << endl;
  cout << "  --noforcing             : disable forcing term" << endl;
  cout << "  --tf <time>             : final time" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --atol <atol>           : absoltue tolerance" << endl;
  cout << "  --gmres                 : use GMRES linear solver" << endl;
  cout << "  --liniters <iters>      : max number of iterations" << endl;
  cout << "  --epslin <factor>       : linear tolerance factor" << endl;
  cout << "  --noprec                : disable preconditioner" << endl;
  cout << "  --msbp <steps>          : max steps between prec setups" << endl;
  cout << "  --pfmg_relax <types>    : relaxtion type in PFMG" << endl;
  cout << "  --pfmg_nrelax <iters>   : pre/post relaxtion sweeps in PFMG" << endl;
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
  cout << "  nprocs         = " << udata->nprocs_w        << endl;
  cout << "  npx            = " << udata->npx             << endl;
  cout << "  npy            = " << udata->npy             << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  kx             = " << udata->kx              << endl;
  cout << "  ky             = " << udata->ky              << endl;
  cout << "  forcing        = " << udata->forcing         << endl;
  cout << "  tf             = " << udata->tf              << endl;
  cout << "  xu             = " << udata->xu              << endl;
  cout << "  yu             = " << udata->yu              << endl;
  cout << "  nx             = " << udata->nx              << endl;
  cout << "  ny             = " << udata->ny              << endl;
  cout << "  nxl (proc 0)   = " << udata->nx_loc          << endl;
  cout << "  nyl (proc 0)   = " << udata->ny_loc          << endl;
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
  cout << "  pfmg_relax     = " << udata->pfmg_relax      << endl;
  cout << "  pfmg_nrelax    = " << udata->pfmg_nrelax     << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  output         = " << udata->output          << endl;
  cout << " --------------------------------- "           << endl;
  cout << endl;

  return 0;
}

// Initialize output
static int OpenOutput(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  // Header for status output
  if (udata->output > 0 && outproc)
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
    stringstream fname;
    fname << "heat2d_info." << setfill('0') << setw(5) << udata->myid_c
          << ".txt";

    ofstream dout;
    dout.open(fname.str());
    dout <<  "xu  " << udata->xu       << endl;
    dout <<  "yu  " << udata->yu       << endl;
    dout <<  "nx  " << udata->nx       << endl;
    dout <<  "ny  " << udata->ny       << endl;
    dout <<  "px  " << udata->npx      << endl;
    dout <<  "py  " << udata->npy      << endl;
    dout <<  "np  " << udata->nprocs_w << endl;
    dout <<  "is  " << udata->is       << endl;
    dout <<  "ie  " << udata->ie       << endl;
    dout <<  "js  " << udata->js       << endl;
    dout <<  "je  " << udata->je       << endl;
    dout <<  "nt  " << udata->nout + 1 << endl;
    dout.close();

    // Open output streams for solution and error
    fname.str("");
    fname.clear();
    fname << "heat2d_solution." << setfill('0') << setw(5) << udata->myid_c
          << ".txt";
    udata->uout.open(fname.str());

    udata->uout << scientific;
    udata->uout << setprecision(numeric_limits<realtype>::digits10);

    if (udata->forcing)
    {
      fname.str("");
      fname.clear();
      fname << "heat2d_error." << setfill('0') << setw(5) << udata->myid_c
            << ".txt";
      udata->eout.open(fname.str());

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
  bool     outproc = (udata->myid_c == 0);

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
    if (outproc)
    {
      if (udata->forcing)
      {
        cout << setw(22) << t << setw(25) << urms << setw(25) << max << endl;
      }
      else
      {
        cout << setw(22) << t << setw(25) << urms << endl;
      }
    }

    // Write solution and error to disk
    if (udata->output == 2)
    {
      realtype *uarray = N_VGetArrayPointer(u);
      if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

      udata->uout << t << " ";
      for (sunindextype i = 0; i < udata->nodes_loc; i++)
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
        for (sunindextype i = 0; i < udata->nodes_loc; i++)
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
  bool outproc = (udata->myid_c == 0);

  // Footer for status output
  if (outproc && (udata->output > 0))
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
  long int nst, netf, nf, nni, ncfn, nli, nlcf, nsetups, nJeval;
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
  flag = CVodeGetNumJacEvals(cvode_mem, &nJeval);
  if (check_flag(&flag, "CVodeGetNumJacEvals", 1)) return -1;

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
  cout << "  J evals          = " << nJeval  << endl;
  cout << endl;

  // Compute average nls iters per step attempt and ls iters per nls iter
  realtype avgnli = (realtype) nni / (realtype) nst;
  realtype avgli  = (realtype) nli / (realtype) nni;
  cout << "  Avg NLS iters per step    = " << avgnli << endl;
  cout << "  Avg LS iters per NLS iter = " << avgli  << endl;
  cout << endl;

  return 0;
}

static int OutputTiming(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double maxtime = 0.0;

  MPI_Reduce(&(udata->evolvetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Evolve time   = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->rhstime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  RHS time      = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->exchangetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Exchange time = " << maxtime << " sec" << endl;
    cout << endl;
  }

  MPI_Reduce(&(udata->matfilltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  MatFill time  = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->setuptime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  LS setup time = " << maxtime << " sec" << endl;   }

  MPI_Reduce(&(udata->solvetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  LS solve time = " << maxtime << " sec" << endl;
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

// -----------------------------------------------------------------------------
// SUNMatrix functions
// -----------------------------------------------------------------------------

SUNMatrix Hypre5ptMatrix(UserData *udata)
{
  int flag, result;

  // Check input
  if (udata == NULL) return NULL;

  // Check for valid 2D Cartesian MPI communicator
  flag = MPI_Topo_test(udata->comm_c, &result);
  if ((flag != MPI_SUCCESS) || (result != MPI_CART))  return NULL;

  flag = MPI_Cartdim_get(udata->comm_c, &result);
  if ((flag != MPI_SUCCESS) || (result != 2))  return NULL;

  // Create an empty matrix object
  SUNMatrix A = SUNMatNewEmpty(udata->sunctx);
  if (A == NULL) return NULL;

  // Attach operations
  A->ops->getid     = Hypre5ptMatrix_GetID;
  A->ops->clone     = Hypre5ptMatrix_Clone;
  A->ops->destroy   = Hypre5ptMatrix_Destroy;
  A->ops->copy      = Hypre5ptMatrix_Copy;
  A->ops->scaleaddi = Hypre5ptMatrix_ScaleAddI;

  // Create content
  Hypre5ptMatrixContent *content = new Hypre5ptMatrixContent;
  if (content == NULL) { SUNMatDestroy(A); return NULL; }

  // Attach content
  A->content = content;

  // User data
  content->udata = udata;

  // -----
  // Grid
  // -----

  // Create 2D grid object
  flag = HYPRE_StructGridCreate(udata->comm_c, 2, &(content->grid));
  if (flag != 0) { SUNMatDestroy(A); return NULL; }

  // Set grid extents (lower left and upper right corners)
  content->ilower[0] = udata->is;
  content->ilower[1] = udata->js;

  content->iupper[0] = udata->ie;
  content->iupper[1] = udata->je;

  flag = HYPRE_StructGridSetExtents(content->grid,
                                    content->ilower, content->iupper);
  if (flag != 0) { SUNMatDestroy(A); return NULL; }

  // Assemble the grid
  flag = HYPRE_StructGridAssemble(content->grid);
  if (flag != 0) { SUNMatDestroy(A); return NULL; }

  // --------
  // Stencil
  // --------

  // Create the 2D 5 point stencil object
  flag = HYPRE_StructStencilCreate(2, 5, &(content->stencil));
  if (flag != 0) { SUNMatDestroy(A); return NULL; }

  // Set the stencil entries (center, left, right, bottom, top)
  HYPRE_Int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

  for (int entry = 0; entry < 5; entry++)
  {
    flag = HYPRE_StructStencilSetElement(content->stencil, entry,
                                         offsets[entry]);
    if (flag != 0) { SUNMatDestroy(A); return NULL; }
  }

  // -----------
  // Work array
  // -----------

  content->nwork = 5 * (udata->nodes_loc);
  content->work  = NULL;
  content->work  = new HYPRE_Real[content->nwork];
  if (flag != 0) { SUNMatDestroy(A); return NULL; }

  // ---------
  // A matrix
  // ---------

  flag = HYPRE_StructMatrixCreate(udata->comm_c, content->grid,
                                  content->stencil,
                                  &(content->matrix));
  if (flag != 0) { SUNMatDestroy(A); return NULL; }

  flag = HYPRE_StructMatrixInitialize(content->matrix);
  if (flag != 0) { SUNMatDestroy(A); return NULL; }

  // Return matrix
  return(A);
}

SUNMatrix_ID Hypre5ptMatrix_GetID(SUNMatrix A)
{
  return SUNMATRIX_CUSTOM;
}

SUNMatrix Hypre5ptMatrix_Clone(SUNMatrix A)
{
  return(Hypre5ptMatrix(H5PM_UDATA(A)));
}

void Hypre5ptMatrix_Destroy(SUNMatrix A)
{
  if (A == NULL) return;
  if (A->content != NULL)
  {
    if (H5PM_WORK(A))    delete[] H5PM_WORK(A);
    if (H5PM_MATRIX(A))  HYPRE_StructMatrixDestroy(H5PM_MATRIX(A));
    if (H5PM_GRID(A))    HYPRE_StructGridDestroy(H5PM_GRID(A));
    if (H5PM_STENCIL(A)) HYPRE_StructStencilDestroy(H5PM_STENCIL(A));
    delete ((Hypre5ptMatrixContent*) (A->content));
    A->content = NULL;
  }
  SUNMatFreeEmpty(A);
  return;
}

int Hypre5ptMatrix_Copy(SUNMatrix A, SUNMatrix B)
{
  int flag;

  // Center, left, right, bottom, top
  HYPRE_Int entries[5] = {0, 1, 2, 3, 4};

  // Copy values from A into work array
  flag = HYPRE_StructMatrixGetBoxValues(H5PM_MATRIX(A),
                                        H5PM_ILOWER(A), H5PM_IUPPER(A),
                                        5, entries, H5PM_WORK(A));
  if (flag != 0) return(flag);

  // Insert values into B
  flag = HYPRE_StructMatrixSetBoxValues(H5PM_MATRIX(B),
                                        H5PM_ILOWER(A), H5PM_IUPPER(A),
                                        5, entries, H5PM_WORK(A));
  if (flag != 0) return(flag);

  // Return success
  return(SUNMAT_SUCCESS);
}

int Hypre5ptMatrix_ScaleAddI(realtype c, SUNMatrix A)
{
  int flag;

  // Center, left, right, bottom, top
  HYPRE_Int entries[5] = {0, 1, 2, 3, 4};

  // Copy all matrix values into work array
  flag = HYPRE_StructMatrixGetBoxValues(H5PM_MATRIX(A),
                                        H5PM_ILOWER(A), H5PM_IUPPER(A),
                                        5, entries, H5PM_WORK(A));
  if (flag != 0) return(flag);

  // Scale work array by c
  for (HYPRE_Int i = 0; i < H5PM_NWORK(A); i++)
    H5PM_WORK(A)[i] *= c;

  // Insert scaled values back into A
  flag = HYPRE_StructMatrixSetBoxValues(H5PM_MATRIX(A),
                                        H5PM_ILOWER(A), H5PM_IUPPER(A),
                                        5, entries, H5PM_WORK(A));
  if (flag != 0) return(flag);

  // Set first 1/5 of work array to 1
  for (HYPRE_Int i = 0; i < H5PM_NWORK(A)/5; i++)
    H5PM_WORK(A)[i] = ONE;

  // Insert resulting values back into diagonal of A
  HYPRE_Int entry[1] = {0};
  flag = HYPRE_StructMatrixAddToBoxValues(H5PM_MATRIX(A),
                                          H5PM_ILOWER(A), H5PM_IUPPER(A),
                                          1, entry, H5PM_WORK(A));
  if (flag != 0) return(flag);

  // Return success
  return(SUNMAT_SUCCESS);
}

// -----------------------------------------------------------------------------
// SUNLinearSolver functions
// -----------------------------------------------------------------------------

// Create hypre linear solver
SUNLinearSolver HypreLS(SUNMatrix A, UserData *udata)
{
  int flag;

  // Check input
  if (udata == NULL) return NULL;

  // Create an empty linear solver
  SUNLinearSolver LS = SUNLinSolNewEmpty(udata->sunctx);
  if (LS == NULL) return NULL;

  // Attach operations
  LS->ops->gettype    = HypreLS_GetType;
  LS->ops->setup      = HypreLS_Setup;
  LS->ops->solve      = HypreLS_Solve;
  LS->ops->numiters   = HypreLS_NumIters;
  LS->ops->free       = HypreLS_Free;

  // Create content
  HypreLSContent *content = new HypreLSContent;
  if (content == NULL) { SUNLinSolFree(LS); return NULL; }

  // Attach the content
  LS->content = content;

  // Initialize iteration count and solver choice
  content->iters = 0;
  content->pcg   = udata->pcg;

  // ---------
  // x vector
  // ---------

  flag = HYPRE_StructVectorCreate(udata->comm_c, H5PM_GRID(A),
                                  &(content->xvec));
  if (flag != 0) { SUNLinSolFree(LS); return NULL; }

  flag = HYPRE_StructVectorInitialize(content->xvec);
  if (flag != 0) { SUNLinSolFree(LS); return NULL; }

  // ---------
  // b vector
  // ---------

  flag = HYPRE_StructVectorCreate(udata->comm_c, H5PM_GRID(A),
                                  &(content->bvec));
  if (flag != 0) { SUNLinSolFree(LS); return NULL; }

  flag = HYPRE_StructVectorInitialize(content->bvec);
  if (flag != 0) { SUNLinSolFree(LS); return NULL; }

  // --------------
  // linear solver
  // --------------

  if (content->pcg)
  {
    // Create the struct PCG solver
    flag = HYPRE_StructPCGCreate(udata->comm_c, &(content->solver));
    if (flag != 0) { HypreLS_Free(LS); return(NULL); }

    // Max number of iterations
    flag = HYPRE_StructPCGSetMaxIter(content->solver, udata->liniters);
    if (flag != 0) { HypreLS_Free(LS); return(NULL); }
  }
  else
  {
    // Create the struct GMRES solver
    flag = HYPRE_StructGMRESCreate(udata->comm_c, &(content->solver));
    if (flag != 0) { HypreLS_Free(LS); return(NULL); }

    // Max Krylov space size and number of iterations (no restarts)
    flag = HYPRE_StructGMRESSetKDim(content->solver, udata->liniters);
    if (flag != 0) { HypreLS_Free(LS); return(NULL); }

    flag = HYPRE_StructGMRESSetMaxIter(content->solver, udata->liniters);
    if (flag != 0) { HypreLS_Free(LS); return(NULL); }
  }

  // --------------------
  // PFMG preconditioner
  // --------------------

  // Note a new PFMG preconditioner must be created and attached each time the
  // linear system is updated. As such it is constructed in the linear solver
  // setup function.
  content->precond = NULL;

  // Return solver
  return(LS);
}

SUNLinearSolver_Type HypreLS_GetType(SUNLinearSolver LS)
{
  return(SUNLINEARSOLVER_MATRIX_ITERATIVE);
}

int HypreLS_Setup(SUNLinearSolver LS, SUNMatrix A)
{
  int flag;

  // Start timer
  double t1 = MPI_Wtime();

  // Shortcut to user data
  UserData *udata = H5PM_UDATA(A);

  // Assemble the matrix
  flag = HYPRE_StructMatrixAssemble(H5PM_MATRIX(A));
  if (flag != 0) return(flag);

  // Set rhs/solution vectors as all zero for now
  flag = HYPRE_StructVectorSetConstantValues(HLS_B(LS), ZERO);
  if (flag != 0) return(flag);

  flag = HYPRE_StructVectorAssemble(HLS_B(LS));
  if (flag != 0) return(flag);

  flag = HYPRE_StructVectorSetConstantValues(HLS_X(LS), ZERO);
  if (flag != 0) return(flag);

  flag = HYPRE_StructVectorAssemble(HLS_X(LS));
  if (flag != 0) return(flag);

  // Setup the preconditioner
  if (udata->prec)
  {
    // Free the existing preconditioner if necessary
    if (HLS_PRECOND(LS)) HYPRE_StructPFMGDestroy(HLS_PRECOND(LS));

    // Create the new preconditioner
    flag = HYPRE_StructPFMGCreate(udata->comm_c, &(HLS_PRECOND(LS)));
    if (flag != 0) return(flag);

    // Signal that the inital guess is zero
    flag = HYPRE_StructPFMGSetZeroGuess(HLS_PRECOND(LS));
    if (flag != 0) return(flag);

    // tol <= 0.0 means do the max number of iterations
    flag = HYPRE_StructPFMGSetTol(HLS_PRECOND(LS), ZERO);
    if (flag != 0) return(flag);

    // Use one v-cycle
    flag = HYPRE_StructPFMGSetMaxIter(HLS_PRECOND(LS), 1);
    if (flag != 0) return(flag);

    // Use non-Galerkin corase grid operator
    flag = HYPRE_StructPFMGSetRAPType(HLS_PRECOND(LS), 1);
    if (flag != 0) return(flag);

    // Set the relaxation type
    flag = HYPRE_StructPFMGSetRelaxType(HLS_PRECOND(LS), udata->pfmg_relax);
    if (flag != 0) return(flag);

    // Set the number of pre and post relaxation sweeps
    flag = HYPRE_StructPFMGSetNumPreRelax(HLS_PRECOND(LS), udata->pfmg_nrelax);
    if (flag != 0) return(flag);

    flag = HYPRE_StructPFMGSetNumPostRelax(HLS_PRECOND(LS), udata->pfmg_nrelax);
    if (flag != 0) return(flag);

    // Set preconditioner
    if (HLS_PCG(LS))
    {
      flag = HYPRE_StructPCGSetPrecond(HLS_SOLVER(LS),
                                       HYPRE_StructPFMGSolve,
                                       HYPRE_StructPFMGSetup,
                                       HLS_PRECOND(LS));
    }
    else
    {
      flag = HYPRE_StructGMRESSetPrecond(HLS_SOLVER(LS),
                                         HYPRE_StructPFMGSolve,
                                         HYPRE_StructPFMGSetup,
                                         HLS_PRECOND(LS));
    }
    if (flag != 0) return(flag);
  }

  // Set up the solver
  if (HLS_PCG(LS))
  {
    flag = HYPRE_StructPCGSetup(HLS_SOLVER(LS), H5PM_MATRIX(A),
                                HLS_B(LS), HLS_X(LS));
  }
  else
  {
    flag = HYPRE_StructGMRESSetup(HLS_SOLVER(LS), H5PM_MATRIX(A),
                                  HLS_B(LS), HLS_X(LS));
  }
  if (flag != 0) return(flag);

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->setuptime += t2 - t1;

  // Return success
  return(SUNLS_SUCCESS);
}

int HypreLS_Solve(SUNLinearSolver LS, SUNMatrix A,
                  N_Vector x, N_Vector b, realtype tol)
{
  int flag;

  // Start timer
  double t1 = MPI_Wtime();

  // Shortcut to user data
  UserData *udata = H5PM_UDATA(A);

  // Insert rhs N_Vector entries into HYPRE vector b and assemble
  flag = HYPRE_StructVectorSetBoxValues(HLS_B(LS),
                                        H5PM_ILOWER(A), H5PM_IUPPER(A),
                                        N_VGetArrayPointer(b));
  if (flag != 0) return -1;

  flag = HYPRE_StructVectorAssemble(HLS_B(LS));
  if (flag != 0) return -1;

  // Insert solution N_Vector entries into HYPRE vector x and assemble
  flag = HYPRE_StructVectorSetBoxValues(HLS_X(LS),
                                        H5PM_ILOWER(A), H5PM_IUPPER(A),
                                        N_VGetArrayPointer(x));
  if (flag != 0) return -1;

  flag = HYPRE_StructVectorAssemble(HLS_X(LS));
  if (flag != 0) return -1;

  if (HLS_PCG(LS))
  {
    // Relative tolerance
    flag = HYPRE_StructPCGSetTol(HLS_SOLVER(LS), ZERO);
    if (flag != 0) return -1;

    // Absolute tolerance
    flag = HYPRE_StructPCGSetAbsoluteTol(HLS_SOLVER(LS), tol);
    if (flag != 0) return -1;

    // Use two norm
    flag = HYPRE_StructPCGSetTwoNorm(HLS_SOLVER(LS), 1);
    if (flag != 0) return -1;

    // Solve the linear system
    flag = HYPRE_StructPCGSolve(HLS_SOLVER(LS), H5PM_MATRIX(A),
                                HLS_B(LS), HLS_X(LS));
  }
  else
  {
    // Relative tolerance
    flag = HYPRE_StructGMRESSetTol(HLS_SOLVER(LS), ZERO);
    if (flag != 0) return -1;

    // Absolute tolerance
    flag = HYPRE_StructGMRESSetAbsoluteTol(HLS_SOLVER(LS), tol);
    if (flag != 0) return -1;

    // Solve the linear system
    flag = HYPRE_StructGMRESSolve(HLS_SOLVER(LS), H5PM_MATRIX(A),
                                  HLS_B(LS), HLS_X(LS));
  }

  // If a convergence error occured, clear the error, and return with a
  // recoverable error.
  if (flag == HYPRE_ERROR_CONV)
  {
    HYPRE_ClearError(HYPRE_ERROR_CONV);
    return SUNLS_CONV_FAIL;
  }
  // If any other error occured return with an unrecoverable error.
  else if (flag != 0)
  {
    return SUNLS_PACKAGE_FAIL_UNREC;
  }

  // Update iteration count
  if (HLS_PCG(LS))
  {
    flag = HYPRE_StructPCGGetNumIterations(HLS_SOLVER(LS), &(HLS_ITERS(LS)));
    if (flag != 0) return -1;
  }
  else
  {
    flag = HYPRE_StructGMRESGetNumIterations(HLS_SOLVER(LS), &(HLS_ITERS(LS)));
    if (flag != 0) return -1;
  }

  // Extract solution values
  flag = HYPRE_StructVectorGetBoxValues(HLS_X(LS),
                                        H5PM_ILOWER(A), H5PM_IUPPER(A),
                                        N_VGetArrayPointer(x));
  if (flag != 0) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->solvetime += t2 - t1;

  // Return success
  return(SUNLS_SUCCESS);
}

int HypreLS_NumIters(SUNLinearSolver LS)
{
  return((int) HLS_ITERS(LS));
}

int HypreLS_Free(SUNLinearSolver LS)
{
  if (LS == NULL) return(SUNLS_SUCCESS);
  if (LS->content != NULL)
  {
    if (HLS_SOLVER(LS))
    {
      if (HLS_PCG(LS))
      {
        HYPRE_StructPCGDestroy(HLS_SOLVER(LS));
      }
      else
      {
        HYPRE_StructGMRESDestroy(HLS_SOLVER(LS));
      }
    }
    if (HLS_PRECOND(LS)) HYPRE_StructPFMGDestroy(HLS_PRECOND(LS));
    if (HLS_B(LS))       HYPRE_StructVectorDestroy(HLS_B(LS));
    if (HLS_X(LS))       HYPRE_StructVectorDestroy(HLS_X(LS));
    delete ((HypreLSContent*) (LS->content));
    LS->content = NULL;
  }
  SUNLinSolFreeEmpty(LS);
  return(SUNLS_SUCCESS);
}

//---- end of file ----
