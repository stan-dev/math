/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 LLNS/SMU Copyright Start
 Copyright (c) 2018, Southern Methodist University and
 Lawrence Livermore National Security
 
 This work was performed under the auspices of the U.S. Department
 of Energy by Southern Methodist University and Lawrence Livermore
 National Laboratory under Contract DE-AC52-07NA27344.
 Produced at Southern Methodist University and the Lawrence
 Livermore National Laboratory.
 
 All rights reserved.
 For details, see the LICENSE file.
 LLNS/SMU Copyright End
 ----------------------------------------------------------------
 Example problem:

 The following test simulates a simple anisotropic 2D heat
 equation,
    u_t = kx*u_xx + ky*u_yy + h,
 for t in [0, 10], (x,y) in [0, 1]^2, with initial conditions
    u(0,x,y) =  0,
 stationary boundary conditions, i.e.
    u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
 and a heat source of the form
    h(x,y) = sin(pi*x)*sin(2*pi*y).

 Under this setup, the problem has an analytical solution:
    u(t,x,y) = a(t)*sin(pi*x)*sin(2*pi*y), where
    a(t) = (1 - exp(-(kx+4*ky)*pi^2*t)) / ((kx+4*ky)*pi^2).

 The spatial derivatives are computed using second-order
 centered differences, with the data distributed over nx*ny
 points on a uniform spatial grid.

 This program solves the problem with a DIRK method.  This
 employs a Newton iteration with a custom HYPRE-based,
 structured-grid, PCG+PFMG iterative linear solver.

 The custom SUNMatrix implementation is non-optimal, in that it
 must allocate an extraneous work array to use for adjusting
 matrix entries, since HYPRE's StructMatrix API contains *no*
 routines for zeroing out, scaling, adding an identity, or
 adding two matrices.

 20 outputs are printed at equal intervals, and run statistics
 are printed at the end.
---------------------------------------------------------------*/

// Header files
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>           // prototypes for ARKStep fcts., consts.
#include "nvector/nvector_parallel.h"        // parallel N_Vector types, fcts., macros
#include <sundials/sundials_linearsolver.h>  // definition of generic SUNLinearSolver object
#include <sundials/sundials_matrix.h>        // definition of generic SUNMatrix object
#include <sundials/sundials_nvector.h>       // definition of generic N_Vector object
#include "sundials/sundials_types.h"         // def. of type 'realtype'
#include "mpi.h"                             // MPI header file
#include "HYPRE_struct_ls.h"                 // HYPRE structured grid solver interface
#include "HYPRE_krylov.h"                    // other HYPRE utility routines

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

using namespace std;

// accessor macros between (x,y) location and 1D NVector array
#define IDX(x,y,n) ((n)*(y)+(x))
#define PI   RCONST(3.141592653589793238462643383279502884197169)
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define FOUR RCONST(4.0)

// user data structure
typedef struct {
  sunindextype nx;      // global number of x grid points
  sunindextype ny;      // global number of y grid points
  sunindextype is;      // global x indices of this subdomain
  sunindextype ie;
  sunindextype js;      // global y indices of this subdomain
  sunindextype je;
  sunindextype nxl;     // local number of x grid points
  sunindextype nyl;     // local number of y grid points
  realtype dx;          // x-directional mesh spacing
  realtype dy;          // y-directional mesh spacing
  realtype kx;          // x-directional diffusion coefficient
  realtype ky;          // y-directional diffusion coefficient
  N_Vector h;           // heat source vector
  MPI_Comm comm;        // communicator object
  int myid;             // MPI process ID
  int nprocs;           // total number of MPI processes
  bool HaveBdry[2][2];  // flags denoting if on physical boundary
  realtype *Erecv;      // receive buffers for neighbor exchange
  realtype *Wrecv;
  realtype *Nrecv;
  realtype *Srecv;
  realtype *Esend;      // send buffers for neighbor exchange
  realtype *Wsend;
  realtype *Nsend;
  realtype *Ssend;
} UserData;

//------------------------------------------------------------
// User-supplied functions called by the solver

//   ODE RHS function
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int J(realtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac,
             void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


//------------------------------------------------------------
// custom HYPRE 5-point Structured grid matrix definition

typedef struct _Hypre5ptMatrixContent {
  MPI_Comm           *comm;
  HYPRE_Int           ilower[2], iupper[2];
  HYPRE_StructGrid    grid;
  HYPRE_StructStencil stencil;
  HYPRE_StructMatrix  matrix;
  HYPRE_Real         *work;
  HYPRE_Int           nwork;
} *Hypre5ptMatrixContent;

#define H5PM_CONTENT(A)    ( (Hypre5ptMatrixContent)(A->content) )
#define H5PM_COMM(A)       ( *(H5PM_CONTENT(A)->comm) )
#define H5PM_ILOWER(A)     ( H5PM_CONTENT(A)->ilower )
#define H5PM_IUPPER(A)     ( H5PM_CONTENT(A)->iupper )
#define H5PM_GRID(A)       ( H5PM_CONTENT(A)->grid )
#define H5PM_STENCIL(A)    ( H5PM_CONTENT(A)->stencil )
#define H5PM_MATRIX(A)     ( H5PM_CONTENT(A)->matrix )
#define H5PM_WORK(A)       ( H5PM_CONTENT(A)->work )
#define H5PM_NWORK(A)      ( H5PM_CONTENT(A)->nwork )

//--- member function prototypes
SUNMatrix Hypre5ptMatrix(MPI_Comm &comm, sunindextype is, sunindextype ie,
                         sunindextype js, sunindextype je);
SUNMatrix_ID Hypre5ptMatrix_GetID(SUNMatrix A);
SUNMatrix Hypre5ptMatrix_Clone(SUNMatrix A);
void Hypre5ptMatrix_Destroy(SUNMatrix A);
int Hypre5ptMatrix_Zero(SUNMatrix A);
int Hypre5ptMatrix_Copy(SUNMatrix A, SUNMatrix B);
int Hypre5ptMatrix_ScaleAdd(realtype c, SUNMatrix A, SUNMatrix B);
int Hypre5ptMatrix_ScaleAddI(realtype c, SUNMatrix A);
int Hypre5ptMatrix_Matvec(SUNMatrix A, N_Vector x, N_Vector y);

//------------------------------------------------------------
// custom HYPRE PCG+PFMG solver definition

typedef struct _HyprePcgPfmgContent {
  HYPRE_StructVector bvec;
  HYPRE_StructVector xvec;
  HYPRE_StructSolver precond;
  HYPRE_StructSolver solver;
  realtype           resnorm;
  int                PCGits;
  int                PFMGits;
  long int           last_flag;
} *HyprePcgPfmgContent;

#define HPP_CONTENT(S)  ( (HyprePcgPfmgContent)(S->content) )
#define HPP_RESNORM(S)  ( HPP_CONTENT(S)->resnorm )
#define HPP_PCGITS(S)   ( HPP_CONTENT(S)->PCGits )
#define HPP_PFMGITS(S)  ( HPP_CONTENT(S)->PFMGits )
#define HPP_LASTFLAG(S) ( HPP_CONTENT(S)->last_flag )
#define HPP_X(S)        ( HPP_CONTENT(S)->xvec )
#define HPP_B(S)        ( HPP_CONTENT(S)->bvec )
#define HPP_PRECOND(S)  ( HPP_CONTENT(S)->precond )
#define HPP_SOLVER(S)   ( HPP_CONTENT(S)->solver )

//--- member function prototypes
SUNLinearSolver HyprePcgPfmg(SUNMatrix A, int PCGmaxit, int PFMGmaxit,
                             int relch, int rlxtype, int npre, int npost);
SUNLinearSolver_Type HyprePcgPfmg_GetType(SUNLinearSolver S);
int HyprePcgPfmg_Initialize(SUNLinearSolver S);
int HyprePcgPfmg_Setup(SUNLinearSolver S, SUNMatrix A);
int HyprePcgPfmg_Solve(SUNLinearSolver S, SUNMatrix A, 
                       N_Vector x, N_Vector b, realtype tol);
int HyprePcgPfmg_NumIters(SUNLinearSolver S);
realtype HyprePcgPfmg_ResNorm(SUNLinearSolver S);
long int HyprePcgPfmg_LastFlag(SUNLinearSolver S);
int HyprePcgPfmg_Free(SUNLinearSolver S);



//------------------------------------------------------------
// Private functions
//    checks function return values
static int check_flag(void *flagvalue, const string funcname, int opt);
//    sets default values into UserData structure
static int InitUserData(UserData *udata);
//    sets up parallel decomposition
static int SetupDecomp(UserData *udata);
//    performs neighbor exchange
static int Exchange(N_Vector y, UserData *udata);
//    frees memory allocated within UserData
static int FreeUserData(UserData *udata);
//    fills analytical solution
static int AnalyticalSolution(N_Vector ytrue, realtype t, UserData *udata);


//------------------------------------------------------------
// Main Program
int main(int argc, char* argv[]) {

  // general problem parameters
  realtype T0 = RCONST(0.0);   // initial time
  realtype Tf = RCONST(0.3);   // final time
  int Nt = 20;                 // total number of output times
  sunindextype nx = 60;        // spatial mesh size
  sunindextype ny = 120;
  realtype kx = 0.5;           // heat conductivity coefficients
  realtype ky = 0.75;
  realtype rtol = 1.e-5;       // relative and absolute tolerances
  realtype atol = 1.e-10;
  UserData *udata = NULL;
  realtype *data;
  sunindextype N, Ntot, i, j;

  // general problem variables
  int flag;                      // reusable error-checking flag
  int myid;                      // MPI process ID
  N_Vector y = NULL;             // empty vector for storing solution
  N_Vector ytrue = NULL;         // empty vector for storing reference solution
  N_Vector yerr = NULL;          // empty vector for storing solution error
  SUNMatrix A = NULL;            // empty matrix object for Newton method
  SUNLinearSolver LS = NULL;     // empty linear solver for Newton method
  void *arkode_mem = NULL;       // empty ARKode memory structure
  int PCGmaxit = 5;              // PCG+PFMG solver parameters
  int PFMGmaxit = 1;
  int relch = 1;
  int rlxtype = 3;
  int npre = 3;
  int npost = 3;

  // initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;
  flag = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

  // allocate and fill udata structure
  udata = new UserData;
  flag = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) return 1;
  udata->nx = nx;
  udata->ny = ny;
  udata->kx = kx;
  udata->ky = ky;
  udata->dx = ONE/(ONE*nx-ONE);   // x mesh spacing
  udata->dy = ONE/(ONE*ny-ONE);   // y mesh spacing

  // Set up parallel decomposition
  flag = SetupDecomp(udata);
  if (check_flag(&flag, "SetupDecomp", 1)) return 1;

  // Initial problem output
  bool outproc = (udata->myid == 0);
  if (outproc) {
    cout << "\n2D Heat PDE test problem:\n";
    cout << "   nprocs = " << udata->nprocs << "\n";
    cout << "   nx = " << udata->nx << "\n";
    cout << "   ny = " << udata->ny << "\n";
    cout << "   kx = " << udata->kx << "\n";
    cout << "   ky = " << udata->ky << "\n";
    cout << "   rtol = " << rtol << "\n";
    cout << "   atol = " << atol << "\n";
    cout << "   nxl (proc 0) = " << udata->nxl << "\n";
    cout << "   nyl (proc 0) = " << udata->nyl << "\n\n";
  }

  // open solver diagnostics output file for writing
  FILE *DFID;
  if (outproc)
    DFID=fopen("diags_ark_heat2D_hypre.txt","w");

  // Initialize data structures
  N = (udata->nxl)*(udata->nyl);
  Ntot = nx*ny;
  y = N_VNew_Parallel(udata->comm, N, Ntot);         // Create parallel vector for solution
  if (check_flag((void *) y, "N_VNew_Parallel", 0)) return 1;
  N_VConst(0.0, y);                                  // Set initial conditions
  ytrue = N_VNew_Parallel(udata->comm, N, Ntot);     // Create parallel vector for reference
  if (check_flag((void *) y, "N_VNew_Parallel", 0)) return 1;
  yerr = N_VNew_Parallel(udata->comm, N, Ntot);      // Create parallel vector for error
  if (check_flag((void *) y, "N_VNew_Parallel", 0)) return 1;
  udata->h = N_VNew_Parallel(udata->comm, N, Ntot);  // Create vector for heat source
  if (check_flag((void *) udata->h, "N_VNew_Parallel", 0)) return 1;

  // fill in the heat source array
  data = N_VGetArrayPointer(udata->h);
  for (j=0; j<udata->nyl; j++)
    for (i=0; i<udata->nxl; i++)
      data[IDX(i,j,udata->nxl)] = sin(PI*(udata->is+i)*udata->dx)
                                * sin(TWO*PI*(udata->js+j)*udata->dy);

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y);
  if (check_flag((void *) arkode_mem, "ARKStepCreate", 0)) return 1;

  // Set routines
  flag = ARKStepSetUserData(arkode_mem, (void *) udata);   // Pass udata to user functions
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
  flag = ARKStepSStolerances(arkode_mem, rtol, atol);
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  // create custom matrix and linear solver objects
  A = Hypre5ptMatrix(udata->comm, udata->is, udata->ie, udata->js, udata->je);
  if (check_flag((void *) A, "Hypre5ptMatrix", 0)) return 1;
  LS = HyprePcgPfmg(A, PCGmaxit, PFMGmaxit, relch, rlxtype, npre, npost);
  if (check_flag((void *) LS, "HyprePcgPfmg", 0)) return 1;

  // attach matrix, solver to ARKStep; set Jacobian construction routine
  flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
  flag = ARKStepSetJacFn(arkode_mem, J);
  if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;

  // Specify linearly implicit RHS, with non-time-dependent preconditioner
  flag = ARKStepSetLinear(arkode_mem, 0);
  if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;

  // Each processor outputs subdomain information
  char outname[100];
  sprintf(outname, "heat2d_subdomain.%03i.txt", udata->myid);
  FILE *UFID = fopen(outname,"w");
  fprintf(UFID, "%li  %li  %li  %li  %li  %li\n",
          (long int) udata->nx, (long int) udata->ny, (long int) udata->is,
          (long int) udata->ie, (long int) udata->js, (long int) udata->je);
  fclose(UFID);

  // Open output streams for results, access data array
  sprintf(outname, "heat2d.%03i.txt", udata->myid);
  UFID = fopen(outname,"w");
  data = N_VGetArrayPointer(y);

  // output initial condition to disk
  for (i=0; i<N; i++)  fprintf(UFID," %.16" ESYM, data[i]);
  fprintf(UFID,"\n");

  /* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  realtype t = T0;
  realtype dTout = (Tf-T0)/Nt;
  realtype tout = T0+dTout;
  realtype urms = sqrt(N_VDotProd(y,y)/nx/ny);
  realtype uerr, errI=0.0, err2=0.0;
  realtype hcur;
  if (outproc) {
    cout << "        t      ||u||_rms    ||uerr||      hcur\n\n";
    cout << "   ----------------------------------------------\n";
    printf("  %10.6" FSYM"  %10.6" FSYM"  %12.5" ESYM"    ----\n", t, urms, ZERO);
  }
  int iout;
  for (iout=0; iout<Nt; iout++) {

    flag = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKStepEvolve", 1)) {
      // unsuccessful step: break
      if (outproc)
        cerr << "Solver failure, stopping integration\n";
      break;
    } else {
      // successful step: update output time
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    urms = sqrt(N_VDotProd(y,y)/nx/ny);
    flag = AnalyticalSolution(ytrue, t, udata);
    if (check_flag(&flag, "AnalyticalSolution", 1)) break;
    N_VLinearSum( 1.0, ytrue, -1.0, y, yerr );
    uerr = N_VDotProd(yerr,yerr);
    errI = (errI > N_VMaxNorm(yerr)) ? errI : N_VMaxNorm(yerr);
    err2 += uerr;
    flag = ARKStepGetCurrentStep(arkode_mem, &hcur);
    if (check_flag(&flag, "ARKStepGetCurrentStep", 1)) break;
    if (outproc) printf("  %10.6" FSYM"  %10.6" FSYM"  %12.5" ESYM"  %8.1" ESYM"\n",
                        t, urms, sqrt(uerr/nx/ny), hcur);

    // output results to disk
    for (i=0; i<N; i++)  fprintf(UFID," %.16" ESYM"", data[i]);
    fprintf(UFID,"\n");
  }
  err2 = sqrt(err2 / nx / ny / Nt);
  if (outproc)  cout << "   ----------------------------------------------\n";
  fclose(UFID);

  // Print some final statistics
  long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nJe, nlcf, nni, ncfn, netf;
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return 1;
  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(&flag, "ARKStepGetNumStepAttempts", 1)) return 1;
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return 1;
  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1)) return 1;
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(&flag, "ARKStepGetNumErrTestFails", 1)) return 1;
  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return 1;
  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1)) return 1;
  flag = ARKStepGetNumLinIters(arkode_mem, &nli);
  if (check_flag(&flag, "ARKStepGetNumLinIters", 1)) return 1;
  flag = ARKStepGetNumJtimesEvals(arkode_mem, &nJv);
  if (check_flag(&flag, "ARKStepGetNumJtimesEvals", 1)) return 1;
  flag = ARKStepGetNumJacEvals(arkode_mem, &nJe);
  if (check_flag(&flag, "ARKStepGetNumJacEvals", 1)) return 1;
  flag = ARKStepGetNumLinConvFails(arkode_mem, &nlcf);
  if (check_flag(&flag, "ARKStepGetNumLinConvFails", 1)) return 1;

  if (outproc) {
    cout << "\nFinal Solver Statistics:\n";
    cout << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
    cout << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";
    cout << "   Total linear solver setups = " << nsetups << "\n";
    cout << "   Total linear iterations = " << nli << "\n";
    cout << "   Total number of Jacobian-vector products = " << nJv << "\n";
    cout << "   Total number of Jacobian evaluations = " << nJe << "\n";
    cout << "   Total number of linear solver convergence failures = " << nlcf << "\n";
    cout << "   Total number of Nonlinear iterations = " << nni << "\n";
    cout << "   Total number of nonlinear solver convergence failures = " << ncfn << "\n";
    cout << "   Total number of error test failures = " << netf << "\n";
    cout << "   Error: max = " << errI << ", rms = " << err2 << "\n";
  }

  // Clean up and return with successful completion
  ARKStepFree(&arkode_mem);    // Free integrator memory
  SUNLinSolFree(LS);           // Free linear solver memory
  SUNMatDestroy(A);            // Free matrix memory
  N_VDestroy(y);               // Free vectors
  N_VDestroy(ytrue);
  N_VDestroy(yerr);
  N_VDestroy(udata->h);
  FreeUserData(udata);         // Free user data
  delete udata;
  flag = MPI_Finalize();       // Finalize MPI
  return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

// f routine to compute the ODE RHS function f(t,y).
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(0.0, ydot);                           // Initialize ydot to zero
  UserData *udata = (UserData *) user_data;      // access problem data
  sunindextype nxl = udata->nxl;                 // set variable shortcuts
  sunindextype nyl = udata->nyl;
  realtype kx = udata->kx;
  realtype ky = udata->ky;
  realtype dx = udata->dx;
  realtype dy = udata->dy;
  realtype *Y = N_VGetArrayPointer(y);           // access data arrays
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;
  realtype *Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void *) Ydot, "N_VGetArrayPointer", 0)) return -1;

  // Exchange boundary data with neighbors
  int ierr = Exchange(y, udata);
  if (check_flag(&ierr, "Exchange", 1)) return -1;

  // iterate over subdomain interior, computing approximation to RHS
  realtype c1 = kx/dx/dx;
  realtype c2 = ky/dy/dy;
  realtype c3 = -TWO*(c1 + c2);
  sunindextype i, j;
  for (j=1; j<nyl-1; j++)                        // diffusive terms
    for (i=1; i<nxl-1; i++)
      Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + Y[IDX(i+1,j,nxl)])
                         + c2*(Y[IDX(i,j-1,nxl)] + Y[IDX(i,j+1,nxl)])
                         + c3*Y[IDX(i,j,nxl)];

  // iterate over subdomain boundaries (if not at overall domain boundary)
  if (!udata->HaveBdry[0][0]) {    // West face
    i=0;
    for (j=1; j<nyl-1; j++)
      Ydot[IDX(i,j,nxl)] = c1*(udata->Wrecv[j]   + Y[IDX(i+1,j,nxl)])
                         + c2*(Y[IDX(i,j-1,nxl)] + Y[IDX(i,j+1,nxl)])
                         + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][1]) {    // East face
    i=nxl-1;
    for (j=1; j<nyl-1; j++)
      Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + udata->Erecv[j])
                         + c2*(Y[IDX(i,j-1,nxl)] + Y[IDX(i,j+1,nxl)])
                         + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[1][0]) {    // South face
    j=0;
    for (i=1; i<nxl-1; i++)
      Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + Y[IDX(i+1,j,nxl)])
                         + c2*(udata->Srecv[i]   + Y[IDX(i,j+1,nxl)])
                         + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[1][1]) {    // West face
    j=nyl-1;
    for (i=1; i<nxl-1; i++)
      Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + Y[IDX(i+1,j,nxl)])
                         + c2*(Y[IDX(i,j-1,nxl)] + udata->Nrecv[i])
                         + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][0] && !udata->HaveBdry[1][0]) {  // South-West corner
    i = 0;
    j = 0;
    Ydot[IDX(i,j,nxl)] = c1*(udata->Wrecv[j] + Y[IDX(i+1,j,nxl)])
                       + c2*(udata->Srecv[i] + Y[IDX(i,j+1,nxl)])
                       + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][0] && !udata->HaveBdry[1][1]) {  // North-West corner
    i = 0;
    j = nyl-1;
    Ydot[IDX(i,j,nxl)] = c1*(udata->Wrecv[j]   + Y[IDX(i+1,j,nxl)])
                       + c2*(Y[IDX(i,j-1,nxl)] + udata->Nrecv[i])
                       + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][1] && !udata->HaveBdry[1][0]) {  // South-East corner
    i = nxl-1;
    j = 0;
    Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + udata->Erecv[j])
                       + c2*(udata->Srecv[i]   + Y[IDX(i,j+1,nxl)])
                       + c3*Y[IDX(i,j,nxl)];
  }
  if (!udata->HaveBdry[0][1] && !udata->HaveBdry[1][1]) {  // North-East corner
    i = nxl-1;
    j = nyl-1;
    Ydot[IDX(i,j,nxl)] = c1*(Y[IDX(i-1,j,nxl)] + udata->Erecv[j])
                       + c2*(Y[IDX(i,j-1,nxl)] + udata->Nrecv[i])
                       + c3*Y[IDX(i,j,nxl)];
  }

  // add in heat source
  N_VLinearSum(ONE, ydot, ONE, udata->h, ydot);
  return 0;     // Return with success
}


// J routine to compute the ODE RHS function Jacobian, (df/dy)(t,y).
static int J(realtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac,
             void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // access problem data, and set variable shortcuts
  UserData    *udata = (UserData *) user_data;
  sunindextype nxl = udata->nxl;
  sunindextype nyl = udata->nyl;
  realtype     kx = udata->kx;
  realtype     ky = udata->ky;
  realtype     dx = udata->dx;
  realtype     dy = udata->dy;
  realtype     c1 = kx/dx/dx;
  realtype     c2 = ky/dy/dy;
  realtype     c3 = -TWO*(c1 + c2);
  HYPRE_Real   vals[5] = {c2, c1, c3, c1, c2};
  HYPRE_Int    entries[5] = {0, 1, 2, 3, 4};
  HYPRE_Int    ix, iy, is, ie, js, je, index[2];
  int          ierr;

  // iterate over subdomain interior setting stencil entries
  is = (udata->HaveBdry[0][0]) ? 1 : 0;
  ie = (udata->HaveBdry[0][1]) ? nxl-1 : nxl;
  js = (udata->HaveBdry[1][0]) ? 1 : 0;
  je = (udata->HaveBdry[1][1]) ? nyl-1 : nyl;
  for (iy=js; iy<je; iy++) {
    index[1] = H5PM_ILOWER(Jac)[1] + iy;
    for (ix=is; ix<ie; ix++) {
      index[0] = H5PM_ILOWER(Jac)[0] + ix;
      ierr = HYPRE_StructMatrixSetValues(H5PM_MATRIX(Jac), index, 5, entries, vals);
      if (ierr != 0)  return(ierr);
    }
  }

  // assemble matrix
  ierr = HYPRE_StructMatrixAssemble(H5PM_MATRIX(Jac));

  // Return with success
  return 0;
}


/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void *flagvalue, const string funcname, int opt)
{
  int *errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL) {
    cerr << "\nSUNDIALS_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1; }

  // Check if flag < 0
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      cerr << "\nSUNDIALS_ERROR: " << funcname << " failed with flag = " << *errflag << "\n\n";
      return 1;
    }
  }

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && flagvalue == NULL) {
    cerr << "\nMEMORY_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1; }

  return 0;
}

// Set up parallel decomposition
static int SetupDecomp(UserData *udata)
{
  // check that this has not been called before
  if (udata->Erecv != NULL || udata->Wrecv != NULL ||
      udata->Srecv != NULL || udata->Nrecv != NULL) {
    cerr << "SetupDecomp warning: parallel decomposition already set up\n";
    return 1;
  }

  // get suggested parallel decomposition
  int ierr, dims[] = {0, 0};
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &(udata->nprocs));
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Comm_size = " << ierr << "\n";
    return -1;
  }
  ierr = MPI_Dims_create(udata->nprocs, 2, dims);
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Dims_create = " << ierr << "\n";
    return -1;
  }

  // set up 2D Cartesian communicator
  int periods[] = {0, 0};
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &(udata->comm));
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Cart_create = " << ierr << "\n";
    return -1;
  }
  ierr = MPI_Comm_rank(udata->comm, &(udata->myid));
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Comm_rank = " << ierr << "\n";
    return -1;
  }

  // determine local extents
  int coords[2];
  ierr = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Cart_get = " << ierr << "\n";
    return -1;
  }
  udata->is = (udata->nx)*(coords[0])/(dims[0]);
  udata->ie = (udata->nx)*(coords[0]+1)/(dims[0])-1;
  udata->js = (udata->ny)*(coords[1])/(dims[1]);
  udata->je = (udata->ny)*(coords[1]+1)/(dims[1])-1;
  udata->nxl = (udata->ie)-(udata->is)+1;
  udata->nyl = (udata->je)-(udata->js)+1;

  // determine if I have neighbors, and allocate exchange buffers
  udata->HaveBdry[0][0] = (udata->is == 0);
  udata->HaveBdry[0][1] = (udata->ie == udata->nx-1);
  udata->HaveBdry[1][0] = (udata->js == 0);
  udata->HaveBdry[1][1] = (udata->je == udata->ny-1);
  if (!udata->HaveBdry[0][0]) {
    udata->Wrecv = new realtype[udata->nyl];
    udata->Wsend = new realtype[udata->nyl];
  }
  if (!udata->HaveBdry[0][1]) {
    udata->Erecv = new realtype[udata->nyl];
    udata->Esend = new realtype[udata->nyl];
  }
  if (!udata->HaveBdry[1][0]) {
    udata->Srecv = new realtype[udata->nxl];
    udata->Ssend = new realtype[udata->nxl];
  }
  if (!udata->HaveBdry[1][1]) {
    udata->Nrecv = new realtype[udata->nxl];
    udata->Nsend = new realtype[udata->nxl];
  }

  return 0;     // return with success flag
}

// Perform neighbor exchange
static int Exchange(N_Vector y, UserData *udata)
{
  // local variables
  MPI_Request reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN;
  MPI_Status stat;
  int ierr, i, ipW=-1, ipE=-1, ipS=-1, ipN=-1;
  int coords[2], dims[2], periods[2], nbcoords[2];
  int nyl = udata->nyl;
  int nxl = udata->nxl;

  // access data array
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;

  // MPI equivalent of realtype type
  #if defined(SUNDIALS_SINGLE_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_FLOAT
  #elif defined(SUNDIALS_DOUBLE_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_DOUBLE
  #elif defined(SUNDIALS_EXTENDED_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_LONG_DOUBLE
  #endif

  // MPI neighborhood information
  ierr = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
  if (ierr != MPI_SUCCESS) {
    cerr << "Error in MPI_Cart_get = " << ierr << "\n";
    return -1;
  }
  if (!udata->HaveBdry[0][0]) {
    nbcoords[0] = coords[0]-1;
    nbcoords[1] = coords[1];
    ierr = MPI_Cart_rank(udata->comm, nbcoords, &ipW);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    nbcoords[0] = coords[0]+1;
    nbcoords[1] = coords[1];
    ierr = MPI_Cart_rank(udata->comm, nbcoords, &ipE);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]-1;
    ierr = MPI_Cart_rank(udata->comm, nbcoords, &ipS);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]+1;
    ierr = MPI_Cart_rank(udata->comm, nbcoords, &ipN);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }

  // open Irecv buffers
  if (!udata->HaveBdry[0][0]) {
    ierr = MPI_Irecv(udata->Wrecv, udata->nyl, REALTYPE_MPI_TYPE, ipW,
                   MPI_ANY_TAG, udata->comm, &reqRW);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    ierr = MPI_Irecv(udata->Erecv, udata->nyl, REALTYPE_MPI_TYPE, ipE,
                   MPI_ANY_TAG, udata->comm, &reqRE);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    ierr = MPI_Irecv(udata->Srecv, udata->nxl, REALTYPE_MPI_TYPE, ipS,
                   MPI_ANY_TAG, udata->comm, &reqRS);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    ierr = MPI_Irecv(udata->Nrecv, udata->nxl, REALTYPE_MPI_TYPE, ipN,
                   MPI_ANY_TAG, udata->comm, &reqRN);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }

  // send data
  if (!udata->HaveBdry[0][0]) {
    for (i=0; i<nyl; i++)  udata->Wsend[i] = Y[IDX(0,i,nxl)];
    ierr = MPI_Isend(udata->Wsend, udata->nyl, REALTYPE_MPI_TYPE, ipW, 0,
              udata->comm, &reqSW);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    for (i=0; i<nyl; i++)  udata->Esend[i] = Y[IDX(nxl-1,i,nxl)];
    ierr = MPI_Isend(udata->Esend, udata->nyl, REALTYPE_MPI_TYPE, ipE, 1,
              udata->comm, &reqSE);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    for (i=0; i<nxl; i++)  udata->Ssend[i] = Y[IDX(i,0,nxl)];
    ierr = MPI_Isend(udata->Ssend, udata->nxl, REALTYPE_MPI_TYPE, ipS, 2,
              udata->comm, &reqSS);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    for (i=0; i<nxl; i++)  udata->Nsend[i] = Y[IDX(i,nyl-1,nxl)];
    ierr = MPI_Isend(udata->Nsend, udata->nxl, REALTYPE_MPI_TYPE, ipN, 3,
              udata->comm, &reqSN);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }

  // wait for messages to finish
  if (!udata->HaveBdry[0][0]) {
    ierr = MPI_Wait(&reqRW, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSW, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    ierr = MPI_Wait(&reqRE, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSE, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    ierr = MPI_Wait(&reqRS, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSS, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    ierr = MPI_Wait(&reqRN, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSN, &stat);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }

  return 0;     // return with success flag
}

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
  udata->nx = 0;
  udata->ny = 0;
  udata->is = 0;
  udata->ie = 0;
  udata->js = 0;
  udata->je = 0;
  udata->nxl = 0;
  udata->nyl = 0;
  udata->dx = ZERO;
  udata->dy = ZERO;
  udata->kx = ZERO;
  udata->ky = ZERO;
  udata->h = NULL;
  udata->comm = MPI_COMM_WORLD;
  udata->myid = 0;
  udata->nprocs = 0;
  udata->HaveBdry[0][0] = 1;
  udata->HaveBdry[0][1] = 1;
  udata->HaveBdry[1][0] = 1;
  udata->HaveBdry[1][1] = 1;
  udata->Erecv = NULL;
  udata->Wrecv = NULL;
  udata->Nrecv = NULL;
  udata->Srecv = NULL;
  udata->Esend = NULL;
  udata->Wsend = NULL;
  udata->Nsend = NULL;
  udata->Ssend = NULL;

  return 0;     // return with success flag
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{
  // free exchange buffers
  if (udata->Wrecv != NULL)  delete[] udata->Wrecv;
  if (udata->Wsend != NULL)  delete[] udata->Wsend;
  if (udata->Erecv != NULL)  delete[] udata->Erecv;
  if (udata->Esend != NULL)  delete[] udata->Esend;
  if (udata->Srecv != NULL)  delete[] udata->Srecv;
  if (udata->Ssend != NULL)  delete[] udata->Ssend;
  if (udata->Nrecv != NULL)  delete[] udata->Nrecv;
  if (udata->Nsend != NULL)  delete[] udata->Nsend;

  return 0;     // return with success flag
}

// Fills analytical solution vector:
//     u(t,x,y) = a(t)*sin(pi*x)*sin(2*pi*y), where
//     a(t) = (1 - exp(-(kx+4*ky)*pi^2*t)) / ((kx+4*ky)*pi^2).
static int AnalyticalSolution(N_Vector ytrue, realtype t, UserData *udata)
{
  long int i, j;
  realtype *y, at;

  // set time-dependent portion
  at = (ONE - exp(-(udata->kx + FOUR*udata->ky)*PI*PI*t))
     / ((udata->kx + FOUR*udata->ky)*PI*PI);

  // iterate over the domain, filling solution
  y = N_VGetArrayPointer(ytrue);
  for (j=0; j<udata->nyl; j++)
    for (i=0; i<udata->nxl; i++)
      y[IDX(i,j,udata->nxl)] = at * sin(PI*(udata->is+i)*udata->dx)
                                  * sin(TWO*PI*(udata->js+j)*udata->dy);

}


//------------------------------------------------------------
// custom HYPRE 5-point Structured grid matrix implementation

SUNMatrix Hypre5ptMatrix(MPI_Comm &comm, sunindextype is, sunindextype ie,
                         sunindextype js, sunindextype je) {
  SUNMatrix A;
  SUNMatrix_Ops ops;
  Hypre5ptMatrixContent content;
  HYPRE_Int offset[2];
  int ierr, result;

  // return with NULL matrix on illegal input
  if ((je<js) || (ie<is))  return(NULL);

  // check for valid 2D Cartesian MPI communicator
  ierr = MPI_Topo_test(comm, &result);
  if ((ierr != MPI_SUCCESS) || (result != MPI_CART))  return(NULL);
  ierr = MPI_Cartdim_get(comm, &result);
  if ((ierr != MPI_SUCCESS) || (result != 2))  return(NULL);

  // Create matrix
  A = NULL;
  A = (SUNMatrix) malloc(sizeof *A);
  if (A == NULL) return(NULL);
  memset(A, 0, sizeof(struct _generic_SUNMatrix));

  // Create matrix operation structure
  ops = NULL;
  ops = (SUNMatrix_Ops) malloc(sizeof(struct _generic_SUNMatrix_Ops));
  if (ops == NULL) { free(A); return(NULL); }
  memset(ops, 0, sizeof(struct _generic_SUNMatrix_Ops));

  // Attach operations
  ops->getid     = Hypre5ptMatrix_GetID;
  ops->clone     = Hypre5ptMatrix_Clone;
  ops->destroy   = Hypre5ptMatrix_Destroy;
  ops->zero      = Hypre5ptMatrix_Zero;
  ops->copy      = Hypre5ptMatrix_Copy;
  ops->scaleadd  = Hypre5ptMatrix_ScaleAdd;
  ops->scaleaddi = Hypre5ptMatrix_ScaleAddI;
  ops->matvec    = Hypre5ptMatrix_Matvec;
  ops->space     = NULL;

  // Create content
  content = NULL;
  content = (Hypre5ptMatrixContent) malloc(sizeof(struct _Hypre5ptMatrixContent));
  if (content == NULL) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  memset(content, 0, sizeof(struct _Hypre5ptMatrixContent));

  // Fill content

  //    ilower/iupper
  content->ilower[0] = is;  content->iupper[0] = ie;
  content->ilower[1] = js;  content->iupper[1] = je;

  //    comm
  content->comm = &comm;

  //    grid
  ierr = HYPRE_StructGridCreate(comm, 2, &(content->grid));
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  ierr = HYPRE_StructGridSetExtents(content->grid, content->ilower, content->iupper);
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  ierr = HYPRE_StructGridAssemble(content->grid);
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }

  //    stencil
  ierr = HYPRE_StructStencilCreate(2, 5, &(content->stencil));
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  offset[0] = 0;  offset[1] = -1;    // dependency to bottom
  ierr = HYPRE_StructStencilSetElement(content->stencil, 0, offset);
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  offset[0] = -1;  offset[1] = 0;    // dependency to left
  ierr = HYPRE_StructStencilSetElement(content->stencil, 1, offset);
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  offset[0] = 0;  offset[1] = 0;     // dependency to self
  ierr = HYPRE_StructStencilSetElement(content->stencil, 2, offset);
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  offset[0] = 1;  offset[1] = 0;     // dependency to right
  ierr = HYPRE_StructStencilSetElement(content->stencil, 3, offset);
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  offset[0] = 0;  offset[1] = 1;     // dependency to top
  ierr = HYPRE_StructStencilSetElement(content->stencil, 4, offset);
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }

  //    matrix
  ierr = HYPRE_StructMatrixCreate(comm, content->grid, content->stencil,
                                  &(content->matrix));
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }
  ierr = HYPRE_StructMatrixInitialize(content->matrix);
  if (ierr != 0) { Hypre5ptMatrix_Destroy(A); return(NULL); }

  //    work array
  content->nwork = 5 * (ie-is+1) * (je-js+1);
  content->work = NULL;
  content->work = new HYPRE_Real[content->nwork];
  if (content->work == NULL) { Hypre5ptMatrix_Destroy(A); return(NULL); }

  // Attach content and ops
  A->content = content;
  A->ops     = ops;

  return(A);
}

SUNMatrix_ID Hypre5ptMatrix_GetID(SUNMatrix A) {
  return SUNMATRIX_CUSTOM;
}

SUNMatrix Hypre5ptMatrix_Clone(SUNMatrix A) {
  SUNMatrix B = Hypre5ptMatrix(H5PM_COMM(A), H5PM_ILOWER(A)[0], H5PM_IUPPER(A)[0], 
                               H5PM_ILOWER(A)[1], H5PM_IUPPER(A)[1]);
  return(B);
}

void Hypre5ptMatrix_Destroy(SUNMatrix A) {
  if (A == NULL)  return;
  if (A->ops)  free(A->ops);
  if (A->content == NULL) {
    free(A);
    return;
  }
  if (H5PM_WORK(A))    delete[] H5PM_WORK(A);
  if (H5PM_MATRIX(A))  HYPRE_StructMatrixDestroy(H5PM_MATRIX(A));
  if (H5PM_GRID(A))    HYPRE_StructGridDestroy(H5PM_GRID(A));
  if (H5PM_STENCIL(A)) HYPRE_StructStencilDestroy(H5PM_STENCIL(A));
  free(A->content);
  free(A);
}

int Hypre5ptMatrix_Zero(SUNMatrix A) {
  int ierr, i;
  HYPRE_Int entries[5] = {0,1,2,3,4};

  // set work array to all zeros
  for (i=0; i<H5PM_NWORK(A); i++)
    H5PM_WORK(A)[i] = ZERO;

  // set values into matrix
  ierr = HYPRE_StructMatrixSetBoxValues(H5PM_MATRIX(A), H5PM_ILOWER(A),
                                        H5PM_IUPPER(A), 5, entries, H5PM_WORK(A));
  return(ierr);
}

int Hypre5ptMatrix_Copy(SUNMatrix A, SUNMatrix B) {
  int ierr, i;
  HYPRE_Int entries[5] = {0,1,2,3,4};

  // copy values from A into work array
  ierr = HYPRE_StructMatrixGetBoxValues(H5PM_MATRIX(A), H5PM_ILOWER(A),
                                        H5PM_IUPPER(A), 5, entries, H5PM_WORK(A));
  if (ierr != 0)  return(ierr);

  // insert values into B
  ierr = HYPRE_StructMatrixSetBoxValues(H5PM_MATRIX(B), H5PM_ILOWER(A),
                                        H5PM_IUPPER(A), 5, entries, H5PM_WORK(A));
  return(ierr);
}

int Hypre5ptMatrix_ScaleAdd(realtype c, SUNMatrix A, SUNMatrix B) {
  int ierr, i;
  HYPRE_Int entries[5] = {0,1,2,3,4};

  // copy values from A into work array
  ierr = HYPRE_StructMatrixGetBoxValues(H5PM_MATRIX(A), H5PM_ILOWER(A),
                                        H5PM_IUPPER(A), 5, entries, H5PM_WORK(A));
  if (ierr != 0)  return(ierr);

  // scale work array by c
  for (i=0; i<H5PM_NWORK(A); i++)
    H5PM_WORK(A)[i] *= c;

  // add resulting values to B
  ierr = HYPRE_StructMatrixAddToBoxValues(H5PM_MATRIX(B), H5PM_ILOWER(A),
                                          H5PM_IUPPER(A), 5, entries, H5PM_WORK(A));
  return(ierr);
}

int Hypre5ptMatrix_ScaleAddI(realtype c, SUNMatrix A) {
  int ierr, i;
  HYPRE_Int entries[5] = {0,1,2,3,4};

  // set work array to all ones
  ierr = HYPRE_StructMatrixGetBoxValues(H5PM_MATRIX(A), H5PM_ILOWER(A),
                                        H5PM_IUPPER(A), 5, entries, H5PM_WORK(A));
  if (ierr != 0)  return(ierr);

  // scale work array by c
  for (i=0; i<H5PM_NWORK(A); i++)
    H5PM_WORK(A)[i] *= c;

  // insert resulting values back into A
  ierr = HYPRE_StructMatrixSetBoxValues(H5PM_MATRIX(A),
                                        H5PM_ILOWER(A),
                                        H5PM_IUPPER(A), 5,
                                        entries, H5PM_WORK(A));

  // set first 1/5 of work array to 1
  for (i=0; i<H5PM_NWORK(A)/5; i++)
    H5PM_WORK(A)[i] = ONE;

  // insert resulting values back into diagonal of A
  entries[0] = 2;
  ierr = HYPRE_StructMatrixAddToBoxValues(H5PM_MATRIX(A),
                                          H5PM_ILOWER(A),
                                          H5PM_IUPPER(A), 1,
                                          entries, H5PM_WORK(A));
  return(ierr);
}

int Hypre5ptMatrix_Matvec(SUNMatrix A, N_Vector x, N_Vector y) {
  HYPRE_StructVector x_hypre, y_hypre;
  int ierr;

  // create temporary HYPRE_StructVectors to hold x and y
  ierr = HYPRE_StructVectorCreate(H5PM_COMM(A), H5PM_GRID(A), &x_hypre);
  if (ierr != 0)  return(ierr);
  ierr = HYPRE_StructVectorInitialize(x_hypre);
  if (ierr != 0)  return(ierr);
  ierr = HYPRE_StructVectorCreate(H5PM_COMM(A), H5PM_GRID(A), &y_hypre);
  if (ierr != 0)  return(ierr);
  ierr = HYPRE_StructVectorInitialize(y_hypre);
  if (ierr != 0)  return(ierr);

  // insert x into x_hypre
  ierr = HYPRE_StructVectorSetBoxValues(x_hypre, H5PM_ILOWER(A),
                                        H5PM_IUPPER(A),
                                        N_VGetArrayPointer(x));
  if (ierr != 0)  return(ierr);

  // Perform operation
  ierr = HYPRE_StructMatrixMatvec(ONE, H5PM_MATRIX(A), x_hypre, ZERO, y_hypre);
  if (ierr != 0)  return(ierr);

  // insert y_hypre into y and return
  ierr = HYPRE_StructVectorGetBoxValues(y_hypre, H5PM_ILOWER(A),
                                        H5PM_IUPPER(A),
                                        N_VGetArrayPointer(y));
  return(ierr);
}

//------------------------------------------------------------
// custom HYPRE PCG+PFMG solver implementation

SUNLinearSolver HyprePcgPfmg(SUNMatrix A, int PCGmaxit, int PFMGmaxit,
                             int relch, int rlxtype, int npre, int npost) {
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  HyprePcgPfmgContent content;
  int ierr, result;

  // return with NULL solver on illegal inputs
  if ((PCGmaxit < 1) || (PFMGmaxit < 1) || (npre < 1) || (npost < 1))
    return(NULL);

  // check for valid 2D Cartesian MPI communicator
  ierr = MPI_Topo_test(H5PM_COMM(A), &result);
  if ((ierr != MPI_SUCCESS) || (result != MPI_CART))  return(NULL);
  ierr = MPI_Cartdim_get(H5PM_COMM(A), &result);
  if ((ierr != MPI_SUCCESS) || (result != 2))  return(NULL);

  // Create linear solver
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  memset(S, 0, sizeof(struct _generic_SUNLinearSolver));

  // Create linear solver operation structure
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { delete S; return(NULL); }
  memset(ops, 0, sizeof(struct _generic_SUNLinearSolver_Ops));

  // Attach operations
  ops->gettype           = HyprePcgPfmg_GetType;
  ops->initialize        = HyprePcgPfmg_Initialize;
  ops->setatimes         = NULL;
  ops->setpreconditioner = NULL;
  ops->setscalingvectors = NULL;
  ops->setup             = HyprePcgPfmg_Setup;
  ops->solve             = HyprePcgPfmg_Solve;
  ops->numiters          = HyprePcgPfmg_NumIters;
  ops->resnorm           = HyprePcgPfmg_ResNorm;
  ops->resid             = NULL;
  ops->lastflag          = HyprePcgPfmg_LastFlag;
  ops->space             = NULL;
  ops->free              = HyprePcgPfmg_Free;

  // Create content
  content = NULL;
  content = (HyprePcgPfmgContent) malloc(sizeof(struct _HyprePcgPfmgContent));
  if (content == NULL) { HyprePcgPfmg_Free(S); return(NULL); }
  memset(content, 0, sizeof(struct _HyprePcgPfmgContent));

  //     xvec
  ierr = HYPRE_StructVectorCreate(H5PM_COMM(A), H5PM_GRID(A), &(content->xvec));
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructVectorInitialize(content->xvec);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }

  //     bvec
  ierr = HYPRE_StructVectorCreate(H5PM_COMM(A), H5PM_GRID(A), &(content->bvec));
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructVectorInitialize(content->bvec);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }

  //     precond
  ierr = HYPRE_StructPFMGCreate(H5PM_COMM(A), &(content->precond));
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructPFMGSetMaxIter(content->precond, PFMGmaxit);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructPFMGSetRelaxType(content->precond, rlxtype);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructPFMGSetNumPreRelax(content->precond, npre);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructPFMGSetNumPostRelax(content->precond, npost);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }

  //     solver
  ierr = HYPRE_StructPCGCreate(H5PM_COMM(A), &(content->solver));
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructPCGSetMaxIter(content->solver, PCGmaxit);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructPCGSetRelChange(content->solver, relch);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructPCGSetLogging(content->solver, 1);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }
  ierr = HYPRE_StructPCGSetPrecond(content->solver,
                                   (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,
                                   (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup,
                                   content->precond);
  if (ierr != 0)  { HyprePcgPfmg_Free(S); return(NULL); }

  //     statistics
  content->resnorm   = ZERO;
  content->PCGits    = 0;
  content->PFMGits   = 0;
  content->last_flag = SUNLS_SUCCESS;

  // Attach content and ops
  S->content = content;
  S->ops     = ops;

  return(S);
}

SUNLinearSolver_Type HyprePcgPfmg_GetType(SUNLinearSolver S) {
  return(SUNLINEARSOLVER_MATRIX_ITERATIVE);
}

int HyprePcgPfmg_Initialize(SUNLinearSolver S) {
  // no additional memory to allocate; return with success
  HPP_LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUNLS_SUCCESS);
}

int HyprePcgPfmg_Setup(SUNLinearSolver S, SUNMatrix A) {
  int ierr;

  // set rhs/solution vectors as all zero for now
  ierr = HYPRE_StructVectorSetConstantValues(HPP_B(S), ZERO);
  if (ierr != 0)  return(ierr);
  ierr = HYPRE_StructVectorAssemble(HPP_B(S));
  if (ierr != 0)  return(ierr);
  ierr = HYPRE_StructVectorSetConstantValues(HPP_X(S), ZERO);
  if (ierr != 0)  return(ierr);
  ierr = HYPRE_StructVectorAssemble(HPP_X(S));
  if (ierr != 0)  return(ierr);

  // set up the solver
  ierr = HYPRE_StructPCGSetup(HPP_SOLVER(S), H5PM_MATRIX(A),
                              HPP_B(S), HPP_X(S));
  if (ierr != 0)  return(ierr);

  // return with success
  HPP_LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUNLS_SUCCESS);
}

int HyprePcgPfmg_Solve(SUNLinearSolver S, SUNMatrix A,
                       N_Vector x, N_Vector b, realtype tol) {
  HYPRE_Real finalresid;
  HYPRE_Int PCGits, PFMGits, converged;
  int ierr;

  // supply the desired [absolute] linear solve tolerance to HYPRE
  ierr = HYPRE_StructPCGSetAbsoluteTol(HPP_SOLVER(S), tol);
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }
  ierr = HYPRE_StructPCGSetTol(HPP_SOLVER(S), ZERO);
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }


  // insert rhs N_Vector entries into HYPRE vector b and assemble
  ierr = HYPRE_StructVectorSetBoxValues(HPP_B(S), H5PM_ILOWER(A),
                                        H5PM_IUPPER(A),
                                        N_VGetArrayPointer(b));
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }
  ierr = HYPRE_StructVectorAssemble(HPP_X(S));
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }

  // insert initial guess N_Vector entries into HYPRE vector x and assemble
  ierr = HYPRE_StructVectorSetBoxValues(HPP_X(S), H5PM_ILOWER(A),
                                        H5PM_IUPPER(A),
                                        N_VGetArrayPointer(x));
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }
  ierr = HYPRE_StructVectorAssemble(HPP_B(S));
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }

  // solve the linear system
  ierr = HYPRE_StructPCGSolve(HPP_SOLVER(S), H5PM_MATRIX(A),
                              HPP_B(S), HPP_X(S));
  if (ierr == 0) {
    HPP_LASTFLAG(S) = SUNLS_SUCCESS;
  } else if (ierr == 256) {
    HPP_LASTFLAG(S) = SUNLS_RES_REDUCED;
  } else {
    HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S));
  }

  // extract solver statistics, and store for later
  ierr = HYPRE_StructPCGGetFinalRelativeResidualNorm(HPP_SOLVER(S), &finalresid);
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }
  ierr = HYPRE_StructPCGGetNumIterations(HPP_SOLVER(S), &PCGits);
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }
  ierr = HYPRE_StructPFMGGetNumIterations(HPP_PRECOND(S), &PFMGits);
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }
  HPP_RESNORM(S) = finalresid;
  HPP_PCGITS(S)  = PCGits;
  HPP_PFMGITS(S) = PFMGits;

  // extract solution values
  ierr = HYPRE_StructVectorGetBoxValues(HPP_X(S), H5PM_ILOWER(A),
                                        H5PM_IUPPER(A),
                                        N_VGetArrayPointer(x));
  if (ierr != 0)  { HPP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;  return(HPP_LASTFLAG(S)); }

  // solve finished, return with solver result (stored in HPP_LASTFLAG(S))
  return(HPP_LASTFLAG(S));
}


int HyprePcgPfmg_NumIters(SUNLinearSolver S) {
  // return the stored number of outer PCG iterations
  if (S == NULL) return(-1);
  return (HPP_PCGITS(S));
}


realtype HyprePcgPfmg_ResNorm(SUNLinearSolver S) {
  // return the stored 'resnorm' value
  if (S == NULL) return(-ONE);
  return (HPP_RESNORM(S));
}


long int HyprePcgPfmg_LastFlag(SUNLinearSolver S) {
  // return the stored 'last_flag' value
  if (S == NULL) return(-1);
  return (HPP_LASTFLAG(S));
}


int HyprePcgPfmg_Free(SUNLinearSolver S) {

  // check for easy return
  if (S == NULL) return(SUNLS_SUCCESS);
  if (S->ops)  free(S->ops);
  if (S->content == NULL) {
    free(S);
    return(SUNLS_SUCCESS);
  }

  // delete items from within the content structure
  if (HPP_SOLVER(S))   HYPRE_StructPCGDestroy(HPP_SOLVER(S));
  if (HPP_PRECOND(S))  HYPRE_StructPFMGDestroy(HPP_PRECOND(S));
  if (HPP_B(S))        HYPRE_StructVectorDestroy(HPP_B(S));
  if (HPP_X(S))        HYPRE_StructVectorDestroy(HPP_X(S));

  // delete content structure itself, as well as overall object, and return
  free(S->content);
  free(S);
  return(SUNLS_SUCCESS);
}


//---- end of file ----
