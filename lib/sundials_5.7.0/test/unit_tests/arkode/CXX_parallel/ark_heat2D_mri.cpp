/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Routine to check that MRIStep and ARKStep exhibit the same
 * solver statistics when both run with fixed-steps, the same
 * solver parameters, and MRIStep runs using a solve-decoupled
 * DIRK method at the slow time scale.
 *
 * This routine will switch between the default Newton nonlinear
 * solver and the 'linear' version based on a 0/1 command-line
 * argument (1 => linear).
 *---------------------------------------------------------------*/

// Header files
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include "arkode/arkode_arkstep.h"    // prototypes for ARKStep fcts., consts
#include "arkode/arkode_mristep.h"    // prototypes for MRIStep fcts., consts
#include "nvector/nvector_parallel.h" // parallel N_Vector types, fcts., macros
#include "sunlinsol/sunlinsol_pcg.h"  // access to PCG SUNLinearSolver
#include "sundials/sundials_types.h"  // def. of type 'realtype'
#include "mpi.h"                      // MPI header file

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

// user data structure
typedef struct {
  sunindextype nx;          // global number of x grid points
  sunindextype ny;          // global number of y grid points
  sunindextype is;          // global x indices of this subdomain
  sunindextype ie;
  sunindextype js;          // global y indices of this subdomain
  sunindextype je;
  sunindextype nxl;         // local number of x grid points
  sunindextype nyl;         // local number of y grid points
  realtype dx;          // x-directional mesh spacing
  realtype dy;          // y-directional mesh spacing
  realtype kx;          // x-directional diffusion coefficient
  realtype ky;          // y-directional diffusion coefficient
  N_Vector h;           // heat source vector
  N_Vector d;           // inverse of Jacobian diagonal
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

// User-supplied Functions Called by the Solver
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int PSet(realtype t, N_Vector y, N_Vector fy, booleantype jok,
                booleantype *jcurPtr, realtype gamma, void *user_data);
static int PSol(realtype t, N_Vector y, N_Vector fy, N_Vector r,
                N_Vector z, realtype gamma, realtype delta, int lr,
                void *user_data);

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

// Main Program
int main(int argc, char* argv[]) {

  // general problem parameters
  realtype T0 = RCONST(0.0);     // initial time
  realtype Tf = RCONST(0.3);     // final time
  int Nt = 1000;                 // total number of internal steps
  sunindextype nx = 60;          // spatial mesh size
  sunindextype ny = 120;
  realtype kx = RCONST(0.5);     // heat conductivity coefficients
  realtype ky = RCONST(0.75);
  realtype rtol = RCONST(1.e-5); // relative and absolute tolerances
  realtype atol = RCONST(1.e-10);
  UserData *udata = NULL;
  realtype *data;
  sunindextype N, Ntot, i, j;
  int numfails;
  booleantype linear;
  realtype t;
  long int ark_nst, ark_nfe, ark_nfi, ark_nsetups, ark_nli, ark_nJv, ark_nlcf, ark_nni, ark_ncfn, ark_npe, ark_nps;
  long int mri_nst, mri_nfs, mri_nsetups, mri_nli, mri_nJv, mri_nlcf, mri_nni, mri_ncfn, mri_npe, mri_nps;

  // general problem variables
  int flag;                       // reusable error-checking flag
  int myid;                       // MPI process ID
  N_Vector y = NULL;              // empty vector for storing solution
  SUNLinearSolver LSa = NULL;     // empty linear solver memory structures
  SUNLinearSolver LSm = NULL;
  void *arkstep_mem = NULL;       // empty ARKStep memory structure
  void *mristep_mem = NULL;       // empty MRIStep memory structure
  void *inner_mem = NULL;         // empty inner ARKStep memory structure
  ARKodeButcherTable B = NULL;    // empty Butcher table

  // initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;
  flag = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

  // if an argument supplied, set linear (otherwise use SUNFALSE)
  linear = SUNFALSE;
  if (argc > 1)  linear = stoi(argv[1], NULL);

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
    cout << "   nyl (proc 0) = " << udata->nyl << "\n";
    if (linear) {
      cout << "   Linearly implicit solver\n\n";
    } else {
      cout << "   Nonlinear implicit solver\n\n";
    }
  }

  // Initialize vector data structures
  N = (udata->nxl)*(udata->nyl);
  Ntot = nx*ny;
  y = N_VNew_Parallel(udata->comm, N, Ntot);         // Create parallel vector for solution
  if (check_flag((void *) y, "N_VNew_Parallel", 0)) return 1;
  N_VConst(ZERO, y);                                 // Set initial conditions
  udata->h = N_VNew_Parallel(udata->comm, N, Ntot);  // Create vector for heat source
  if (check_flag((void *) udata->h, "N_VNew_Parallel", 0)) return 1;
  udata->d = N_VNew_Parallel(udata->comm, N, Ntot);  // Create vector for Jacobian diagonal
  if (check_flag((void *) udata->d, "N_VNew_Parallel", 0)) return 1;

  // Initialize linear solver data structures
  LSa = SUNLinSol_PCG(y, 1, 20);
  if (check_flag((void *) LSa, "SUNLinSol_PCG", 0)) return 1;
  LSm = SUNLinSol_PCG(y, 1, 20);
  if (check_flag((void *) LSm, "SUNLinSol_PCG", 0)) return 1;

  // fill in the heat source array
  data = N_VGetArrayPointer(udata->h);
  for (j=0; j<udata->nyl; j++)
    for (i=0; i<udata->nxl; i++)
      data[IDX(i,j,udata->nxl)] = sin(PI*(udata->is+i)*udata->dx)
                                * sin(TWO*PI*(udata->js+j)*udata->dy);

  /* Call ARKStepCreate and MRIStepCreate to initialize the timesteppers */
  arkstep_mem = ARKStepCreate(NULL, f, T0, y);
  if (check_flag((void *) arkstep_mem, "ARKStepCreate", 0)) return 1;
  inner_mem = ARKStepCreate(f0, NULL, T0, y);
  if (check_flag((void *) inner_mem, "ARKStepCreate", 0)) return 1;
  mristep_mem = MRIStepCreate(f, T0, y, MRISTEP_ARKSTEP, inner_mem);
  if (check_flag((void *) mristep_mem, "MRIStepCreate", 0)) return 1;

  // Create solve-decoupled DIRK2 (trapezoidal) Butcher table
  B = ARKodeButcherTable_Alloc(3, SUNFALSE);
  if (check_flag((void *)B, "ARKodeButcherTable_Alloc", 0)) return 1;
  B->A[1][0] = ONE;
  B->A[2][0] = RCONST(0.5);
  B->A[2][2] = RCONST(0.5);
  B->b[0] = RCONST(0.5);
  B->b[2] = RCONST(0.5);
  B->c[1] = ONE;
  B->c[2] = ONE;
  B->q=2;

  // Set routines
  flag = ARKStepSetUserData(arkstep_mem, (void *) udata);      // Pass udata to user functions
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
  flag = ARKStepSetNonlinConvCoef(arkstep_mem, RCONST(1.e-7)); // Update solver convergence coeff.
  if (check_flag(&flag, "ARKStepSetNonlinConvCoef", 1)) return 1;
  flag = ARKStepSStolerances(arkstep_mem, rtol, atol);         // Specify tolerances
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;
  flag = ARKStepSetFixedStep(arkstep_mem, Tf/Nt);              // Specify fixed time step size
  if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;
  flag = ARKStepSetTables(arkstep_mem, 2, 0, B, NULL);         // Specify Butcher table
  if (check_flag(&flag, "ARKStepSetTables", 1)) return 1;
  flag = ARKStepSetMaxNumSteps(arkstep_mem, 2*Nt);             // Increase num internal steps
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;
  flag = MRIStepSetUserData(mristep_mem, (void *) udata);      // Pass udata to user functions
  if (check_flag(&flag, "MRIStepSetUserData", 1)) return 1;
  flag = MRIStepSetNonlinConvCoef(mristep_mem, RCONST(1.e-7)); // Update solver convergence coeff.
  if (check_flag(&flag, "MRIStepSetNonlinConvCoef", 1)) return 1;
  flag = MRIStepSStolerances(mristep_mem, rtol, atol);         // Specify tolerances
  if (check_flag(&flag, "MRIStepSStolerances", 1)) return 1;
  flag = MRIStepSetFixedStep(mristep_mem, Tf/Nt);              // Specify fixed time step sizes
  if (check_flag(&flag, "MRIStepSetFixedStep", 1)) return 1;
  flag = ARKStepSetFixedStep(inner_mem, Tf/Nt/10);
  if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;
  flag = MRIStepSetTable(mristep_mem, 2, B);                   // Specify Butcher table
  if (check_flag(&flag, "MRIStepSetTable", 1)) return 1;
  flag = MRIStepSetMaxNumSteps(mristep_mem, 2*Nt);             // Increase num internal steps
  if (check_flag(&flag, "MRIStepSetMaxNumSteps", 1)) return 1;

  // Linear solver interface
  flag = ARKStepSetLinearSolver(arkstep_mem, LSa, NULL);      // Attach linear solver
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
  flag = ARKStepSetPreconditioner(arkstep_mem, PSet, PSol);   // Specify the Preconditoner
  if (check_flag(&flag, "ARKStepSetPreconditioner", 1)) return 1;
  flag = MRIStepSetLinearSolver(mristep_mem, LSm, NULL);      // Attach linear solver
  if (check_flag(&flag, "MRIStepSetLinearSolver", 1)) return 1;
  flag = MRIStepSetPreconditioner(mristep_mem, PSet, PSol);   // Specify the Preconditoner
  if (check_flag(&flag, "MRIStepSetPreconditioner", 1)) return 1;

  // Optionally specify linearly implicit RHS, with non-time-dependent preconditioner
  if (linear) {
    flag = ARKStepSetLinear(arkstep_mem, 0);
    if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;
    flag = MRIStepSetLinear(mristep_mem, 0);
    if (check_flag(&flag, "MRIStepSetLinear", 1)) return 1;
  }

  // First call ARKStep to evolve the full problem, and print results
  t = T0;
  N_VConst(ZERO, y);
  flag = ARKStepEvolve(arkstep_mem, Tf, y, &t, ARK_NORMAL);
  if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;
  flag = ARKStepGetNumSteps(arkstep_mem, &ark_nst);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return 1;
  flag = ARKStepGetNumRhsEvals(arkstep_mem, &ark_nfe, &ark_nfi);
  if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return 1;
  flag = ARKStepGetNumLinSolvSetups(arkstep_mem, &ark_nsetups);
  if (check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1)) return 1;
  flag = ARKStepGetNumNonlinSolvIters(arkstep_mem, &ark_nni);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return 1;
  flag = ARKStepGetNumNonlinSolvConvFails(arkstep_mem, &ark_ncfn);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1)) return 1;
  flag = ARKStepGetNumLinIters(arkstep_mem, &ark_nli);
  if (check_flag(&flag, "ARKStepGetNumLinIters", 1)) return 1;
  flag = ARKStepGetNumJtimesEvals(arkstep_mem, &ark_nJv);
  if (check_flag(&flag, "ARKStepGetNumJtimesEvals", 1)) return 1;
  flag = ARKStepGetNumLinConvFails(arkstep_mem, &ark_nlcf);
  if (check_flag(&flag, "ARKStepGetNumLinConvFails", 1)) return 1;
  flag = ARKStepGetNumPrecEvals(arkstep_mem, &ark_npe);
  if (check_flag(&flag, "ARKStepGetNumPrecEvals", 1)) return 1;
  flag = ARKStepGetNumPrecSolves(arkstep_mem, &ark_nps);
  if (check_flag(&flag, "ARKStepGetNumPrecSolves", 1)) return 1;
  if (outproc) {
    cout << "\nARKStep Solver Statistics:\n";
    cout << "   Internal solver steps = " << ark_nst << "\n";
    cout << "   Total RHS evals:  Fe = " << ark_nfe << ",  Fi = " << ark_nfi << "\n";
    cout << "   Total linear solver setups = " << ark_nsetups << "\n";
    cout << "   Total linear iterations = " << ark_nli << "\n";
    cout << "   Total number of Jacobian-vector products = " << ark_nJv << "\n";
    cout << "   Total number of Preconditioner setups = " << ark_npe << "\n";
    cout << "   Total number of Preconditioner solves = " << ark_nps << "\n";
    cout << "   Total number of linear solver convergence failures = " << ark_nlcf << "\n";
    cout << "   Total number of Newton iterations = " << ark_nni << "\n";
    cout << "   Total number of nonlinear solver convergence failures = " << ark_ncfn << "\n";
  }


  // Second call MRIStep to evolve the full problem, and print results
  t = T0;
  N_VConst(ZERO, y);
  flag = MRIStepEvolve(mristep_mem, Tf, y, &t, ARK_NORMAL);
  if (check_flag(&flag, "MRIStepEvolve", 1)) return 1;
  flag = MRIStepGetNumSteps(mristep_mem, &mri_nst);
  if (check_flag(&flag, "MRIStepGetNumSteps", 1)) return 1;
  flag = MRIStepGetNumRhsEvals(mristep_mem, &mri_nfs);
  if (check_flag(&flag, "MRIStepGetNumRhsEvals", 1)) return 1;
  flag = MRIStepGetNumLinSolvSetups(mristep_mem, &mri_nsetups);
  if (check_flag(&flag, "MRIStepGetNumLinSolvSetups", 1)) return 1;
  flag = MRIStepGetNumNonlinSolvIters(mristep_mem, &mri_nni);
  if (check_flag(&flag, "MRIStepGetNumNonlinSolvIters", 1)) return 1;
  flag = MRIStepGetNumNonlinSolvConvFails(mristep_mem, &mri_ncfn);
  if (check_flag(&flag, "MRIStepGetNumNonlinSolvConvFails", 1)) return 1;
  flag = MRIStepGetNumLinIters(mristep_mem, &mri_nli);
  if (check_flag(&flag, "MRIStepGetNumLinIters", 1)) return 1;
  flag = MRIStepGetNumJtimesEvals(mristep_mem, &mri_nJv);
  if (check_flag(&flag, "MRIStepGetNumJtimesEvals", 1)) return 1;
  flag = MRIStepGetNumLinConvFails(mristep_mem, &mri_nlcf);
  if (check_flag(&flag, "MRIStepGetNumLinConvFails", 1)) return 1;
  flag = MRIStepGetNumPrecEvals(mristep_mem, &mri_npe);
  if (check_flag(&flag, "MRIStepGetNumPrecEvals", 1)) return 1;
  flag = MRIStepGetNumPrecSolves(mristep_mem, &mri_nps);
  if (check_flag(&flag, "MRIStepGetNumPrecSolves", 1)) return 1;
  if (outproc) {
    cout << "\nMRIStep Solver Statistics:\n";
    cout << "   Internal solver steps = " << mri_nst << "\n";
    cout << "   Total RHS evals:  Fe = " << mri_nfs << "\n";
    cout << "   Total linear solver setups = " << mri_nsetups << "\n";
    cout << "   Total linear iterations = " << mri_nli << "\n";
    cout << "   Total number of Jacobian-vector products = " << mri_nJv << "\n";
    cout << "   Total number of Preconditioner setups = " << mri_npe << "\n";
    cout << "   Total number of Preconditioner solves = " << mri_nps << "\n";
    cout << "   Total number of linear solver convergence failures = " << mri_nlcf << "\n";
    cout << "   Total number of Newton iterations = " << mri_nni << "\n";
    cout << "   Total number of nonlinear solver convergence failures = " << mri_ncfn << "\n";
  }


  // Compare solver statistics
  numfails = 0;
  if (outproc)
    cout << "\nComparing Solver Statistics:\n";
  if (ark_nst != mri_nst) {
    numfails += 1;
    if (outproc)
      cout << "  Internal solver steps error: " << ark_nst << " vs " << mri_nst << "\n";
  }
  if (ark_nfi != (mri_nfs+mri_nst)) {
    numfails += 1;
    if (outproc)
      cout << "  RHS evals error: " << ark_nfi << " vs " << mri_nfs << "\n";
  }
  if (ark_nsetups != mri_nsetups) {
    numfails += 1;
    if (outproc)
      cout << "  Linear solver setups error: " << ark_nsetups << " vs " << mri_nsetups << "\n";
  }
  if (ark_nli < mri_nli) {
    numfails += 1;
    if (outproc)
      cout << "  Linear iterations error: " << ark_nli << " vs " << mri_nli << "\n";
  }
  if (ark_nJv < mri_nJv) {
    numfails += 1;
    if (outproc)
      cout << "  Jacobian-vector products error: " << ark_nJv << " vs " << mri_nJv << "\n";
  }
  if (ark_nps < mri_nps) {
    numfails += 1;
    if (outproc)
      cout << "  Preconditioner solves error: " << ark_nps << " vs " << mri_nps << "\n";
  }
  if (ark_nlcf != mri_nlcf) {
    numfails += 1;
    if (outproc)
      cout << "  Linear convergence failures error: " << ark_nlcf << " vs " << mri_nlcf << "\n";
  }
  if (ark_nni != mri_nni) {
    numfails += 1;
    if (outproc)
      cout << "  Newton iterations error: " << ark_nni << " vs " << mri_nni << "\n";
  }
  if (ark_ncfn != mri_ncfn) {
    numfails += 1;
    if (outproc)
      cout << "  Nonlinear convergence failures error: " << ark_ncfn << " vs " << mri_ncfn << "\n";
  }
  if (outproc) {
    if (numfails) {
      cout << "Failed " << numfails << " tests\n";
    } else {
      cout << "All tests pass!\n";
    }
  }

  // Clean up and return with successful completion
  ARKodeButcherTable_Free(B);  // Free Butcher table
  ARKStepFree(&arkstep_mem);   // Free integrator memory
  MRIStepFree(&mristep_mem);
  ARKStepFree(&inner_mem);
  SUNLinSolFree(LSa);          // Free linear solver
  SUNLinSolFree(LSm);
  N_VDestroy(y);               // Free vectors
  N_VDestroy(udata->h);
  N_VDestroy(udata->d);
  FreeUserData(udata);         // Free user data
  delete udata;
  flag = MPI_Finalize();       // Finalize MPI
  return (numfails);
}


/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

// f routine to compute the ODE RHS function f(t,y).
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(ZERO, ydot);                          // Initialize ydot to zero
  UserData *udata = (UserData *) user_data;      // access problem data
  sunindextype nxl = udata->nxl;                     // set variable shortcuts
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
  return 0;                                      // Return with success
}

// f0 routine to compute a zero-valued ODE RHS function f(t,y).
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  // Initialize ydot to zero and return
  N_VConst(ZERO, ydot);
  return 0;
}

// Preconditioner setup routine (fills inverse of Jacobian diagonal)
static int PSet(realtype t, N_Vector y, N_Vector fy, booleantype jok,
                booleantype *jcurPtr, realtype gamma, void *user_data)
{
  UserData *udata = (UserData *) user_data;      // variable shortcuts
  realtype kx = udata->kx;
  realtype ky = udata->ky;
  realtype dx = udata->dx;
  realtype dy = udata->dy;
  realtype *diag = N_VGetArrayPointer(udata->d);  // access data arrays
  if (check_flag((void *) diag, "N_VGetArrayPointer", 0)) return -1;

  // set all entries of d to the diagonal values of interior
  // (since boundary RHS is 0, set boundary diagonals to the same)
  realtype c = ONE + gamma*TWO*(kx/dx/dx + ky/dy/dy);
  N_VConst(c, udata->d);
  N_VInv(udata->d, udata->d);  // invert diagonal
  return 0;                    // Return with success
}

// Preconditioner solve routine
static int PSol(realtype t, N_Vector y, N_Vector fy, N_Vector r,
                N_Vector z, realtype gamma, realtype delta, int lr,
                void *user_data)
{
  UserData *udata = (UserData *) user_data;  // access user_data structure
  N_VProd(r, udata->d, z);                   // perform Jacobi iteration
  return 0;                                  // Return with success
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
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
  sunindextype nyl = udata->nyl;
  sunindextype nxl = udata->nxl;

  // access data array
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;

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
    ierr = MPI_Irecv(udata->Wrecv, (int) udata->nyl, MPI_SUNREALTYPE, ipW,
                     MPI_ANY_TAG, udata->comm, &reqRW);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    ierr = MPI_Irecv(udata->Erecv, (int) udata->nyl, MPI_SUNREALTYPE, ipE,
                     MPI_ANY_TAG, udata->comm, &reqRE);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    ierr = MPI_Irecv(udata->Srecv, (int) udata->nxl, MPI_SUNREALTYPE, ipS,
                     MPI_ANY_TAG, udata->comm, &reqRS);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    ierr = MPI_Irecv(udata->Nrecv, (int) udata->nxl, MPI_SUNREALTYPE, ipN,
                     MPI_ANY_TAG, udata->comm, &reqRN);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }

  // send data
  if (!udata->HaveBdry[0][0]) {
    for (i=0; i<nyl; i++)  udata->Wsend[i] = Y[IDX(0,i,nxl)];
    ierr = MPI_Isend(udata->Wsend, (int) udata->nyl, MPI_SUNREALTYPE, ipW, 0,
                     udata->comm, &reqSW);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1]) {
    for (i=0; i<nyl; i++)  udata->Esend[i] = Y[IDX(nxl-1,i,nxl)];
    ierr = MPI_Isend(udata->Esend, (int) udata->nyl, MPI_SUNREALTYPE, ipE, 1,
                     udata->comm, &reqSE);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0]) {
    for (i=0; i<nxl; i++)  udata->Ssend[i] = Y[IDX(i,0,nxl)];
    ierr = MPI_Isend(udata->Ssend, (int) udata->nxl, MPI_SUNREALTYPE, ipS, 2,
                     udata->comm, &reqSS);
    if (ierr != MPI_SUCCESS) {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1]) {
    for (i=0; i<nxl; i++)  udata->Nsend[i] = Y[IDX(i,nyl-1,nxl)];
    ierr = MPI_Isend(udata->Nsend, (int) udata->nxl, MPI_SUNREALTYPE, ipN, 3,
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
  udata->d = NULL;
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

  if (udata->comm != MPI_COMM_WORLD)
    MPI_Comm_free(&(udata->comm));

  return 0;     // return with success flag
}


//---- end of file ----
