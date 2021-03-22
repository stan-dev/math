/*-----------------------------------------------------------------
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
 * solver and the fixed-point solver based on a 0/1 command-line
 * argument (1 => fixed-point).
 *-----------------------------------------------------------------*/

// Header files
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sundials/sundials_types.h>

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

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

// User-supplied Functions Called by the Solver
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private function to perform matrix-matrix product
static int dense_MM(SUNMatrix A, SUNMatrix B, SUNMatrix C);

// Private function to check function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

// Main Program
int main(int argc, char* argv[])
{
  // general problem parameters
  realtype T0 = RCONST(0.0);         // initial time
  realtype Tf = RCONST(0.05);        // final time
  int Nt = 1000;                     // total number of internal steps
  sunindextype NEQ = 3;              // number of dependent vars.
  realtype reltol = RCONST(1.0e-6);  // tolerances
  realtype abstol = RCONST(1.0e-10);
  realtype lamda  = RCONST(-100.0);  // stiffness parameter

  // general problem variables
  int flag;                       // reusable error-checking flag
  N_Vector y = NULL;              // empty vector for storing solution
  SUNMatrix Aa = NULL;            // empty dense matrices for solvers
  SUNMatrix Am = NULL;
  SUNNonlinearSolver NLSa = NULL; // empty nonlinear solvers
  SUNNonlinearSolver NLSm = NULL;
  SUNLinearSolver LSa = NULL;     // empty dense linear solvers
  SUNLinearSolver LSm = NULL;
  void *arkstep_mem = NULL;       // empty ARKStep memory structure
  void *mristep_mem = NULL;       // empty MRIStep memory structure
  void *inner_mem = NULL;         // empty inner ARKStep memory structure
  ARKodeButcherTable B = NULL;
  int numfails;
  booleantype fixedpoint;
  realtype t;
  long int ark_nst, ark_nfe, ark_nfi, ark_nsetups, ark_nje, ark_nfeLS, ark_nni, ark_ncfn;
  long int mri_nst, mri_nfs, mri_nsetups, mri_nje, mri_nfeLS, mri_nni, mri_ncfn;

  // if an argument supplied, set fixedpoint (otherwise use SUNFALSE)
  fixedpoint = SUNFALSE;
  if (argc > 1)  fixedpoint = stoi(argv[1], NULL);

  // Initial problem output
  cout << "\nAnalytical ODE test problem:\n";
  cout << "   lamda  = " << lamda << "\n";
  cout << "   reltol = " << reltol << "\n";
  cout << "   abstol = " << abstol << "\n\n";
  if (fixedpoint) {
    cout << "   Fixed-point nonlinear solver\n";
  } else {
    cout << "   Newton nonlinear solver\n";
  }

  // Initialize vector data structure and specify initial condition
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  N_VConst(ONE, y);

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
  flag = ARKStepSetUserData(arkstep_mem, (void *) &lamda);   // Pass lamda to user functions
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
  flag = ARKStepSStolerances(arkstep_mem, reltol, abstol);   // Specify tolerances
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;
  flag = ARKStepSetFixedStep(arkstep_mem, Tf/Nt);            // Specify fixed time step size
  if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;
  flag = ARKStepSetTables(arkstep_mem, 2, 0, B, NULL);       // Specify Butcher table
  if (check_flag(&flag, "ARKStepSetTables", 1)) return 1;
  flag = ARKStepSetMaxNumSteps(arkstep_mem, 2*Nt);           // Increase num internal steps
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;
  flag = MRIStepSetUserData(mristep_mem, (void *) &lamda);   // Pass lamda to user functions
  if (check_flag(&flag, "MRIStepSetUserData", 1)) return 1;
  flag = MRIStepSStolerances(mristep_mem, reltol, abstol);   // Specify tolerances
  if (check_flag(&flag, "MRIStepSStolerances", 1)) return 1;
  flag = MRIStepSetFixedStep(mristep_mem, Tf/Nt);            // Specify fixed time step sizes
  if (check_flag(&flag, "MRIStepSetFixedStep", 1)) return 1;
  flag = ARKStepSetFixedStep(inner_mem, Tf/Nt/10);
  if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;
  flag = MRIStepSetTable(mristep_mem, 2, B);                 // Specify Butcher table
  if (check_flag(&flag, "MRIStepSetTable", 1)) return 1;
  flag = MRIStepSetMaxNumSteps(mristep_mem, 2*Nt);           // Increase num internal steps
  if (check_flag(&flag, "MRIStepSetMaxNumSteps", 1)) return 1;

  // Initialize implicit solver data structures
  if (fixedpoint) {

    // Initialize fixed-point solvers and attach to integrators
    NLSa = SUNNonlinSol_FixedPoint(y, 50);
    if (check_flag((void *)NLSa, "SUNNonlinSol_FixedPoint", 0)) return 1;
    flag = ARKStepSetNonlinearSolver(arkstep_mem, NLSa);
    if (check_flag(&flag, "ARKStepSetNonlinearSolver", 1)) return 1;
    NLSm = SUNNonlinSol_FixedPoint(y, 50);
    if (check_flag((void *)NLSm, "SUNNonlinSol_FixedPoint", 0)) return 1;
    flag = MRIStepSetNonlinearSolver(mristep_mem, NLSm);
    if (check_flag(&flag, "MRIStepSetNonlinearSolver", 1)) return 1;

  } else {

    // Initialize dense matrix data structures and solvers
    Aa = SUNDenseMatrix(NEQ, NEQ);
    if (check_flag((void *)Aa, "SUNDenseMatrix", 0)) return 1;
    LSa = SUNLinSol_Dense(y, Aa);
    if (check_flag((void *)LSa, "SUNLinSol_Dense", 0)) return 1;
    Am = SUNDenseMatrix(NEQ, NEQ);
    if (check_flag((void *)Am, "SUNDenseMatrix", 0)) return 1;
    LSm = SUNLinSol_Dense(y, Am);
    if (check_flag((void *)LSm, "SUNLinSol_Dense", 0)) return 1;

    // Linear solver interface
    flag = ARKStepSetLinearSolver(arkstep_mem, LSa, Aa);
    if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
    flag = ARKStepSetJacFn(arkstep_mem, Jac);
    if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;
    flag = MRIStepSetLinearSolver(mristep_mem, LSm, Am);
    if (check_flag(&flag, "MRIStepSetLinearSolver", 1)) return 1;
    flag = MRIStepSetJacFn(mristep_mem, Jac);
    if (check_flag(&flag, "MRIStepSetJacFn", 1)) return 1;

    // Specify linearly implicit RHS, with non-time-dependent Jacobian
    flag = ARKStepSetLinear(arkstep_mem, 0);
    if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;
    flag = MRIStepSetLinear(mristep_mem, 0);
    if (check_flag(&flag, "MRIStepSetLinear", 1)) return 1;

  }


  // First call ARKStep to evolve the full problem, and print results
  t = T0;
  N_VConst(ONE, y);
  flag = ARKStepEvolve(arkstep_mem, Tf, y, &t, ARK_NORMAL);
  if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;
  flag = ARKStepGetNumSteps(arkstep_mem, &ark_nst);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return 1;
  flag = ARKStepGetNumRhsEvals(arkstep_mem, &ark_nfe, &ark_nfi);
  if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return 1;
  flag = ARKStepGetNumNonlinSolvIters(arkstep_mem, &ark_nni);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return 1;
  flag = ARKStepGetNumNonlinSolvConvFails(arkstep_mem, &ark_ncfn);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1)) return 1;
  if (!fixedpoint) {
    flag = ARKStepGetNumLinSolvSetups(arkstep_mem, &ark_nsetups);
    if (check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1)) return 1;
    flag = ARKStepGetNumJacEvals(arkstep_mem, &ark_nje);
    if (check_flag(&flag, "ARKStepGetNumJacEvals", 1)) return 1;
    flag = ARKStepGetNumLinRhsEvals(arkstep_mem, &ark_nfeLS);
    check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1);
  }
  cout << "\nARKStep Solver Statistics:\n";
  cout << "   Internal solver steps = " << ark_nst << "\n";
  cout << "   Total RHS evals:  Fe = " << ark_nfe << ",  Fi = " << ark_nfi << "\n";
  cout << "   Total number of nonlinear iterations = " << ark_nni << "\n";
  cout << "   Total number of nonlinear solver convergence failures = " << ark_ncfn << "\n";
  if (!fixedpoint) {
    cout << "   Total linear solver setups = " << ark_nsetups << "\n";
    cout << "   Total RHS evals for setting up the linear system = " << ark_nfeLS << "\n";
    cout << "   Total number of Jacobian evaluations = " << ark_nje << "\n";
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
  flag = MRIStepGetNumNonlinSolvIters(mristep_mem, &mri_nni);
  if (check_flag(&flag, "MRIStepGetNumNonlinSolvIters", 1)) return 1;
  flag = MRIStepGetNumNonlinSolvConvFails(mristep_mem, &mri_ncfn);
  if (check_flag(&flag, "MRIStepGetNumNonlinSolvConvFails", 1)) return 1;
  if (!fixedpoint) {
    flag = MRIStepGetNumLinSolvSetups(mristep_mem, &mri_nsetups);
    if (check_flag(&flag, "MRIStepGetNumLinSolvSetups", 1)) return 1;
    flag = MRIStepGetNumJacEvals(mristep_mem, &mri_nje);
    if (check_flag(&flag, "MRIStepGetNumJacEvals", 1)) return 1;
    flag = MRIStepGetNumLinRhsEvals(mristep_mem, &mri_nfeLS);
    check_flag(&flag, "MRIStepGetNumLinRhsEvals", 1);
  }
  cout << "\nMRIStep Solver Statistics:\n";
  cout << "   Internal solver steps = " << mri_nst << "\n";
  cout << "   Total RHS evals:  Fs = " << mri_nfs << "\n";
  cout << "   Total number of nonlinear iterations = " << mri_nni << "\n";
  cout << "   Total number of nonlinear solver convergence failures = " << mri_ncfn << "\n";
  if (!fixedpoint) {
    cout << "   Total linear solver setups = " << mri_nsetups << "\n";
    cout << "   Total RHS evals for setting up the linear system = " << mri_nfeLS << "\n";
    cout << "   Total number of Jacobian evaluations = " << mri_nje << "\n";
  }


  // Compare solver statistics
  numfails = 0;
  cout << "\nComparing Solver Statistics:\n";
  if (ark_nst != mri_nst) {
    numfails += 1;
    cout << "  Internal solver steps error: " << ark_nst << " vs " << mri_nst << "\n";
  }
  if (ark_nfi != (mri_nfs+mri_nst)) {
    numfails += 1;
    cout << "  RHS evals error: " << ark_nfi << " vs " << mri_nfs << "\n";
  }
  if (ark_nni != mri_nni) {
    numfails += 1;
    cout << "  Nonlinear iterations error: " << ark_nni << " vs " << mri_nni << "\n";
  }
  if (ark_ncfn != mri_ncfn) {
    numfails += 1;
    cout << "  Nonlinear convergence failures error: " << ark_ncfn << " vs " << mri_ncfn << "\n";
  }
  if (!fixedpoint) {
    if (ark_nsetups != mri_nsetups) {
      numfails += 1;
      cout << "  Linear solver setups error: " << ark_nsetups << " vs " << mri_nsetups << "\n";
    }
    if (ark_nfeLS != mri_nfeLS) {
      numfails += 1;
      cout << "  RHS evals for LS error: " << ark_nfeLS << " vs " << mri_nfeLS << "\n";
    }
    if (ark_nje != mri_nje) {
      numfails += 1;
      cout << "  Jacobian evaluations error: " << ark_nje << " vs " << mri_nje << "\n";
    }
  }
  if (numfails) {
    cout << "Failed " << numfails << " tests\n";
  } else {
    cout << "All tests pass!\n";
  }

  // Clean up and return with successful completion
  ARKodeButcherTable_Free(B);  // Free Butcher table
  ARKStepFree(&arkstep_mem);   // Free integrator memory
  MRIStepFree(&mristep_mem);
  ARKStepFree(&inner_mem);
  if (fixedpoint) {
    SUNNonlinSolFree(NLSa);    // Free nonlinear solvers
    SUNNonlinSolFree(NLSm);
  } else {
    SUNLinSolFree(LSa);        // Free linear solvers
    SUNLinSolFree(LSm);
    SUNMatDestroy(Aa);         // Free A matrices
    SUNMatDestroy(Am);
  }
  N_VDestroy(y);               // Free y vector
  return (numfails);
}


/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

// f routine to compute the ODE RHS function f(t,y).
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  realtype y0 = NV_Ith_S(y,0);                // access current solution values
  realtype y1 = NV_Ith_S(y,1);
  realtype y2 = NV_Ith_S(y,2);
  realtype yd0, yd1, yd2;

  // fill in the RHS function: f(t,y) = V*D*Vi*y
  yd0 = RCONST(0.25)*(RCONST(5.0)*y0 + RCONST(1.0)*y1 - RCONST(3.0)*y2);  // yd = Vi*y
  yd1 = RCONST(0.25)*(RCONST(2.0)*y0 + RCONST(2.0)*y1 - RCONST(2.0)*y2);
  yd2 = RCONST(0.25)*(RCONST(1.0)*y0 + RCONST(1.0)*y1 + RCONST(1.0)*y2);
  y0  = -RCONST(0.5)*yd0;                                                 //  y = D*yd
  y1  = -RCONST(0.1)*yd1;
  y2  =  lam*yd2;
  yd0 =  RCONST(1.0)*y0 - RCONST(1.0)*y1 + RCONST(1.0)*y2;                // yd = V*y
  yd1 = -RCONST(1.0)*y0 + RCONST(2.0)*y1 + RCONST(1.0)*y2;
  yd2 =  RCONST(0.0)*y0 - RCONST(1.0)*y1 + RCONST(2.0)*y2;
  NV_Ith_S(ydot,0) = yd0;
  NV_Ith_S(ydot,1) = yd1;
  NV_Ith_S(ydot,2) = yd2;

  return 0;                                  // Return with success
}

// f0 routine to compute a zero-valued ODE RHS function f(t,y).
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  // Initialize ydot to zero and return
  N_VConst(ZERO, ydot);
  return 0;
}

// Jacobian routine to compute J(t,y) = df/dy.
static int Jac(realtype t, N_Vector y, N_Vector fy,
               SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  SUNMatrix V  = SUNDenseMatrix(3,3);          // create temporary SUNMatrix objects
  SUNMatrix D  = SUNDenseMatrix(3,3);          // create temporary SUNMatrix objects
  SUNMatrix Vi = SUNDenseMatrix(3,3);          // create temporary SUNMatrix objects

  SUNMatZero(V);        // initialize temporary matrices to zero
  SUNMatZero(D);        // (not technically required)
  SUNMatZero(Vi);

  // Fill in temporary matrices:
  //    V = [1 -1 1; -1 2 1; 0 -1 2]
  SM_ELEMENT_D(V,0,0) =  RCONST(1.0);
  SM_ELEMENT_D(V,0,1) = -RCONST(1.0);
  SM_ELEMENT_D(V,0,2) =  RCONST(1.0);
  SM_ELEMENT_D(V,1,0) = -RCONST(1.0);
  SM_ELEMENT_D(V,1,1) =  RCONST(2.0);
  SM_ELEMENT_D(V,1,2) =  RCONST(1.0);
  SM_ELEMENT_D(V,2,0) =  RCONST(0.0);
  SM_ELEMENT_D(V,2,1) = -RCONST(1.0);
  SM_ELEMENT_D(V,2,2) =  RCONST(2.0);

  //    Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
  SM_ELEMENT_D(Vi,0,0) =  RCONST(0.25)*RCONST(5.0);
  SM_ELEMENT_D(Vi,0,1) =  RCONST(0.25)*RCONST(1.0);
  SM_ELEMENT_D(Vi,0,2) = -RCONST(0.25)*RCONST(3.0);
  SM_ELEMENT_D(Vi,1,0) =  RCONST(0.25)*RCONST(2.0);
  SM_ELEMENT_D(Vi,1,1) =  RCONST(0.25)*RCONST(2.0);
  SM_ELEMENT_D(Vi,1,2) = -RCONST(0.25)*RCONST(2.0);
  SM_ELEMENT_D(Vi,2,0) =  RCONST(0.25)*RCONST(1.0);
  SM_ELEMENT_D(Vi,2,1) =  RCONST(0.25)*RCONST(1.0);
  SM_ELEMENT_D(Vi,2,2) =  RCONST(0.25)*RCONST(1.0);

  //    D = [-0.5 0 0; 0 -0.1 0; 0 0 lam]
  SM_ELEMENT_D(D,0,0) = -RCONST(0.5);
  SM_ELEMENT_D(D,1,1) = -RCONST(0.1);
  SM_ELEMENT_D(D,2,2) = lam;

  // Compute J = V*D*Vi
  if (dense_MM(D,Vi,J) != 0) {     // J = D*Vi
    cerr << "matmul error\n";
    return 1;
  }
  if (dense_MM(V,J,D) != 0) {      // D = V*J [= V*D*Vi]
    cerr << "matmul error\n";
    return 1;
  }
  SUNMatCopy(D, J);

  SUNMatDestroy(V);                // Free V matrix
  SUNMatDestroy(D);                // Free D matrix
  SUNMatDestroy(Vi);               // Free Vi matrix

  return 0;                        // Return with success
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

// SUNDenseMatrix matrix-multiply utility routine: C = A*B
static int dense_MM(SUNMatrix A, SUNMatrix B, SUNMatrix C)
{
  // check for legal dimensions
  if ( (SUNDenseMatrix_Columns(A) != SUNDenseMatrix_Rows(B)) ||
       (SUNDenseMatrix_Rows(C) != SUNDenseMatrix_Rows(A)) ||
       (SUNDenseMatrix_Columns(C) != SUNDenseMatrix_Columns(B)) ) {
    cerr << "\n matmul error: dimension mismatch\n\n";
    return 1;
  }

  realtype **adata = SUNDenseMatrix_Cols(A);     // access data and extents
  realtype **bdata = SUNDenseMatrix_Cols(B);
  realtype **cdata = SUNDenseMatrix_Cols(C);
  sunindextype m = SUNDenseMatrix_Rows(C);
  sunindextype n = SUNDenseMatrix_Columns(C);
  sunindextype l = SUNDenseMatrix_Columns(A);
  sunindextype i, j, k;
  SUNMatZero(C);                                 // initialize output

  // perform multiply (not optimal, but fine for 3x3 matrices)
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      for (k=0; k<l; k++)
        cdata[i][j] += adata[i][k] * bdata[k][j];

  return 0;                       // Return with success
}

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



//---- end of file ----
