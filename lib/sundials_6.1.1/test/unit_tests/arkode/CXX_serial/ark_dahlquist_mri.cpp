/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * IMEX multirate Dahlquist problem:
 *
 * y' = lambda_e * y + lambda_i * y + lambda_f * y
 * ---------------------------------------------------------------------------*/

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

#include "arkode/arkode_mri_tables_impl.h"

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

// User data structure
struct UserData
{
  realtype lambda_e = RCONST(-1.0);
  realtype lambda_i = RCONST(-1.0);
  realtype lambda_f = RCONST(-1.0);
};

// User-supplied Functions called by the solver
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int ff(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Ji(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private function to check function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

// Test drivers
static int run_tests(MRISTEP_METHOD_TYPE type, realtype t0, int nsteps,
                     realtype hs, realtype hf, realtype reltol, realtype abstol,
                     UserData* udata, SUNContext ctx);


// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------


int main(int argc, char* argv[])
{
  // Initial time
  realtype t0 = RCONST(0.0);

  // Number of time steps
  int nsteps = 1;

  // Relative and absolute tolerances
  realtype reltol = RCONST(1.0e-4);
  realtype abstol = RCONST(1.0e-6);

  // Slow and fast step sizes
  realtype hs = RCONST(0.01);
  realtype hf = RCONST(0.01);

  // User data structure
  UserData udata;

  // Check for inputs
  if (argc > 1) udata.lambda_e = stod(argv[1]);
  if (argc > 2) udata.lambda_i = stod(argv[2]);
  if (argc > 3) udata.lambda_f = stod(argv[3]);
  if (argc > 4) hs = stod(argv[4]);
  if (argc > 5) hf = stod(argv[5]);
  if (argc > 5) nsteps = stoi(argv[6]);

  // Output problem setup
  cout << "\nDahlquist ODE test problem:\n";
  cout << "   lambda expl  = " << udata.lambda_e << "\n";
  cout << "   lambda impl  = " << udata.lambda_i << "\n";
  cout << "   lambda fast  = " << udata.lambda_f << "\n";
  cout << "   h slow       = " << hs << "\n";
  cout << "   h fast       = " << hf << "\n";
  cout << "   relative tol = " << reltol << "\n";
  cout << "   absolute tol = " << abstol << "\n";

  // Create SUNDIALS context
  sundials::Context sunctx;

  // Test methods
  int numfails = 0;

  numfails += run_tests(MRISTEP_EXPLICIT,
                        t0, nsteps, hs, hf, reltol, abstol, &udata, sunctx);

  numfails += run_tests(MRISTEP_IMPLICIT,
                        t0, nsteps, hs, hf, reltol, abstol, &udata, sunctx);

  numfails += run_tests(MRISTEP_IMEX,
                        t0, nsteps, hs, hf, reltol, abstol, &udata, sunctx);

  if (numfails)
  {
    cout << "\n\nFailed " << numfails << " tests!\n";
  }
  else
  {
    cout << "\n\nAll tests passed!\n";
  }

  // Return test status
  return numfails;
}


// -----------------------------------------------------------------------------
// Test drivers
// -----------------------------------------------------------------------------


int run_tests(MRISTEP_METHOD_TYPE type, realtype t0, int nsteps,
              realtype hs, realtype hf, realtype reltol, realtype abstol,
              UserData* udata, SUNContext sunctx)
{
  // Reusable error-checking flag
  int flag;

  // Test failure counter
  int numfails = 0;

  // Create initial condition vector
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;

  N_VConst(RCONST(1.0), y);

  // Create matrix and linear solver (if necessary)
  SUNMatrix       A  = nullptr;
  SUNLinearSolver LS = nullptr;

  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
  {
    // Initialize dense matrix data structures and solvers
    A = SUNDenseMatrix(1, 1, sunctx);
    if (check_flag((void *)A, "SUNDenseMatrix", 0)) return 1;

    LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_flag((void *)LS, "SUNLinSol_Dense", 0)) return 1;
  }

  // ----------------------
  // Create fast integrator
  // ----------------------

  // Create explicit fast integrator
  void* arkstep_mem = ARKStepCreate(ff, nullptr, t0, y, sunctx);
  if (check_flag((void *) arkstep_mem, "ARKStepCreate", 0)) return 1;

  // Set user data
  flag = ARKStepSetUserData(arkstep_mem, udata);
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

  // Specify tolerances
  flag = ARKStepSStolerances(arkstep_mem, reltol, abstol);
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  // Specify fixed time step size
  flag = ARKStepSetFixedStep(arkstep_mem, hf);
  if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;

  // Wrap ARKStep integrator as fast integrator object
  MRIStepInnerStepper inner_stepper = nullptr;
  flag = ARKStepCreateMRIStepInnerStepper(arkstep_mem, &inner_stepper);
  if (check_flag(&flag, "ARKStepCreateMRIStepInnerStepper", 1)) return 1;

  // ----------------------
  // Create slow integrator
  // ----------------------

  // Create slow integrator based on MRI type
  void* mristep_mem = nullptr;

  if (type == MRISTEP_EXPLICIT)
  {
    mristep_mem = MRIStepCreate(fe, nullptr, t0, y, inner_stepper, sunctx);
  }
  else if (type == MRISTEP_IMPLICIT)
  {
    mristep_mem = MRIStepCreate(nullptr, fi, t0, y, inner_stepper, sunctx);
  }
  else if (type == MRISTEP_IMEX)
  {
    mristep_mem = MRIStepCreate(fe, fi, t0, y, inner_stepper, sunctx);
  }
  else
  {
    return 1;
  }
  if (check_flag((void *) mristep_mem, "MRIStepCreate", 0)) return 1;

  // Set user data
  flag = MRIStepSetUserData(mristep_mem, udata);
  if (check_flag(&flag, "MRIStepSetUserData", 1)) return 1;

  // Specify tolerances
  flag = MRIStepSStolerances(mristep_mem, reltol, abstol);
  if (check_flag(&flag, "MRIStepSStolerances", 1)) return 1;

  // Specify fixed time step sizes
  flag = MRIStepSetFixedStep(mristep_mem, hs);
  if (check_flag(&flag, "MRIStepSetFixedStep", 1)) return 1;

  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
  {
    // Attach linear solver
    flag = MRIStepSetLinearSolver(mristep_mem, LS, A);
    if (check_flag(&flag, "MRIStepSetLinearSolver", 1)) return 1;

    // Set Jacobian function
    flag = MRIStepSetJacFn(mristep_mem, Ji);
    if (check_flag(&flag, "MRIStepSetJacFn", 1)) return 1;

    // Specify linearly implicit RHS, with non-time-dependent Jacobian
    flag = MRIStepSetLinear(mristep_mem, 0);
    if (check_flag(&flag, "MRIStepSetLinear", 1)) return 1;
  }

  // ------------------------------------
  // Evolve with various IMEX MRI methods
  // ------------------------------------

  // Methods to test (order most stages to least since reinit does not realloc)
  int   num_methods;
  ARKODE_MRITableID*  methods = nullptr;
  bool* stiffly_accurate = nullptr;

  if (type == MRISTEP_EXPLICIT)
  {
    cout << "\n=========================\n";
    cout << "Test explicit MRI methods\n";
    cout << "=========================\n";

    num_methods = 3;
    methods = new ARKODE_MRITableID[num_methods];

    methods[0] = ARKODE_MRI_GARK_ERK45a;
    methods[1] = ARKODE_MRI_GARK_ERK33a;
    methods[2] = ARKODE_MIS_KW3;
  }
  else if (type == MRISTEP_IMPLICIT)
  {
    cout << "\n=========================\n";
    cout << "Test implicit MRI methods\n";
    cout << "=========================\n";

    num_methods = 3;
    methods = new ARKODE_MRITableID[num_methods];
    stiffly_accurate = new bool[num_methods];

    methods[0] = ARKODE_MRI_GARK_ESDIRK46a;
    stiffly_accurate[0] = true;

    methods[1] = ARKODE_MRI_GARK_ESDIRK34a;
    stiffly_accurate[1] = true;

    methods[2] = ARKODE_MRI_GARK_IRK21a;
    stiffly_accurate[2] = true;
  }
  else if (type == MRISTEP_IMEX)
  {
    cout << "\n=====================\n";
    cout << "Test IMEX MRI methods\n";
    cout << "=====================\n";

    num_methods = 3;
    methods = new ARKODE_MRITableID[num_methods];
    stiffly_accurate = new bool[num_methods];

    methods[0] = ARKODE_IMEX_MRI_GARK4;
    stiffly_accurate[0] = false;

    methods[1] = ARKODE_IMEX_MRI_GARK3b;
    stiffly_accurate[1] = false;

    methods[2] = ARKODE_IMEX_MRI_GARK3a;
    stiffly_accurate[2] = false;
  }
  else
  {
    return 1;
  }

  for (int i = 0; i < num_methods; i++)
  {
    cout << "\nTesting method " << i << "\n";

    // -------------
    // Select method
    // -------------

    // Load method table
    MRIStepCoupling C = MRIStepCoupling_LoadTable(methods[i]);
    if (check_flag((void *)C, "MRIStepCoupling_LoadTable", 0)) return 1;

    MRIStepCoupling_Write(C, stdout);

    // Get the number of stored stages
    int* stage_map = new int[C->stages];
    int  nstages_stored;

    flag = mriStepCoupling_GetStageMap(C, stage_map, &nstages_stored);
    if (check_flag(&flag, "mriStepCoupling_GetStageMap", 1)) return 1;

    cout << "  Stored stages = " << nstages_stored << "\n";
    delete[] stage_map;

    // Set coupling table
    flag = MRIStepSetCoupling(mristep_mem, C);
    if (check_flag(&flag, "MRIStepSetCoupling", 1)) return 1;

    // -----------------
    // Output statistics
    // -----------------

    realtype t  = t0;
    realtype tf = nsteps * hs;

    for (int i = 0; i < nsteps; i++)
    {
      // Advance in time
      flag = MRIStepEvolve(mristep_mem, tf, y, &t, ARK_ONE_STEP);
      if (check_flag(&flag, "MRIStepEvolve", 1)) return 1;

      // Update output time
      tf += hs;
    }

    // -----------------
    // Output statistics
    // -----------------

    long int mri_nst, mri_nfse, mri_nfsi;      // integrator
    long int mri_nni, mri_ncfn;                // nonlinear solver
    long int mri_nsetups, mri_nje, mri_nfeLS;  // linear solver

    flag = MRIStepGetNumSteps(mristep_mem, &mri_nst);
    if (check_flag(&flag, "MRIStepGetNumSteps", 1)) return 1;

    flag = MRIStepGetNumRhsEvals(mristep_mem, &mri_nfse, &mri_nfsi);
    if (check_flag(&flag, "MRIStepGetNumRhsEvals", 1)) return 1;

    if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
    {
      flag = MRIStepGetNumNonlinSolvIters(mristep_mem, &mri_nni);
      if (check_flag(&flag, "MRIStepGetNumNonlinSolvIters", 1)) return 1;

      flag = MRIStepGetNumNonlinSolvConvFails(mristep_mem, &mri_ncfn);
      if (check_flag(&flag, "MRIStepGetNumNonlinSolvConvFails", 1)) return 1;

      flag = MRIStepGetNumLinSolvSetups(mristep_mem, &mri_nsetups);
      if (check_flag(&flag, "MRIStepGetNumLinSolvSetups", 1)) return 1;

      flag = MRIStepGetNumJacEvals(mristep_mem, &mri_nje);
      if (check_flag(&flag, "MRIStepGetNumJacEvals", 1)) return 1;

      flag = MRIStepGetNumLinRhsEvals(mristep_mem, &mri_nfeLS);
      check_flag(&flag, "MRIStepGetNumLinRhsEvals", 1);
    }

    realtype  pow   = udata->lambda_f;
    if (type == MRISTEP_EXPLICIT || type == MRISTEP_IMEX)
    {
      pow += udata->lambda_e;
    }
    if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
    {
      pow += udata->lambda_i;
    }
    realtype  ytrue = exp(pow * t);

    realtype* ydata = N_VGetArrayPointer(y);
    realtype  error = ytrue - ydata[0];

    cout << "\nMRIStep Statistics:\n";
    cout << "   Time        = " << t           << "\n";
    cout << "   y(t)        = " << ytrue       << "\n";
    cout << "   y_n         = " << ydata[0]    << "\n";
    cout << "   Error       = " << error       << "\n";
    cout << "   Steps       = " << mri_nst     << "\n";
    cout << "   Fe evals    = " << mri_nfse    << "\n";
    cout << "   Fi evals    = " << mri_nfsi    << "\n";

    if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
    {
      cout << "   NLS iters   = " << mri_nni     << "\n";
      cout << "   NLS fails   = " << mri_ncfn    << "\n";
      cout << "   LS setups   = " << mri_nsetups << "\n";
      cout << "   LS Fi evals = " << mri_nfeLS   << "\n";
      cout << "   Ji evals    = " << mri_nje     << "\n";
    }

    // ----------------
    // Check statistics
    // ----------------

    cout << "\nComparing Solver Statistics:\n";

    long int fe_evals = 0;
    if (type == MRISTEP_EXPLICIT || type == MRISTEP_IMEX)
    {
      fe_evals = mri_nst * nstages_stored + 1;
    }

    if (mri_nfse != fe_evals)
    {
      numfails++;
      cout << "Fe RHS evals: " << mri_nfse << " vs " << fe_evals << "\n";
    }

    long int fi_evals = 0;
    if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
    {
      fi_evals = mri_nst * nstages_stored + mri_nni;
      if (stiffly_accurate && !stiffly_accurate[i]) fi_evals++;
    }

    if (mri_nfsi != fi_evals)
    {
      numfails++;
      cout << "Fi RHS evals: " << mri_nfsi << " vs " << fi_evals << "\n";
    }

    if (numfails)
    {
      cout << "Failed " << numfails << " tests\n";
    }
    else
    {
      cout << "All checks passed\n";
    }

    // -------------------
    // Setup for next test
    // -------------------

    // Free coupling table
    MRIStepCoupling_Free(C);

    // Reset state vector to the initial condition
    N_VConst(RCONST(1.0), y);

    // Re-initialize fast integrator
    flag = ARKStepReInit(arkstep_mem, ff, nullptr, t0, y);
    if (check_flag(&flag, "ARKStepReInit", 1)) return 1;

    // Re-initialize slow integrator based on MRI type
    if (type == MRISTEP_EXPLICIT)
    {
      flag = MRIStepReInit(mristep_mem, fe, nullptr, t0, y);
    }
    else if (type == MRISTEP_IMPLICIT)
    {
      flag = MRIStepReInit(mristep_mem, nullptr, fi, t0, y);
    }
    else if (type == MRISTEP_IMEX)
    {
      flag = MRIStepReInit(mristep_mem, fe, fi, t0, y);
    }
    else
    {
      return 1;
    }
    if (check_flag(&flag, "MRIStepReInit", 1)) return 1;
  }

  // Clean up
  MRIStepInnerStepper_Free(&inner_stepper);
  MRIStepFree(&mristep_mem);
  ARKStepFree(&arkstep_mem);
  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX)
  {
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
  }
  N_VDestroy(y);
  delete[] methods;
  delete[] stiffly_accurate;

  return numfails;
}


// -----------------------------------------------------------------------------
// Functions called by the solver
// -----------------------------------------------------------------------------


// Explicit ODE RHS function fe(t,y)
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype* y_data  = N_VGetArrayPointer(y);
  realtype* yd_data = N_VGetArrayPointer(ydot);
  UserData* udata   = static_cast<UserData*>(user_data);

  yd_data[0] = udata->lambda_e * y_data[0];

  return 0;
}

// Implicit ODE RHS function fi(t,y)
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype* y_data  = N_VGetArrayPointer(y);
  realtype* yd_data = N_VGetArrayPointer(ydot);
  UserData* udata   = static_cast<UserData*>(user_data);

  yd_data[0] = udata->lambda_i * y_data[0];

  return 0;
}


// Fast ODE RHS function ff(t,y)
static int ff(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype* y_data  = N_VGetArrayPointer(y);
  realtype* yd_data = N_VGetArrayPointer(ydot);
  UserData* udata   = static_cast<UserData*>(user_data);

  yd_data[0] = udata->lambda_f * y_data[0];

  return 0;
}


// Jacobian routine to compute J(t,y) = dfi/dy.
static int Ji(realtype t, N_Vector y, N_Vector fy,
              SUNMatrix J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype* J_data = SUNDenseMatrix_Data(J);
  UserData* udata  = static_cast<UserData*>(user_data);

  J_data[0] = udata->lambda_i;

  return 0;
}


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------


// Check function return value
static int check_flag(void *flagvalue, const string funcname, int opt)
{
  int *errflag;

  // Check if function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == nullptr)
  {

    cerr << "\nMEMORY_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1;
  }
  // Check if flag < 0
  else if (opt == 1)
  {
    errflag = (int *) flagvalue;
    if (*errflag < 0)
    {
      cerr << "\nSUNDIALS_ERROR: " << funcname << " failed with flag = " << *errflag << "\n\n";
      return 1;
    }
  }

  return 0;
}
