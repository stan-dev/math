/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is an ODE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions
 * Y0 = [u0,v0,w0].
 *
 * We have 3 different testing scenarios:
 *
 * Test 1:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
 *    Here, all three components exhibit a rapid transient change
 *    during the first 0.2 time units, followed by a slow and
 *    smooth evolution.
 *
 * Test 2:  u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6
 *    Here, w experiences a fast initial transient, jumping 0.5
 *    within a few steps.  All values proceed smoothly until
 *    around t=6.5, when both u and v undergo a sharp transition,
 *    with u increaseing from around 0.5 to 5 and v decreasing
 *    from around 6 to 1 in less than 0.5 time units.  After this
 *    transition, both u and v continue to evolve somewhat
 *    rapidly for another 1.4 time units, and finish off smoothly.
 *
 * Test 3:  u0=3,  v0=3,  w0=3.5,  a=0.5,  b=3,  ep=5.0e-4
 *    Here, all components undergo very rapid initial transients
 *    during the first 0.3 time units, and all then proceed very
 *    smoothly for the remainder of the simulation.
 *
 * This file is hard-coded to use test 2.
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the SUNDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>      /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_serial.h>     /* serial N_Vector types, fcts., macros */
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>    /* def. of type 'realtype' */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);     /* initial time */
  realtype Tf = RCONST(10.0);    /* final time */
  realtype dTout = RCONST(1.0);  /* time between outputs */
  sunindextype NEQ = 3;          /* number of dependent vars. */
  int Nt = ceil(Tf/dTout);       /* number of output times */
  int test = 2;                  /* test problem to run */
  realtype reltol = 1.0e-6;      /* tolerances */
  realtype abstol = 1.0e-10;
  realtype a, b, ep, u0, v0, w0;

  /* general problem variables */
  int flag;                      /* reusable error-checking flag */
  N_Vector y = NULL;             /* empty vector for storing solution */
  SUNMatrix A = NULL;            /* empty matrix for solver */
  SUNLinearSolver LS = NULL;     /* empty linear solver object */
  void *arkode_mem = NULL;       /* empty ARKode memory structure */
  realtype rdata[3];
  FILE *UFID;
  realtype t, tout;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;

  /* set up the test problem according to the desired test */
  if (test == 1) {
    u0 = RCONST(3.9);
    v0 = RCONST(1.1);
    w0 = RCONST(2.8);
    a  = RCONST(1.2);
    b  = RCONST(2.5);
    ep = RCONST(1.0e-5);
  } else if (test == 3) {
    u0 = RCONST(3.0);
    v0 = RCONST(3.0);
    w0 = RCONST(3.5);
    a  = RCONST(0.5);
    b  = RCONST(3.0);
    ep = RCONST(5.0e-4);
  } else {
    u0 = RCONST(1.2);
    v0 = RCONST(3.1);
    w0 = RCONST(3.0);
    a  = RCONST(1.0);
    b  = RCONST(3.5);
    ep = RCONST(5.0e-6);
  }

  /* Initial problem output */
  printf("\nBrusselator ODE test problem:\n");
  printf("    initial conditions:  u0 = %"GSYM",  v0 = %"GSYM",  w0 = %"GSYM"\n",u0,v0,w0);
  printf("    problem parameters:  a = %"GSYM",  b = %"GSYM",  ep = %"GSYM"\n",a,b,ep);
  printf("    reltol = %.1"ESYM",  abstol = %.1"ESYM"\n\n",reltol,abstol);

  /* Initialize data structures */
  rdata[0] = a;     /* set user data  */
  rdata[1] = b;
  rdata[2] = ep;
  y = N_VNew_Serial(NEQ);           /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  NV_Ith_S(y,0) = u0;               /* Set initial conditions */
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y);
  if (check_flag((void *)arkode_mem, "ARKStepCreate", 0)) return 1;

  /* Set routines */
  flag = ARKStepSetUserData(arkode_mem, (void *) rdata);     /* Pass rdata to user functions */
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
  flag = ARKStepSStolerances(arkode_mem, reltol, abstol);    /* Specify tolerances */
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  /* Initialize dense matrix data structure and solver */
  A = SUNDenseMatrix(NEQ, NEQ);
  if (check_flag((void *)A, "SUNDenseMatrix", 0)) return 1;
  LS = SUNLinSol_Dense(y, A);
  if (check_flag((void *)LS, "SUNLinSol_Dense", 0)) return 1;

  /* Linear solver interface */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, A);        /* Attach matrix and linear solver */
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
  flag = ARKStepSetJacFn(arkode_mem, Jac);                 /* Set Jacobian routine */
  if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt","w");
  fprintf(UFID,"# t u v w\n");

  /* output initial condition to disk */
  fprintf(UFID," %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM"\n",
          T0, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));

  /* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  tout = T0+dTout;
  printf("        t           u           v           w\n");
  printf("   -------------------------------------------\n");
  printf("  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"\n",
         t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));

  for (iout=0; iout<Nt; iout++) {

    flag = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);      /* call integrator */
    if (check_flag(&flag, "ARKStepEvolve", 1)) break;
    printf("  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"\n",             /* access/print solution */
           t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
    fprintf(UFID," %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM"\n",
            t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
    if (flag >= 0) {                                         /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                 /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   -------------------------------------------\n");
  fclose(UFID);

  /* Print some final statistics */
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKStepGetNumSteps", 1);
  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKStepGetNumStepAttempts", 1);
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKStepGetNumRhsEvals", 1);
  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1);
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKStepGetNumErrTestFails", 1);
  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1);
  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1);
  flag = ARKStepGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKStepGetNumJacEvals", 1);
  flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of linear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy(y);               /* Free y vector */
  ARKStepFree(&arkode_mem);    /* Free integrator memory */
  SUNLinSolFree(LS);           /* Free linear solver */
  SUNMatDestroy(A);            /* Free A matrix */
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
  realtype a  = rdata[0];                     /* access data entries */
  realtype b  = rdata[1];
  realtype ep = rdata[2];
  realtype u = NV_Ith_S(y,0);                 /* access solution values */
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* fill in the RHS function */
  NV_Ith_S(ydot,0) = a - (w+1.0)*u + v*u*u;
  NV_Ith_S(ydot,1) = w*u - v*u*u;
  NV_Ith_S(ydot,2) = (b-w)/ep - w*u;

  return 0;                                  /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
  realtype ep = rdata[2];                     /* access data entries */
  realtype u = NV_Ith_S(y,0);                 /* access solution values */
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* fill in the Jacobian via SUNDenseMatrix macro, SM_ELEMENT_D (see sunmatrix_dense.h) */
  SM_ELEMENT_D(J,0,0) = -(w+1.0) + 2.0*u*v;
  SM_ELEMENT_D(J,0,1) = u*u;
  SM_ELEMENT_D(J,0,2) = -u;

  SM_ELEMENT_D(J,1,0) = w - 2.0*u*v;
  SM_ELEMENT_D(J,1,1) = -u*u;
  SM_ELEMENT_D(J,1,2) = u;

  SM_ELEMENT_D(J,2,0) = -w;
  SM_ELEMENT_D(J,2,1) = 0.0;
  SM_ELEMENT_D(J,2,2) = -1.0/ep - u;

  return 0;                                   /* Return with success */
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
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


/*---- end of file ----*/
