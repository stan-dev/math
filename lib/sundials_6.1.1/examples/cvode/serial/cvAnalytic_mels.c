/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with analytical solution,
 * solution,
 *    dy/dt = lamda*y + 1/(1+t^2) - lamda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lamda".  The value of lamda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the BDF method, Newton
 * iteration, and a custom 'matrix-embedded' SUNLinearSolver. Output
 * is printed every 1.0 units of time (10 total).  Run statistics
 * (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied functions called by CVode */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Custom linear solver data structure, accessor macros, and routines */
static SUNLinearSolver MatrixEmbeddedLS(void *cvode_mem);
static SUNLinearSolver_Type MatrixEmbeddedLSType(SUNLinearSolver S);
static int MatrixEmbeddedLSSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                 N_Vector b, realtype tol);
static int MatrixEmbeddedLSFree(SUNLinearSolver S);

/* Private function to check computed solution */
static int check_ans(N_Vector y, realtype t, realtype rtol, realtype atol);

/* Private function to check function return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* SUNDIALS context */
static SUNContext sunctx = NULL;

/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);         /* initial time */
  realtype Tf = RCONST(10.0);        /* final time */
  realtype dTout = RCONST(1.0);      /* time between outputs */
  sunindextype NEQ = 1;              /* number of dependent vars. */
  realtype reltol = RCONST(1.0e-6);  /* tolerances */
  realtype abstol = RCONST(1.0e-10);
  realtype lamda  = RCONST(-100.0);  /* stiffness parameter */

  /* general problem variables */
  int retval;                     /* reusable error-checking flag */
  N_Vector y = NULL;              /* empty vector for storing solution */
  SUNLinearSolver LS = NULL;      /* empty linear solver object */
  void *cvode_mem = NULL;         /* empty CVode memory structure */
  realtype t, tout;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if(check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* Initial diagnostics output */
  printf("\nAnalytical ODE test problem:\n");
  printf("    lamda = %"GSYM"\n",    lamda);
  printf("   reltol = %.1"ESYM"\n",  reltol);
  printf("   abstol = %.1"ESYM"\n\n",abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ, sunctx);          /* Create serial vector for solution */
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return 1;
  N_VConst(RCONST(0.0), y);        /* Specify initial condition */

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSetUserData to specify the stiffness factor */
  retval = CVodeSetUserData(cvode_mem, (void *) &lamda);
  if (check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative and absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Create custom matrix-embedded linear solver */
  LS = MatrixEmbeddedLS(cvode_mem);
  if (check_retval((void *)LS, "MatrixEmbeddedLS", 0)) return 1;

  /* Call CVodeSetLinearSolver to attach the linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  t = T0;
  tout = T0+dTout;
  printf("        t           u\n");
  printf("   ---------------------\n");
  while (Tf - t > 1.0e-15) {

    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);             /* call integrator */
    if (check_retval(&retval, "CVode", 1)) break;
    printf("  %10.6"FSYM"  %10.6"FSYM"\n", t, NV_Ith_S(y,0));      /* access/print solution */
    if (retval >= 0) {                                             /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                       /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   ---------------------\n");

  /* Get/print some final statistics on how the solve progressed */
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);
  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li\n", nst);
  printf("   Total RHS evals = %li\n", nfe);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of linear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* check the solution error */
  retval = check_ans(y, t, reltol, abstol);

  /* Clean up and return */
  N_VDestroy(y);            /* Free y vector */
  CVodeFree(&cvode_mem);    /* Free integrator memory */
  SUNLinSolFree(LS);        /* Free linear solver */
  SUNContext_Free(&sunctx);

  return retval;
}


/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
  realtype lamda = rdata[0];                  /* set shortcut for stiffness parameter */
  realtype u = NV_Ith_S(y,0);                 /* access current solution value */

  /* fill in the RHS function: "NV_Ith_S" accesses the 0th entry of ydot */
  NV_Ith_S(ydot,0) = lamda*u + RCONST(1.0)/(RCONST(1.0)+t*t) - lamda*atan(t);

  return 0;                                   /* return with success */
}

/*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*/

/* constructor */
static SUNLinearSolver MatrixEmbeddedLS(void *cvode_mem)
{
  /* Create an empty linear solver */
  SUNLinearSolver LS = SUNLinSolNewEmpty(sunctx);
  if (LS == NULL) return NULL;

  /* Attach operations */
  LS->ops->gettype = MatrixEmbeddedLSType;
  LS->ops->solve   = MatrixEmbeddedLSSolve;
  LS->ops->free    = MatrixEmbeddedLSFree;

  /* Set content pointer to CVode memory */
  LS->content = cvode_mem;

  /* Return solver */
  return(LS);
}

/* type descriptor */
static SUNLinearSolver_Type MatrixEmbeddedLSType(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_MATRIX_EMBEDDED);
}

/* linear solve routine */
static int MatrixEmbeddedLSSolve(SUNLinearSolver LS, SUNMatrix A, N_Vector x,
                                 N_Vector b, realtype tol)
{
  /* temporary variables */
  int       retval;
  N_Vector  y, ypred, fn, zn1;
  realtype  tcur, gamma, rl1;
  void      *user_data;
  realtype  *rdata;
  realtype  lamda;

  /* retrieve implicit system data from CVode */
  retval = CVodeGetNonlinearSystemData(LS->content, &tcur, &ypred, &y, &fn,
                                       &gamma, &rl1, &zn1, &user_data);
  if (check_retval((void *)&retval, "CVodeGetNonlinearSystemData", 1))
    return(-1);

  /* extract stiffness parameter from user_data */
  rdata = (realtype *) user_data;
  lamda = rdata[0];

  /* perform linear solve: (1-gamma*lamda)*x = b */
  NV_Ith_S(x,0) = NV_Ith_S(b,0) / (1-gamma*lamda);

  /* return with success */
  return(SUNLS_SUCCESS);
}

/* destructor */
static int MatrixEmbeddedLSFree(SUNLinearSolver LS)
{
  if (LS == NULL) return(SUNLS_SUCCESS);
  LS->content = NULL;
  SUNLinSolFreeEmpty(LS);
  return(SUNLS_SUCCESS);
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
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *retval);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}

/* check the computed solution */
static int check_ans(N_Vector y, realtype t, realtype rtol, realtype atol)
{
  int      passfail=0;     /* answer pass (0) or fail (1) flag     */
  realtype ans, err, ewt;  /* answer data, error, and error weight */

  /* compute solution error */
  ans = atan(t);
  ewt = RCONST(1.0) / (rtol * fabs(ans) + atol);
  err = ewt * fabs(NV_Ith_S(y,0) - ans);

  /* is the solution within the tolerances? */
  passfail = (err < RCONST(1.0)) ? 0 : 1;

  if (passfail) {
    fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%"GSYM"\n\n", err);
  }

  return(passfail);
}
