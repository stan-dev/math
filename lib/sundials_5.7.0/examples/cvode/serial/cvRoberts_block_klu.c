/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * The following is a simple example problem based off of
 * cvRoberts_klu.c. We simulate a scenario where a set of independent
 * ODEs are grouped together to form a larger system. For simplicity,
 * each set of ODEs is the same problem. The problem is from chemical
 * kinetics, and consists of the following three rate equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration, the KLU sparse direct linear solver, and a user-supplied
 * Jacobian routine. It uses a scalar relative tolerance and a vector
 * absolute tolerance. Output is printed in decades from t = .4 to t =
 * 4.e10. Run statistics (optional outputs) are printed at the end.
 *
 * The program takes one optional argument, the number of groups
 * of independent ODE systems:
 *
 *    ./cvRoberts_block_klu [number of groups]
 *
 * The problem is comparable to the CUDA version -
 * cvRoberts_block_cusolversp_batchqr.cu. It was based off of the
 * cvRoberts_klu.c example.
 * -----------------------------------------------------------------*/

#include <stdio.h>

#include <cvode/cvode.h>                /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>     /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_klu.h>    /* access to KLU sparse direct solver   */
#include <sundials/sundials_types.h>    /* defs. of realtype, sunindextype      */

/* User-defined vector and matrix accessor macro: Ith */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..neq] and neq is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* Ith numbers components 1..neq */


/* Problem Constants */

#define GROUPSIZE 3            /* number of equations per group */
#define Y1    RCONST(1.0)      /* initial y components */
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1.0e-4)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-8)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-14)
#define ATOL3 RCONST(1.0e-6)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.4)      /* first output time      */
#define TMULT RCONST(10.0)     /* output time factor     */
#define NOUT  12               /* number of output times */

#define ZERO  RCONST(0.0)

/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

/* user data structure */
typedef struct {
  sunindextype ngroups;
  sunindextype neq;
} UserData;

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, char *argv[])
{
  realtype reltol, t, tout;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval, iout, nnz;
  sunindextype neq, ngroups, groupj;
  UserData udata;

  y = abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Parse command line arguments */
  if (argc > 1) {
    ngroups = atoi(argv[1]);
  } else {
    ngroups = 1000;
  }
  neq = ngroups * GROUPSIZE;

  udata.ngroups = ngroups;
  udata.neq = neq;

  /* Create serial vector of length neq for I.C. and abstol */
  y = N_VNew_Serial(neq);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
  abstol = N_VNew_Serial(neq);
  if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);

  /* Initialize y */
  for (groupj = 0; groupj < neq; groupj += GROUPSIZE) {
    Ith(y,1+groupj) = Y1;
    Ith(y,2+groupj) = Y2;
    Ith(y,3+groupj) = Y3;
  }

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  /* Set the vector absolute tolerance */
  for (groupj = 0; groupj < neq; groupj += GROUPSIZE) {
    Ith(abstol,1+groupj) = ATOL1;
    Ith(abstol,2+groupj) = ATOL2;
    Ith(abstol,3+groupj) = ATOL3;
  }

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSetUserData to attach the user data structure */
  retval = CVodeSetUserData(cvode_mem, &udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

  /* Create sparse SUNMatrix for use in linear solves */
  nnz = GROUPSIZE * GROUPSIZE * ngroups;
  A = SUNSparseMatrix(neq, neq, nnz, CSR_MAT);
  if(check_retval((void *)A, "SUNSparseMatrix", 0)) return(1);

  /* Create KLU solver object for use by CVode */
  LS = SUNLinSol_KLU(y, A);
  if(check_retval((void *)LS, "SUNLinSol_KLU", 0)) return(1);

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  printf(" \nGroup of independent 3-species kinetics problems\n\n");
  printf("number of groups = %lld\n\n", (long long int) ngroups);

  iout = 0;  tout = T1;
  while(1) {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    for (groupj = 0; groupj < 1; groupj++) {
      printf("group %lld: ", (long long int) groupj);
      PrintOutput(t, Ith(y,1+GROUPSIZE*groupj),
                     Ith(y,2+GROUPSIZE*groupj),
                     Ith(y,3+GROUPSIZE*groupj));
    }

    if (check_retval(&retval, "CVode", 1)) break;
    if (retval == CV_SUCCESS) {
      iout++;
      tout *= TMULT;
    }

    if (iout == NOUT) break;
  }

  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Free y and abstol vectors */
  N_VDestroy(y);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  SUNMatDestroy(A);

  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData *udata;
  sunindextype groupj;
  realtype y1, y2, y3, yd1, yd3;

  udata = (UserData*) user_data;

  for (groupj = 0; groupj < udata->neq; groupj += GROUPSIZE) {
    y1 = Ith(y,1+groupj); y2 = Ith(y,2+groupj); y3 = Ith(y,3+groupj);

    yd1 = Ith(ydot,1+groupj) = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
    yd3 = Ith(ydot,3+groupj) = RCONST(3.0e7)*y2*y2;
          Ith(ydot,2+groupj) = -yd1 - yd3;
  }

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *udata = (UserData*) user_data;
  sunindextype *rowptrs = SUNSparseMatrix_IndexPointers(J);
  sunindextype *colvals = SUNSparseMatrix_IndexValues(J);
  realtype *data = SUNSparseMatrix_Data(J);
  realtype *ydata;
  realtype y2, y3;
  sunindextype groupj, nnzper;

  ydata = N_VGetArrayPointer(y);
  nnzper = GROUPSIZE * GROUPSIZE;

  SUNMatZero(J);

  rowptrs[0] = 0;
  rowptrs = &rowptrs[1];
  for (groupj = 0; groupj < udata->ngroups; groupj++) {
    /* get y values */
    y2 = ydata[GROUPSIZE*groupj + 1];
    y3 = ydata[GROUPSIZE*groupj + 2];

    /* there are 3 entries per row */
    rowptrs[GROUPSIZE*groupj]     = 3 + nnzper*groupj;
    rowptrs[GROUPSIZE*groupj + 1] = 6 + nnzper*groupj;
    rowptrs[GROUPSIZE*groupj + 2] = 9 + nnzper*groupj;

    /* first row of block */
    data[nnzper*groupj]     = RCONST(-0.04);
    data[nnzper*groupj + 1] = RCONST(1.0e4)*y3;
    data[nnzper*groupj + 2] = RCONST(1.0e4)*y2;
    colvals[nnzper*groupj]     = GROUPSIZE*groupj;
    colvals[nnzper*groupj + 1] = GROUPSIZE*groupj + 1;
    colvals[nnzper*groupj + 2] = GROUPSIZE*groupj + 2;

    /* second row of block */
    data[nnzper*groupj + 3] = RCONST(0.04);
    data[nnzper*groupj + 4] = (RCONST(-1.0e4)*y3) - (RCONST(6.0e7)*y2);
    data[nnzper*groupj + 5] = RCONST(-1.0e4)*y2;
    colvals[nnzper*groupj + 3] = GROUPSIZE*groupj;
    colvals[nnzper*groupj + 4] = GROUPSIZE*groupj + 1;
    colvals[nnzper*groupj + 5] = GROUPSIZE*groupj + 2;

    /* third row of block */
    data[nnzper*groupj + 6] = ZERO;
    data[nnzper*groupj + 7] = RCONST(6.0e7)*y2;
    data[nnzper*groupj + 8] = ZERO;
    colvals[nnzper*groupj + 6] = GROUPSIZE*groupj;
    colvals[nnzper*groupj + 7] = GROUPSIZE*groupj + 1;
    colvals[nnzper*groupj + 8] = GROUPSIZE*groupj + 2;
  }

  return(0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}

/*
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nni, ncfn, netf, nge;
  int retval;

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

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld    nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
