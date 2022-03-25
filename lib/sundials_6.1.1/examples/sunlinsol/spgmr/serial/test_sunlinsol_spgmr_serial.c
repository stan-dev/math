/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
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
 * This is the testing routine to check the SUNLinSol SPGMR module
 * implementation.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* constants */
#define FIVE      RCONST(5.0)
#define THOUSAND  RCONST(1000.0)

/* user data structure */
typedef struct {
  sunindextype N; /* problem size */
  N_Vector d;     /* matrix diagonal */
  N_Vector s1;    /* scaling vectors supplied to SPGMR */
  N_Vector s2;
} UserData;

/* private functions */
/*    matrix-vector product  */
int ATimes(void* ProbData, N_Vector v, N_Vector z);
/*    preconditioner setup */
int PSetup(void* ProbData);
/*    preconditioner solve */
int PSolve(void* ProbData, N_Vector r, N_Vector z, realtype tol, int lr);
/*    checks function return values  */
static int check_flag(void *flagvalue, const char *funcname, int opt);
/*    uniform random number generator in [0,1] */
static realtype urand();

/* global copy of the problem size (for check_vector routine) */
sunindextype problem_size;

/* ----------------------------------------------------------------------
 * SUNLinSol_SPGMR Linear Solver Testing Routine
 *
 * We run multiple tests to exercise this solver:
 * 1. simple tridiagonal system (no preconditioning)
 * 2. simple tridiagonal system (Jacobi preconditioning)
 * 3. tridiagonal system w/ scale vector s1 (no preconditioning)
 * 4. tridiagonal system w/ scale vector s1 (Jacobi preconditioning)
 * 5. tridiagonal system w/ scale vector s2 (no preconditioning)
 * 6. tridiagonal system w/ scale vector s2 (Jacobi preconditioning)
 *
 * Note: We construct a tridiagonal matrix Ahat, a random solution xhat,
 *       and a corresponding rhs vector bhat = Ahat*xhat, such that each
 *       of these is unit-less.  To test row/column scaling, we use the
 *       matrix A = S1-inverse Ahat S2, rhs vector b = S1-inverse bhat,
 *       and solution vector x = (S2-inverse) xhat; hence the linear
 *       system has rows scaled by S1-inverse and columns scaled by S2,
 *       where S1 and S2 are the diagonal matrices with entries from the
 *       vectors s1 and s2, the 'scaling' vectors supplied to SPGMR
 *       having strictly positive entries.  When this is combined with
 *       preconditioning, assume that Phat is the desired preconditioner
 *       for Ahat, then our preconditioning matrix P \approx A should be
 *         left prec:  P-inverse \approx S1-inverse Ahat-inverse S1
 *         right prec:  P-inverse \approx S2-inverse Ahat-inverse S2.
 *       Here we use a diagonal preconditioner D, so the S*-inverse
 *       and S* in the product cancel one another.
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int             fails=0;          /* counter for test failures */
  int             passfail=0;       /* overall pass/fail flag    */
  SUNLinearSolver LS;               /* linear solver object      */
  N_Vector        xhat, x, b;       /* test vectors              */
  UserData        ProbData;         /* problem data structure    */
  int             gstype, pretype, maxl, print_timing;
  sunindextype    i;
  realtype        *vecdata;
  double          tol;
  SUNContext      sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  /* check inputs: local problem size, timing flag */
  if (argc < 7) {
    printf("ERROR: SIX (6) Inputs required:\n");
    printf("  Problem size should be >0\n");
    printf("  Gram-Schmidt orthogonalization type should be 1 or 2\n");
    printf("  Preconditioning type should be 1 or 2\n");
    printf("  Maximum Krylov subspace dimension should be >0\n");
    printf("  Solver tolerance should be >0\n");
    printf("  timing output flag should be 0 or 1 \n");
    return 1;
  }
  ProbData.N = (sunindextype) atol(argv[1]);
  problem_size = ProbData.N;
  if (ProbData.N <= 0) {
    printf("ERROR: Problem size must be a positive integer\n");
    return 1;
  }
  gstype = atoi(argv[2]);
  if ((gstype < 1) || (gstype > 2)) {
    printf("ERROR: Gram-Schmidt orthogonalization type must be either 1 or 2\n");
    return 1;
  }
  pretype = atoi(argv[3]);
  if ((pretype < 1) || (pretype > 2)) {
    printf("ERROR: Preconditioning type must be either 1 or 2\n");
    return 1;
  }
  maxl = atoi(argv[4]);
  if (maxl <= 0) {
    printf("ERROR: Maximum Krylov subspace dimension must be a positive integer\n");
    return 1;
  }
  tol = atof(argv[5]);
  if (tol <= ZERO) {
    printf("ERROR: Solver tolerance must be a positive real number\n");
    return 1;
  }
  print_timing = atoi(argv[6]);
  SetTiming(print_timing);

  printf("\nSPGMR linear solver test:\n");
  printf("  Problem size = %ld\n", (long int) ProbData.N);
  printf("  Gram-Schmidt orthogonalization type = %i\n", gstype);
  printf("  Preconditioning type = %i\n", pretype);
  printf("  Maximum Krylov subspace dimension = %i\n", maxl);
  printf("  Solver Tolerance = %g\n", tol);
  printf("  timing output flag = %i\n\n", print_timing);

  /* Create vectors */
  x = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(x, "N_VNew_Serial", 0)) return 1;
  xhat = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(xhat, "N_VNew_Serial", 0)) return 1;
  b = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(b, "N_VNew_Serial", 0)) return 1;
  ProbData.d = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(ProbData.d, "N_VNew_Serial", 0)) return 1;
  ProbData.s1 = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(ProbData.s1, "N_VNew_Serial", 0)) return 1;
  ProbData.s2 = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(ProbData.s2, "N_VNew_Serial", 0)) return 1;

  /* Fill xhat vector with uniform random data in [1,2] */
  vecdata = N_VGetArrayPointer(xhat);
  for (i=0; i<ProbData.N; i++)
    vecdata[i] = ONE + urand();

  /* Fill Jacobi vector with matrix diagonal */
  N_VConst(FIVE, ProbData.d);

  /* Create SPGMR linear solver */
  LS = SUNLinSol_SPGMR(x, pretype, maxl, sunctx);
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_ITERATIVE, 0);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_SPGMR, 0);
  fails += Test_SUNLinSolSetATimes(LS, &ProbData, ATimes, 0);
  fails += Test_SUNLinSolSetPreconditioner(LS, &ProbData, PSetup, PSolve, 0);
  fails += Test_SUNLinSolSetScalingVectors(LS, ProbData.s1, ProbData.s2, 0);
  fails += Test_SUNLinSolSetZeroGuess(LS, 0);
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);
  fails += SUNLinSol_SPGMRSetGSType(LS, gstype);
  if (fails) {
    printf("FAIL: SUNLinSol_SPGMR module failed %i initialization tests\n\n", fails);
    return 1;
  } else {
    printf("SUCCESS: SUNLinSol_SPGMR module passed all initialization tests\n\n");
  }


  /*** Test 1: simple Poisson-like solve (no preconditioning) ***/

  /* set scaling vectors */
  N_VConst(ONE, ProbData.s1);
  N_VConst(ONE, ProbData.s2);

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPGMRSetPrecType(LS, SUN_PREC_NONE);
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPGMR module, problem 1, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_SPGMR module, problem 1, passed all tests\n\n");
  }


  /*** Test 2: simple Poisson-like solve (Jacobi preconditioning) ***/

  /* set scaling vectors */
  N_VConst(ONE,  ProbData.s1);
  N_VConst(ONE,  ProbData.s2);

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPGMRSetPrecType(LS, pretype);
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPGMR module, problem 2, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_SPGMR module, problem 2, passed all tests\n\n");
  }


  /*** Test 3: Poisson-like solve w/ scaled rows (no preconditioning) ***/

  /* set scaling vectors */
  vecdata = N_VGetArrayPointer(ProbData.s1);
  for (i=0; i<ProbData.N; i++)
    vecdata[i] = ONE + THOUSAND*urand();
  N_VConst(ONE, ProbData.s2);

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPGMRSetPrecType(LS, SUN_PREC_NONE);
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPGMR module, problem 3, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_SPGMR module, problem 3, passed all tests\n\n");
  }


  /*** Test 4: Poisson-like solve w/ scaled rows (Jacobi preconditioning) ***/

  /* set scaling vectors */
  vecdata = N_VGetArrayPointer(ProbData.s1);
  for (i=0; i<ProbData.N; i++)
    vecdata[i] = ONE + THOUSAND*urand();
  N_VConst(ONE, ProbData.s2);

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPGMRSetPrecType(LS, pretype);
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPGMR module, problem 4, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_SPGMR module, problem 4, passed all tests\n\n");
  }


  /*** Test 5: Poisson-like solve w/ scaled columns (no preconditioning) ***/

  /* set scaling vectors */
  N_VConst(ONE, ProbData.s1);
  vecdata = N_VGetArrayPointer(ProbData.s2);
  for (i=0; i<ProbData.N; i++)
    vecdata[i] = ONE + THOUSAND*urand();

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPGMRSetPrecType(LS, SUN_PREC_NONE);
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPGMR module, problem 5, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_SPGMR module, problem 5, passed all tests\n\n");
  }


  /*** Test 6: Poisson-like solve w/ scaled columns (Jacobi preconditioning) ***/

  /* set scaling vector, Jacobi solver vector */
  N_VConst(ONE, ProbData.s1);
  vecdata = N_VGetArrayPointer(ProbData.s2);
  for (i=0; i<ProbData.N; i++)
    vecdata[i] = ONE + THOUSAND*urand();

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPGMRSetPrecType(LS, pretype);
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPGMR module, problem 6, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_SPGMR module, problem 6, passed all tests\n\n");
  }


  /* Free solver and vectors */
  SUNLinSolFree(LS);
  N_VDestroy(x);
  N_VDestroy(xhat);
  N_VDestroy(b);
  N_VDestroy(ProbData.d);
  N_VDestroy(ProbData.s1);
  N_VDestroy(ProbData.s2);
  SUNContext_Free(&sunctx);

  return(passfail);
}


/* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*/

/* matrix-vector product  */
int ATimes(void* Data, N_Vector v_vec, N_Vector z_vec)
{
  /* local variables */
  realtype *v, *z, *s1, *s2;
  sunindextype i, N;
  UserData *ProbData;

  /* access user data structure and vector data */
  ProbData = (UserData *) Data;
  v = N_VGetArrayPointer(v_vec);
  if (check_flag(v, "N_VGetArrayPointer", 0)) return 1;
  z = N_VGetArrayPointer(z_vec);
  if (check_flag(z, "N_VGetArrayPointer", 0)) return 1;
  s1 = N_VGetArrayPointer(ProbData->s1);
  if (check_flag(s1, "N_VGetArrayPointer", 0)) return 1;
  s2 = N_VGetArrayPointer(ProbData->s2);
  if (check_flag(s2, "N_VGetArrayPointer", 0)) return 1;
  N = ProbData->N;

  /* perform product at the left domain boundary (note: v is zero at the boundary)*/
  z[0] = (FIVE*v[0]*s2[0] - v[1]*s2[1])/s1[0];

  /* iterate through interior of local domain, performing product */
  for (i=1; i<N-1; i++)
    z[i] = (-v[i-1]*s2[i-1] + FIVE*v[i]*s2[i] - v[i+1]*s2[i+1])/s1[i];

  /* perform product at the right domain boundary (note: v is zero at the boundary)*/
  z[N-1] = (-v[N-2]*s2[N-2] + FIVE*v[N-1]*s2[N-1])/s1[N-1];

  /* return with success */
  return 0;
}

/* preconditioner setup -- nothing to do here since everything is already stored */
int PSetup(void* Data) { return 0; }

/* preconditioner solve */
int PSolve(void* Data, N_Vector r_vec, N_Vector z_vec, realtype tol, int lr)
{
  /* local variables */
  realtype *r, *z, *d;
  sunindextype i;
  UserData *ProbData;

  /* access user data structure and vector data */
  ProbData = (UserData *) Data;
  r = N_VGetArrayPointer(r_vec);
  if (check_flag(r, "N_VGetArrayPointer", 0)) return 1;
  z = N_VGetArrayPointer(z_vec);
  if (check_flag(z, "N_VGetArrayPointer", 0)) return 1;
  d = N_VGetArrayPointer(ProbData->d);
  if (check_flag(d, "N_VGetArrayPointer", 0)) return 1;

  /* iterate through domain, performing Jacobi solve */
  for (i=0; i<ProbData->N; i++)
    z[i] = r[i] / d[i];

  /* return with success */
  return 0;
}

/* uniform random number generator */
static realtype urand()
{
  return ((realtype) rand() / (realtype) RAND_MAX);
}

/* Check function return value based on "opt" input:
     0:  function allocates memory so check for NULL pointer
     1:  function returns a flag so check for flag != 0 */
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if function returned NULL pointer - no memory allocated */
  if (opt==0 && flagvalue==NULL) {
    fprintf(stderr, "\nERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return 1; }

  /* Check if flag != 0 */
  if (opt==1) {
    errflag = (int *) flagvalue;
    if (*errflag != 0) {
      fprintf(stderr, "\nERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return 1; }}

  return 0;
}


/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector X, N_Vector Y, realtype tol)
{
  int failure = 0;
  sunindextype i;
  realtype *Xdata, *Ydata, maxerr;

  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);

  /* check vector data */
  for(i=0; i<problem_size; i++)
    failure += SUNRCompareTol(Xdata[i], Ydata[i], tol);

  if (failure > ZERO) {
    maxerr = ZERO;
    for(i=0; i < problem_size; i++)
      maxerr = SUNMAX(SUNRabs(Xdata[i]-Ydata[i])/SUNRabs(Xdata[i]), maxerr);
    printf("check err failure: maxerr = %"GSYM" (tol = %"GSYM")\n",
	   maxerr, tol);
    return(1);
  }
  else
    return(0);
}

void sync_device()
{
}
