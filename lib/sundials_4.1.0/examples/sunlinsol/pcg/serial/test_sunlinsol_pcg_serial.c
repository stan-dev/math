/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol PCG module 
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_pcg.h>
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
  N_Vector s;     /* scaling vector supplied to PCG */
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
 * SUNOCG Linear Solver Testing Routine
 *
 * We run multiple tests to exercise this solver:
 * 1. simple tridiagonal system (no preconditioning)
 * 2. simple tridiagonal system (Jacobi preconditioning)
 * 3. tridiagonal system w/ scale vector s (no preconditioning)
 * 4. tridiagonal system w/ scale vector s (Jacobi preconditioning)
 *
 * Note: We construct a tridiagonal matrix Ahat, a random solution 
 *       xhat, and a corresponding rhs vector bhat = Ahat*xhat, such 
 *       that each of these is unit-less.  To test scaling, we use 
 *       the matrix 
 *             A = (S-inverse) Ahat (S-inverse), 
 *       solution vector 
 *             x = S xhat; 
 *       and construct b = A*x.  Hence the linear system has both rows 
 *       and columns scaled by (S-inverse), where S is the diagonal 
 *       matrix with entries from the vector s, the 'scaling' vector 
 *       supplied to PCG having strictly positive entries.  
 *
 *       When this is combined with preconditioning, we construct 
 *       P \approx (A-inverse) by taking a unit-less preconditioner 
 *       Phat \approx (Ahat-inverse), and constructing the operator
 *       P via
 *             P = S Phat S \approx S (Ahat-inverse) S = A-inverse
 *       We apply this via the steps:
 *             z = Pr = S Phat S r
 *       Since both S and Phat are diagonal matrices, this is 
 *       equivalent to
 *             z(i) = s(i)^2 Phat(i) r(i)
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int             fails=0;          /* counter for test failures */
  int             passfail=0;       /* overall pass/fail flag    */
  SUNLinearSolver LS;               /* linear solver object      */
  N_Vector        xhat, x, b;       /* test vectors              */
  UserData        ProbData;         /* problem data structure    */
  int             maxl, print_timing;
  sunindextype    i;
  realtype        *vecdata;
  double          tol;

  /* check inputs: local problem size, timing flag */
  if (argc < 5) {
    printf("ERROR: FOUR (4) Inputs required:\n");
    printf("  Problem size should be >0\n");
    printf("  Maximum Krylov subspace dimension should be >0\n");
    printf("  Solver tolerance should be >0\n");
    printf("  timing output flag should be 0 or 1 \n");
    return 1;
  }
  ProbData.N = atol(argv[1]);
  problem_size = ProbData.N;
  if (ProbData.N <= 0) {
    printf("ERROR: Problem size must be a positive integer\n");
    return 1; 
  }
  maxl = atoi(argv[2]);
  if (maxl <= 0) {
    printf("ERROR: Maximum Krylov subspace dimension must be a positive integer\n");
    return 1; 
  }
  tol = atof(argv[3]);
  if (tol <= ZERO) {
    printf("ERROR: Solver tolerance must be a positive real number\n");
    return 1; 
  }
  print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  printf("\nPCG linear solver test:\n");
  printf("  Problem size = %ld\n", (long int) ProbData.N);
  printf("  Maximum Krylov subspace dimension = %i\n", maxl);
  printf("  Solver Tolerance = %"GSYM"\n", tol);
  printf("  timing output flag = %i\n\n", print_timing);
  
  /* Create vectors */
  x = N_VNew_Serial(ProbData.N);
  if (check_flag(x, "N_VNew_Serial", 0)) return 1;
  xhat = N_VNew_Serial(ProbData.N);
  if (check_flag(xhat, "N_VNew_Serial", 0)) return 1;
  b = N_VNew_Serial(ProbData.N);
  if (check_flag(b, "N_VNew_Serial", 0)) return 1;
  ProbData.d = N_VNew_Serial(ProbData.N);
  if (check_flag(ProbData.d, "N_VNew_Serial", 0)) return 1;
  ProbData.s = N_VNew_Serial(ProbData.N);
  if (check_flag(ProbData.s, "N_VNew_Serial", 0)) return 1;

  /* Fill xhat vector with uniform random data in [1,2] */
  vecdata = N_VGetArrayPointer(xhat);
  for (i=0; i<ProbData.N; i++) 
    vecdata[i] = ONE + urand();

  /* Fill Jacobi vector with matrix diagonal */
  N_VConst(FIVE, ProbData.d);
  
  /* Create PCG linear solver */
  LS = SUNLinSol_PCG(x, PREC_RIGHT, maxl);
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_ITERATIVE, 0);
  fails += Test_SUNLinSolSetATimes(LS, &ProbData, ATimes, 0);
  fails += Test_SUNLinSolSetPreconditioner(LS, &ProbData, PSetup, PSolve, 0);
  fails += Test_SUNLinSolSetScalingVectors(LS, ProbData.s, NULL, 0);
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);
  if (fails) {
    printf("FAIL: SUNLinSol_PCG module failed %i initialization tests\n\n", fails);
    return 1;
  } else {
    printf("SUCCESS: SUNLinSol_PCG module passed all initialization tests\n\n");
  }


  
  /*** Test 1: simple Poisson-like solve (no preconditioning) ***/

  /* set scaling vector */
  N_VConst(ONE, ProbData.s);

  /* Fill x vector with scaled version */
  N_VProd(xhat, ProbData.s, x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run test with this setup */
  fails += SUNLinSol_PCGSetPrecType(LS, PREC_NONE);  
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);
  
  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_PCG module, problem 1, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_PCG module, problem 1, passed all tests\n\n");
  }

  
  /*** Test 2: simple Poisson-like solve (Jacobi preconditioning) ***/

  /* set scaling vector */
  N_VConst(ONE, ProbData.s);

  /* Fill x vector with scaled version */
  N_VProd(xhat, ProbData.s, x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_PCGSetPrecType(LS, PREC_RIGHT);  
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_PCG module, problem 2, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_PCG module, problem 2, passed all tests\n\n");
  }

  
  /*** Test 3: Poisson-like solve w/ scaling (no preconditioning) ***/

  /* set scaling vector */
  vecdata = N_VGetArrayPointer(ProbData.s);
  for (i=0; i<ProbData.N; i++)
    vecdata[i] = ONE + THOUSAND*urand();

  /* Fill x vector with scaled version */
  N_VProd(xhat, ProbData.s, x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_PCGSetPrecType(LS, PREC_NONE);  
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_PCG module, problem 3, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_PCG module, problem 3, passed all tests\n\n");
  }

  
  /*** Test 4: Poisson-like solve w/ scaling (Jacobi preconditioning) ***/

  /* set scaling vectors */
  vecdata = N_VGetArrayPointer(ProbData.s);
  for (i=0; i<ProbData.N; i++)
    vecdata[i] = ONE + THOUSAND*urand();

  /* Fill x vector with scaled version */
  N_VProd(xhat, ProbData.s, x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_PCGSetPrecType(LS, PREC_RIGHT);  
  fails += Test_SUNLinSolSetup(LS, NULL, 0);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolResid(LS, 0);

  /* Print result */
  if (fails) { 
    printf("FAIL: SUNLinSol_PCG module, problem 4, failed %i tests\n\n", fails);
    passfail += 1;
  } else {
    printf("SUCCESS: SUNLinSol_PCG module, problem 4, passed all tests\n\n");
  }

  
  /* Free solver and vectors */
  SUNLinSolFree(LS);
  N_VDestroy(x);
  N_VDestroy(xhat);
  N_VDestroy(b);
  N_VDestroy(ProbData.d);
  N_VDestroy(ProbData.s);

  return(passfail);
}


/* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*/

/* matrix-vector product  */
int ATimes(void* Data, N_Vector v_vec, N_Vector z_vec)
{
  /* local variables */
  realtype *v, *z, *s;
  sunindextype i, N;
  UserData *ProbData;
  
  /* access user data structure and vector data */
  ProbData = (UserData *) Data;
  v = N_VGetArrayPointer(v_vec);
  if (check_flag(v, "N_VGetArrayPointer", 0)) return 1;
  z = N_VGetArrayPointer(z_vec);
  if (check_flag(z, "N_VGetArrayPointer", 0)) return 1;
  s = N_VGetArrayPointer(ProbData->s);
  if (check_flag(s, "N_VGetArrayPointer", 0)) return 1;
  N = ProbData->N;  

  /* perform product at left boundary (note: v is zero at the boundary)*/
  z[0] = (FIVE*v[0]/s[0] - v[1]/s[1])/s[0];

  /* iterate through interior of domain, performing product */
  for (i=1; i<N-1; i++) 
    z[i] = (-v[i-1]/s[i-1] + FIVE*v[i]/s[i] - v[i+1]/s[i+1])/s[i];

  /* perform product at right boundary (note: v is zero at the boundary)*/  
  z[N-1] = (-v[N-2]/s[N-2] + FIVE*v[N-1]/s[N-1])/s[N-1];
  
  /* return with success */
  return 0;
}
  
/* preconditioner setup -- nothing to do here since everything is already stored */
int PSetup(void* Data) { return 0; }
  
/* preconditioner solve */
int PSolve(void* Data, N_Vector r_vec, N_Vector z_vec, realtype tol, int lr)
{
  /* local variables */
  realtype *r, *z, *d, *s;
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
  s = N_VGetArrayPointer(ProbData->s);
  if (check_flag(s, "N_VGetArrayPointer", 0)) return 1;
  
  /* iterate through domain, performing Jacobi solve */
  for (i=0; i<ProbData->N; i++) 
    z[i] = s[i] * s[i] * r[i] / d[i];

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
  int      failure = 0;
  long int i;
  realtype *Xdata, *Ydata, maxerr;
  
  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);
  
  /* check vector data */
  for(i=0; i<problem_size; i++)
    failure += FNEQ(Xdata[i], Ydata[i], FIVE*tol*SUNRabs(Xdata[i]));

  if (failure > ZERO) {
    maxerr = ZERO;
    for(i=0; i < problem_size; i++)
      maxerr = SUNMAX(SUNRabs(Xdata[i]-Ydata[i])/SUNRabs(Xdata[i]), maxerr);
    printf("check err failure: maxerr = %g (tol = %g)\n",
	   maxerr, FIVE*tol);
    return(1);
  }
  else
    return(0);
}
