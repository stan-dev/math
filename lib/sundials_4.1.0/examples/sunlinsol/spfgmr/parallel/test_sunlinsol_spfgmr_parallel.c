/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
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
 * This is the testing routine to check the SUNLinSol SPFGMR module 
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"
#include "mpi.h"

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
  sunindextype Nloc;  /* local problem size */
  N_Vector d;         /* matrix diagonal */
  N_Vector s1;        /* scaling vectors supplied to SPFGMR */
  N_Vector s2;
  MPI_Comm comm;      /* communicator object */
  int myid;           /* MPI process ID */
  int nprocs;         /* total number of MPI processes */
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

/* global copy of Nloc (for check_vector routine) */
sunindextype local_problem_size;

/* ----------------------------------------------------------------------
 * SUNLinSol_SPFGMR Linear Solver Testing Routine
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
 *       vectors s1 and s2, the 'scaling' vectors supplied to SPFGMR 
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
  int             gstype, maxl, print_timing;
  sunindextype    i;
  realtype        *vecdata;
  double          tol;

  /* Set up MPI environment */
  fails = MPI_Init(&argc, &argv);
  if (check_flag(&fails, "MPI_Init", 1)) return 1;
  ProbData.comm = MPI_COMM_WORLD;
  fails = MPI_Comm_size(ProbData.comm, &(ProbData.nprocs));
  if (check_flag(&fails, "MPI_Comm_size", 1)) return 1;
  fails = MPI_Comm_rank(ProbData.comm, &(ProbData.myid));
  if (check_flag(&fails, "MPI_Comm_rank", 1)) return 1;

  /* check inputs: local problem size, timing flag */
  if (argc < 6) {
    printf("ERROR: FIVE (5) Inputs required:\n");
    printf("  Local problem size should be >0\n");
    printf("  Gram-Schmidt orthogonalization type should be 1 or 2\n");
    printf("  Maximum Krylov subspace dimension should be >0\n");
    printf("  Solver tolerance should be >0\n");
    printf("  timing output flag should be 0 or 1 \n");
    return 1;
  }
  ProbData.Nloc = atol(argv[1]);
  local_problem_size = ProbData.Nloc;
  if (ProbData.Nloc <= 0) {
    printf("ERROR: local problem size must be a positive integer\n");
    return 1; 
  }
  gstype = atoi(argv[2]);
  if ((gstype < 1) || (gstype > 2)) {
    printf("ERROR: Gram-Schmidt orthogonalization type must be either 1 or 2\n");
    return 1; 
  }
  maxl = atoi(argv[3]);
  if (maxl <= 0) {
    printf("ERROR: Maximum Krylov subspace dimension must be a positive integer\n");
    return 1; 
  }
  tol = atof(argv[4]);
  if (tol <= ZERO) {
    printf("ERROR: Solver tolerance must be a positive real number\n");
    return 1; 
  }
  print_timing = atoi(argv[5]);
  SetTiming(print_timing);

  if (ProbData.myid == 0) {
    printf("\nSPFGMR linear solver test:\n");
    printf("  nprocs = %i\n", ProbData.nprocs);
    printf("  local/global problem sizes = %ld/%ld\n", (long int) ProbData.Nloc,
           (long int) (ProbData.nprocs * ProbData.Nloc));
    printf("  Gram-Schmidt orthogonalization type = %i\n", gstype);
    printf("  Maximum Krylov subspace dimension = %i\n", maxl);
    printf("  Solver Tolerance = %"GSYM"\n", tol);
    printf("  timing output flag = %i\n\n", print_timing);
  }
  
  /* Create vectors */
  x = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                      ProbData.nprocs * ProbData.Nloc);
  if (check_flag(x, "N_VNew_Parallel", 0)) return 1;
  xhat = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                         ProbData.nprocs * ProbData.Nloc);
  if (check_flag(xhat, "N_VNew_Parallel", 0)) return 1;
  b = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                      ProbData.nprocs * ProbData.Nloc);
  if (check_flag(b, "N_VNew_Parallel", 0)) return 1;
  ProbData.d = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                               ProbData.nprocs * ProbData.Nloc);
  if (check_flag(ProbData.d, "N_VNew_Parallel", 0)) return 1;
  ProbData.s1 = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                                ProbData.nprocs * ProbData.Nloc);
  if (check_flag(ProbData.s1, "N_VNew_Parallel", 0)) return 1;
  ProbData.s2 = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                                ProbData.nprocs * ProbData.Nloc);
  if (check_flag(ProbData.s2, "N_VNew_Parallel", 0)) return 1;

  /* Fill xhat vector with uniform random data in [1,2] */
  vecdata = N_VGetArrayPointer(xhat);
  for (i=0; i<ProbData.Nloc; i++) 
    vecdata[i] = ONE + urand();

  /* Fill Jacobi vector with matrix diagonal */
  N_VConst(FIVE, ProbData.d);
  
  /* Create SPFGMR linear solver */
  LS = SUNLinSol_SPFGMR(x, PREC_RIGHT, maxl);
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_ITERATIVE,
                                 ProbData.myid);
  fails += Test_SUNLinSolSetATimes(LS, &ProbData, ATimes, ProbData.myid);
  fails += Test_SUNLinSolSetPreconditioner(LS, &ProbData, PSetup,
                                           PSolve, ProbData.myid);
  fails += Test_SUNLinSolSetScalingVectors(LS, ProbData.s1, ProbData.s2,
                                           ProbData.myid);
  fails += Test_SUNLinSolInitialize(LS, ProbData.myid);
  fails += Test_SUNLinSolSpace(LS, ProbData.myid);
  fails += SUNLinSol_SPFGMRSetGSType(LS, gstype);  
  if (fails) {
    printf("FAIL: SUNLinSol_SPFGMR module failed %i initialization tests\n\n", fails);
    return 1;
  } else if (ProbData.myid == 0)
    printf("SUCCESS: SUNLinSol_SPFGMR module passed all initialization tests\n\n");

  
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
  fails += SUNLinSol_SPFGMRSetPrecType(LS, PREC_NONE);  
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);
  
  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPFGMR module, problem 1, failed %i tests\n\n", fails);
    passfail += 1;
  } else if (ProbData.myid == 0) {
    printf("SUCCESS: SUNLinSol_SPFGMR module, problem 1, passed all tests\n\n");
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
  fails += SUNLinSol_SPFGMRSetPrecType(LS, PREC_RIGHT);  
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPFGMR module, problem 2, failed %i tests\n\n", fails);
    passfail += 1;
  } else if (ProbData.myid == 0) {
    printf("SUCCESS: SUNLinSol_SPFGMR module, problem 2, passed all tests\n\n");
  }

  
  /*** Test 3: Poisson-like solve w/ scaled rows (no preconditioning) ***/

  /* set scaling vectors */
  vecdata = N_VGetArrayPointer(ProbData.s1);
  for (i=0; i<ProbData.Nloc; i++)
    vecdata[i] = ONE + THOUSAND*urand();
  N_VConst(ONE, ProbData.s2);

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPFGMRSetPrecType(LS, PREC_NONE);  
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPFGMR module, problem 3, failed %i tests\n\n", fails);
    passfail += 1;
  } else if (ProbData.myid == 0) {
    printf("SUCCESS: SUNLinSol_SPFGMR module, problem 3, passed all tests\n\n");
  }

  
  /*** Test 4: Poisson-like solve w/ scaled rows (Jacobi preconditioning) ***/

  /* set scaling vectors */
  vecdata = N_VGetArrayPointer(ProbData.s1);
  for (i=0; i<ProbData.Nloc; i++)
    vecdata[i] = ONE + THOUSAND*urand();
  N_VConst(ONE, ProbData.s2);

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPFGMRSetPrecType(LS, PREC_RIGHT);  
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPFGMR module, problem 4, failed %i tests\n\n", fails);
    passfail += 1;
  } else if (ProbData.myid == 0) {
    printf("SUCCESS: SUNLinSol_SPFGMR module, problem 4, passed all tests\n\n");
  }

  
  /*** Test 5: Poisson-like solve w/ scaled columns (no preconditioning) ***/

  /* set scaling vectors */
  N_VConst(ONE, ProbData.s1);
  vecdata = N_VGetArrayPointer(ProbData.s2);
  for (i=0; i<ProbData.Nloc; i++)
    vecdata[i] = ONE + THOUSAND*urand();

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPFGMRSetPrecType(LS, PREC_NONE);
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPFGMR module, problem 5, failed %i tests\n\n", fails);
    passfail += 1;
  } else if (ProbData.myid == 0) {
    printf("SUCCESS: SUNLinSol_SPFGMR module, problem 5, passed all tests\n\n");
  }

  
  /*** Test 6: Poisson-like solve w/ scaled columns (Jacobi preconditioning) ***/

  /* set scaling vector, Jacobi solver vector */
  N_VConst(ONE, ProbData.s1);
  vecdata = N_VGetArrayPointer(ProbData.s2);
  for (i=0; i<ProbData.Nloc; i++)
    vecdata[i] = ONE + THOUSAND*urand();

  /* Fill x vector with scaled version */
  N_VDiv(xhat,ProbData.s2,x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) return 1;

  /* Run tests with this setup */
  fails += SUNLinSol_SPFGMRSetPrecType(LS, PREC_RIGHT);  
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol_SPFGMR module, problem 6, failed %i tests\n\n", fails);
    passfail += 1;
  } else if (ProbData.myid == 0) {
    printf("SUCCESS: SUNLinSol_SPFGMR module, problem 6, passed all tests\n\n");
  }

  /* check if any other process failed */
  (void) MPI_Allreduce(&passfail, &fails, 1, MPI_INT, MPI_MAX, ProbData.comm);
  
  /* Free solver and vectors */
  SUNLinSolFree(LS);
  N_VDestroy(x);
  N_VDestroy(xhat);
  N_VDestroy(b);
  N_VDestroy(ProbData.d);
  N_VDestroy(ProbData.s1);
  N_VDestroy(ProbData.s2);

  MPI_Finalize();
  return(fails);
}


/* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*/

/* matrix-vector product  */
int ATimes(void* Data, N_Vector v_vec, N_Vector z_vec)
{
  /* local variables */
  realtype *v, *z, *s1, *s2, vL, vR, vsL, vsR;
  sunindextype i, Nloc;
  int ierr;
  UserData *ProbData;
  MPI_Request SendReqL, SendReqR, RecvReqL, RecvReqR;
  MPI_Status stat;
  
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
  Nloc = ProbData->Nloc;
  
  /* MPI equivalent of realtype type */
  #if defined(SUNDIALS_SINGLE_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_FLOAT
  #elif defined(SUNDIALS_DOUBLE_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_DOUBLE
  #elif defined(SUNDIALS_EXTENDED_PRECISION)
  #define REALTYPE_MPI_TYPE MPI_LONG_DOUBLE
  #endif
  
  /* send/recv boundary data with neighbors */
  vL = vR = ZERO;
  vsL = v[0]*s2[0];
  vsR = v[Nloc-1]*s2[Nloc-1];
  if (ProbData->myid > 0) {                   /* left neighbor exists */
    ierr = MPI_Irecv(&vL, 1, REALTYPE_MPI_TYPE, ProbData->myid-1,
                     MPI_ANY_TAG, ProbData->comm, &RecvReqL);
    if (ierr != MPI_SUCCESS) return 1;
    ierr = MPI_Isend(&vsL, 1, REALTYPE_MPI_TYPE, ProbData->myid-1,
                     0, ProbData->comm, &SendReqL);
    if (ierr != MPI_SUCCESS) return 1;
  }
  if (ProbData->myid < ProbData->nprocs-1) {  /* right neighbor exists */
    ierr = MPI_Irecv(&vR, 1, REALTYPE_MPI_TYPE, ProbData->myid+1,
                     MPI_ANY_TAG, ProbData->comm, &RecvReqR);
    if (ierr != MPI_SUCCESS) return 1;
    ierr = MPI_Isend(&vsR, 1, REALTYPE_MPI_TYPE, ProbData->myid+1,
                     1, ProbData->comm, &SendReqR);
    if (ierr != MPI_SUCCESS) return 1;
  }
    
  
  /* iterate through interior of local domain, performing product */
  for (i=1; i<Nloc-1; i++) 
    z[i] = (-v[i-1]*s2[i-1] + FIVE*v[i]*s2[i] - v[i+1]*s2[i+1])/s1[i];

  /* wait on neighbor data to arrive */
  if (ProbData->myid > 0) {                   /* left neighbor exists */
    ierr = MPI_Wait(&RecvReqL, &stat);
    if (ierr != MPI_SUCCESS) return 1;
  }
  if (ProbData->myid < ProbData->nprocs-1) {  /* right neighbor exists */
    ierr = MPI_Wait(&RecvReqR, &stat);
    if (ierr != MPI_SUCCESS) return 1;
  }
  
  /* perform product at subdomain boundaries (note: vL/vR are zero at boundary)*/
  z[0] = (-vL + FIVE*v[0]*s2[0] - v[1]*s2[1])/s1[0];
  z[Nloc-1] = (-v[Nloc-2]*s2[Nloc-2] + FIVE*v[Nloc-1]*s2[Nloc-1] - vR)/s1[Nloc-1];
  
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
  for (i=0; i<ProbData->Nloc; i++) 
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
  int      failure = 0;
  sunindextype i;
  realtype *Xdata, *Ydata, maxerr;
  
  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);
  
  /* check vector data */
  for(i=0; i<local_problem_size; i++)
    failure += FNEQ(Xdata[i], Ydata[i], FIVE*tol*SUNRabs(Xdata[i]));

  if (failure > ZERO) {
    maxerr = ZERO;
    for(i=0; i < local_problem_size; i++)
      maxerr = SUNMAX(SUNRabs(Xdata[i]-Ydata[i])/SUNRabs(Xdata[i]), maxerr);
    printf("check err failure: maxerr = %g (tol = %g)\n",
	   maxerr, FIVE*tol);
    return(1);
  }
  else
    return(0);
}


