/*
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ SMU
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
 * These test functions are designed to check a SUNMatrix module
 * implementation.
 * -----------------------------------------------------------------
 */

#include <sundials/sundials_matrix.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <math.h> /* include isnan */
#include <stdio.h>
#include <stdlib.h>

#include "test_sunmatrix.h"

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif

/* private functions */
static double get_time();

int print_time = 0;
int print_all_ranks = 0;

#define FMT "%s Time: %22.15e\n\n"
#define PRINT_TIME(format, time) if(print_time) printf(format, time)

/* ----------------------------------------------------------------------
 * SUNMatGetID Test
 * --------------------------------------------------------------------*/
int Test_SUNMatGetID(SUNMatrix A, SUNMatrix_ID sunid, int myid)
{
  double       start_time, stop_time;
  SUNMatrix_ID mysunid;

  start_time = get_time();
  mysunid = SUNMatGetID(A);
  sync_device(A);
  stop_time  = get_time();

  if (sunid != mysunid) {
    TEST_STATUS(">>> FAILED test -- SUNMatGetID \n", myid);
    PRINT_TIME("    SUNMatGetID Time: %22.15e \n \n", stop_time - start_time);
    return(1);
  } else {
    TEST_STATUS("    PASSED test -- SUNMatGetID \n", myid);
  }

  if (myid == 0)
    PRINT_TIME("    SUNMatGetID Time: %22.15e \n \n", stop_time - start_time);

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNMatClone Test
 *
 * NOTE: This routine depends on SUNMatCopy to check matrix data.
 * --------------------------------------------------------------------*/
int Test_SUNMatClone(SUNMatrix A, int myid)
{
  int       failure;
  double    start_time, stop_time;
  realtype  tol=10*UNIT_ROUNDOFF;
  SUNMatrix B;

  /* clone matrix */
  start_time = get_time();
  B = SUNMatClone(A);
  sync_device(A);
  stop_time = get_time();

  /* check cloned matrix */
  if (B == NULL) {
    TEST_STATUS(">>> FAILED test -- SUNMatClone \n", myid);
    TEST_STATUS("    After SUNMatClone, B == NULL \n \n", myid);
    return(1);
  }

  /* check cloned matrix data */
  if (!has_data(B)) {
    TEST_STATUS(">>> FAILED test -- SUNMatClone \n", myid);
    TEST_STATUS("    Matrix data == NULL \n \n", myid);
    SUNMatDestroy(B);
    return(1);
  }

  failure = SUNMatCopy(A, B);
  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatClone (error using SUNMatCopy) \n", myid);
    SUNMatDestroy(B);
    return(1);
  }

  failure = check_matrix(A, B, tol);
  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatClone \n", myid);
    TEST_STATUS("    Failed SUNMatClone check \n \n", myid);
    SUNMatDestroy(B);
    return(1);
  } else {
    TEST_STATUS("    PASSED test -- SUNMatClone \n", myid);
  }

  if (myid == 0)
    PRINT_TIME("    SUNMatClone Time: %22.15e \n \n", stop_time - start_time);

  SUNMatDestroy(B);
  return(0);
}



/* ----------------------------------------------------------------------
 * SUNMatZero Test
 * --------------------------------------------------------------------*/
int Test_SUNMatZero(SUNMatrix A, int myid)
{
  int       failure;
  double    start_time, stop_time;
  realtype  tol=10*UNIT_ROUNDOFF;
  SUNMatrix B;

  /* protect A */
  B = SUNMatClone(A);

  /* set matrix data to zero */
  start_time = get_time();
  failure = SUNMatZero(B);
  sync_device(B);
  stop_time = get_time();

  if (failure) {
    TEST_STATUS2(">>> FAILED test -- SUNMatZero returned %d \n", failure, myid);
    SUNMatDestroy(B);
    return(1);
  }

  /* A data should be a vector of zeros */
  failure = check_matrix_entry(B, ZERO, tol);

  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatZero check \n", myid);
    PRINT_TIME("    SUNMatZero Time: %22.15e \n \n", stop_time - start_time);
    SUNMatDestroy(B);
    return(1);
  }
  else {
    TEST_STATUS("    PASSED test -- SUNMatZero \n", myid);
  }

  if (myid == 0)
    PRINT_TIME("    SUNMatZero Time: %22.15e \n \n", stop_time - start_time);

  SUNMatDestroy(B);
  return(0);
}


/* ----------------------------------------------------------------------
 * SUNMatCopy Test
 * --------------------------------------------------------------------*/
int Test_SUNMatCopy(SUNMatrix A, int myid)
{
  int       failure;
  double    start_time, stop_time;
  SUNMatrix B;
  realtype  tol=10*UNIT_ROUNDOFF;

  B = SUNMatClone(A);

  /* copy matrix data */
  start_time = get_time();
  failure = SUNMatCopy(A, B);
  sync_device(A);
  stop_time = get_time();

  if (failure) {
    TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned %d \n", failure, myid);
    SUNMatDestroy(B);
    return(1);
  }

  /* check matrix entries */
  failure = check_matrix(B, A, tol);

  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatCopy \n", myid);
    PRINT_TIME("    SUNMatCopy Time: %22.15e \n \n", stop_time - start_time);
    SUNMatDestroy(B);
    return(1);
  }
  else {
    TEST_STATUS("    PASSED test -- SUNMatCopy \n", myid);
  }

  if (myid == 0)
    PRINT_TIME("    SUNMatCopy Time: %22.15e \n \n", stop_time - start_time);

  SUNMatDestroy(B);
  return(0);
}



/* ----------------------------------------------------------------------
 * SUNMatScaleAdd Test: A = c * A + B
 *
 * NOTE: Sparse matrices will need additional testing for possibly
 * different sparsity patterns
 * --------------------------------------------------------------------*/
int Test_SUNMatScaleAdd(SUNMatrix A, SUNMatrix I, int myid)
{
  int       failure;
  double    start_time, stop_time;
  SUNMatrix B, C, D;
  realtype  tol=10*UNIT_ROUNDOFF;

  /*
   * Case 1: same sparsity/bandwith pattern
   */

  /* protect A */
  B = SUNMatClone(A);
  failure = SUNMatCopy(A, B);
  if (failure) {
    TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned %d \n", failure, myid);
    SUNMatDestroy(B);
    return(1);
  }

  /* fill vector data */
  start_time = get_time();
  failure = SUNMatScaleAdd(NEG_ONE, B, B);
  sync_device(B);
  stop_time = get_time();

  if (failure) {
    TEST_STATUS2(">>> FAILED test -- SUNMatScaleAdd returned %d \n", failure, myid);
    SUNMatDestroy(B);
    return(1);
  }

  /* check matrix entries */
  failure = check_matrix_entry(B, ZERO, tol);

  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatScaleAdd case 1 check  \n", myid);
    PRINT_TIME("    SUNMatScaleAdd Time: %22.15e \n \n", stop_time - start_time);
    SUNMatDestroy(B);
    return(1);
  }
  else {
    TEST_STATUS("    PASSED test -- SUNMatScaleAdd case 1 \n", myid);
  }

  if (myid == 0) {
    PRINT_TIME("    SUNMatScaleAdd Time: %22.15e \n \n", stop_time - start_time);
  }

  /*
   * Case 2: different sparsity/bandwidth patterns
   */
  if (is_square(A) && SUNMatGetID(A) != SUNMATRIX_CUSPARSE
      && SUNMatGetID(A) != SUNMATRIX_MAGMADENSE) {

    /* protect A and I */
    D = SUNMatClone(A);
    failure = SUNMatCopy(A, D);
    if (failure) {
      TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned %d \n", failure, myid);
      SUNMatDestroy(B);
      SUNMatDestroy(D);
      return(1);
    }
    C = SUNMatClone(I);
    failure = SUNMatCopy(I, C);
    if (failure) {
      TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned %d \n", failure, myid);
      SUNMatDestroy(B);
      SUNMatDestroy(C);
      SUNMatDestroy(D);
      return(1);
    }

    /* fill B and C */
    start_time = get_time();
    failure = SUNMatScaleAdd(ONE, D, I);
    if (failure) {
      TEST_STATUS2(">>> FAILED test -- SUNMatScaleAdd returned %d \n", failure, myid);
      SUNMatDestroy(B);
      SUNMatDestroy(C);
      SUNMatDestroy(D);
      return(1);
    }
    failure = SUNMatScaleAdd(ONE, C, A);
    if (failure) {
      TEST_STATUS2(">>> FAILED test -- SUNMatScaleAdd returned %d \n", failure, myid);
      SUNMatDestroy(B);
      SUNMatDestroy(C);
      SUNMatDestroy(D);
      return(1);
    }
    sync_device(A);
    stop_time = get_time();

    failure = check_matrix(D, C, tol);
    if (failure) {
      TEST_STATUS(">>> FAILED test -- SUNMatScaleAdd case 2 \n", myid);
      PRINT_TIME("    SUNMatScaleAdd Time: %22.15e \n \n", stop_time - start_time);
      SUNMatDestroy(B);
      SUNMatDestroy(C);
      SUNMatDestroy(D);
      return(1);
    }
    else {
      TEST_STATUS("    PASSED test -- SUNMatScaleAdd case 2 \n", myid);
    }

    if (myid == 0)
      PRINT_TIME("    SUNMatScaleAdd Time: %22.15e \n \n", stop_time - start_time);

    SUNMatDestroy(C);
    SUNMatDestroy(D);
  }

  SUNMatDestroy(B);

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNMatScaleAddI Tests
 *
 * NOTE: Sparse matrices will need additional testing for possibly
 * different sparsity patterns
 * --------------------------------------------------------------------*/
int Test_SUNMatScaleAddI(SUNMatrix A, SUNMatrix I, int myid)
{
  int       failure;
  double    start_time, stop_time;
  SUNMatrix B;
  realtype  tol=10*UNIT_ROUNDOFF;

  /* protect A */
  B = SUNMatClone(A);

  failure = SUNMatCopy(I, B);
  if (failure) {
    TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned %d \n", failure, myid);
    SUNMatDestroy(B);
    return(1);
  }

  /* perform operation */
  start_time = get_time();
  failure = SUNMatScaleAddI(NEG_ONE, B);
  sync_device(B);
  stop_time = get_time();

  if (failure) {
    TEST_STATUS2(">>> FAILED test -- SUNMatScaleAddI returned %d \n", failure, myid);
    SUNMatDestroy(B);
    return(1);
  }

  /* check matrix */
  failure = check_matrix_entry(B, ZERO, tol);

  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatScaleAddI check \n", myid);
    PRINT_TIME("    SUNMatScaleAddI Time: %22.15e \n \n", stop_time - start_time);
    SUNMatDestroy(B);
    return(1);
  }
  else {
    TEST_STATUS("    PASSED test -- SUNMatScaleAddI \n", myid);
  }

  if (myid == 0)
    PRINT_TIME("    SUNMatScaleAddI Time: %22.15e \n \n", stop_time - start_time);

  SUNMatDestroy(B);
  return(0);
}

/* ----------------------------------------------------------------------
 * SUNMatMatvecSetup Test
 * --------------------------------------------------------------------*/
int Test_SUNMatMatvecSetup(SUNMatrix A, int myid)
{
  int      failure;
  double   start_time, stop_time;

  if (A->ops->matvecsetup == NULL) {
    TEST_STATUS("    PASSED test -- SUNMatMatvecSetup not implemented\n", myid);
    return(0);
  }

  start_time = get_time();
  failure = SUNMatMatvecSetup(A);
  sync_device(A);
  stop_time = get_time();

  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatMatvecSetup \n", myid);
    return(1);
  } else {
    TEST_STATUS("    PASSED test -- SUNMatMatvecSetup\n", myid);
  }

  if (myid == 0)
    PRINT_TIME("    SUNMatMatvecSetup Time: %22.15e \n \n", stop_time - start_time);

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNMatMatvec Test (y should be correct A*x product)
 * --------------------------------------------------------------------*/
int Test_SUNMatMatvec(SUNMatrix A, N_Vector x, N_Vector y, int myid)
{
  int      failure;
  double   start_time, stop_time;
  SUNMatrix B;
  N_Vector  z, w;
  realtype  tol=100*UNIT_ROUNDOFF;

  if (A->ops->matvec == NULL) {
    TEST_STATUS("    PASSED test -- SUNMatMatvec not implemented\n", myid);
    return(0);
  }

  /* harder tests for square matrices */
  if (is_square(A)) {

    /* protect A */
    B = SUNMatClone(A);
    failure = SUNMatCopy(A, B);
    if (failure) {
      TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned %d \n", failure, myid);
      SUNMatDestroy(B);
      return(1);
    }

    /* compute matrix vector product */
    failure = SUNMatScaleAddI(THREE,B);
    if (failure) {
      TEST_STATUS2(">>> FAILED test -- SUNMatScaleAddI returned %d \n", failure, myid);
      SUNMatDestroy(B);
      return(1);
    }

    z = N_VClone(y); /* will be computed with matvec */
    w = N_VClone(y); /* will be the reference */

    /* Call the Setup function before the Matvec if it exists */
    if (B->ops->matvecsetup) {
      failure = SUNMatMatvecSetup(B);
      if (failure) {
        TEST_STATUS2(">>> FAILED test -- SUNMatMatvecSetup returned %d \n", failure, myid);
        SUNMatDestroy(B);
        return(1);
      }
    }

    start_time = get_time();
    failure = SUNMatMatvec(B,x,z); /* z = (3A+I)x = 3y + x */
    sync_device(B);
    stop_time = get_time();

    if (failure) {
      TEST_STATUS2(">>> FAILED test -- SUNMatMatvec returned %d \n", failure, myid);
      SUNMatDestroy(B);
      return(1);
    }

    N_VLinearSum(THREE,y,ONE,x,w); /* w = 3y + x */

    failure = check_vector(w,z,tol);

    SUNMatDestroy(B);
    N_VDestroy(z);
    N_VDestroy(w);

  } else {

    z = N_VClone(y); /* will be computed with matvec */

    start_time = get_time();
    failure = SUNMatMatvec(A,x,z); /* z = Ax */
    sync_device(A);
    stop_time = get_time();

    if (failure) {
      TEST_STATUS2(">>> FAILED test -- SUNMatMatvec returned %d \n", failure, myid);
      return(1);
    }

    failure = check_vector(y,z,tol);

    N_VDestroy(z);

  }

  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatMatvec check \n", myid);
    PRINT_TIME("    SUNMatMatvec Time: %22.15e \n \n", stop_time - start_time);
    return(1);
  }
  else {
    TEST_STATUS("    PASSED test -- SUNMatMatvec \n", myid);
  }

  if (myid == 0)
    PRINT_TIME("    SUNMatMatvec Time: %22.15e \n \n", stop_time - start_time);

  return(0);
}



/* ----------------------------------------------------------------------
 * SUNMatSpace Test
 * --------------------------------------------------------------------*/
int Test_SUNMatSpace(SUNMatrix A, int myid)
{
  int      failure;
  double   start_time, stop_time;
  long int lenrw, leniw;

  start_time = get_time();
  failure = SUNMatSpace(A, &lenrw, &leniw);
  sync_device(A);
  stop_time = get_time();

  if (failure) {
    TEST_STATUS(">>> FAILED test -- SUNMatSpace \n", myid);
    PRINT_TIME("    SUNMatSpace Time: %22.15e \n \n", stop_time - start_time);
    return(1);
  } else {
    TEST_STATUS2("    PASSED test -- SUNMatSpace lenrw=%li ", lenrw, myid);
    TEST_STATUS2("leniw=%li\n", leniw, myid);
  }

  if (myid == 0)
    PRINT_TIME("    SUNMatSpace Time: %22.15e \n \n", stop_time - start_time);

  return(0);
}


/* ======================================================================
 * Private functions
 * ====================================================================*/

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
time_t base_time_tv_sec = 0; /* Base time; makes time values returned
                                by get_time easier to read when
                                printed since they will be zero
                                based.
                              */
#endif

void SetTiming(int onoff)
{
   print_time = onoff;

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime( CLOCK_MONOTONIC, &spec );
  base_time_tv_sec = spec.tv_sec;
#endif
}

void SetPrintAllRanks(int onoff)
{
  print_all_ranks = onoff;
}

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
static double get_time()
{
#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime( CLOCK_MONOTONIC, &spec );
  double time = (double)(spec.tv_sec - base_time_tv_sec) + ((double)(spec.tv_nsec) / 1E9);
#else
  double time = 0;
#endif
  return time;
}
