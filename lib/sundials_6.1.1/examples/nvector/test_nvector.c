/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner and Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: These testing routines are based on an
 *                   NVECTOR testing routine by Daniel R. Reynolds
 *                   @ SMU.
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
 * These test functions are designed to check an NVECTOR module
 * implementation.
 * -----------------------------------------------------------------*/

/* Minimum POSIX version needed for struct timespec and clock_monotonic */
#if !defined(_POSIX_C_SOURCE) || (_POSIX_C_SOURCE < 199309L)
#define _POSIX_C_SOURCE 199309L
#endif

#include <sundials/sundials_config.h>

/* POSIX timers */
#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
#include <time.h>
#include <stddef.h>
#include <unistd.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM ".17f"
#endif

/* all tests need a SUNContext */
SUNContext sunctx = NULL;

#if defined(SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
static time_t base_time_tv_sec = 0; /* Base time; makes time values returned
                                       by get_time easier to read when
                                       printed since they will be zero
                                       based.
                                    */
#endif

/* private functions */
static double get_time();

/* private variables */
static int print_time = 0;

/* macro for printing timings */
#define FMT "%s Time: %22.15e\n\n"
#define PRINT_TIME(test, time) if (print_time) printf(FMT, test, time)

int Test_Init(void* comm)
{
  if (sunctx == NULL) {
    if (SUNContext_Create(comm, &sunctx)) {
      printf("ERROR: SUNContext_Create failed\n");
      return -1;
    }
  }
  return 0;
}

int Test_Finalize()
{
  if (sunctx != NULL) {
    if (SUNContext_Free(&sunctx)) {
      printf("ERROR: SUNContext_Create failed\n");
      return -1;
    }
  }
  return 0;
}

void Test_Abort(int code)
{
  Test_Finalize();
  abort();
}

/* ----------------------------------------------------------------------
 * N_VMake Test
 * --------------------------------------------------------------------*/
int Test_N_VMake(N_Vector X, sunindextype local_length, int myid)
{
  /* check if vector is NULL */
  if (X == NULL) {
    printf(">>> FAILED test -- N_VMake, Proc %d \n", myid);
    printf("    After N_VMakeEmpty, X == NULL \n \n");
    return(1);
  }

  /* check for vector data */
  if (!has_data(X)) {
    printf(">>> FAILED test -- N_VMake, Proc %d \n", myid);
    printf("    Vector data == NULL \n \n");
    return(1);
  }

  if (myid == 0)
    printf("PASSED test -- N_VMake \n");

  return(0);
}

/* ----------------------------------------------------------------------
 * N_VGetVectorID Test
 * --------------------------------------------------------------------*/
int Test_N_VGetVectorID(N_Vector X, N_Vector_ID ID, int myid)
{
  if (N_VGetVectorID(X) != ID) {
    printf(">>> FAILED test -- N_VGetVectorID, Proc %d \n", myid);
    printf("    Unrecognized vector type %d \n \n", N_VGetVectorID(X));
    return (1);
  } else if (myid == 0) {
    printf("PASSED test -- N_VGetVectorID \n");
  }
  return (0);
}

/* ----------------------------------------------------------------------
 * N_VCloneVectorArray Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VCloneVectorArray(int count, N_Vector W,
                             sunindextype local_length, int myid)
{
  int      i, failure;
  double   start_time, stop_time, maxt;
  N_Vector *vs;

  /* check if the required operations are implemented */
  if (W->ops->nvconst == NULL) {
    printf(">>> FAILED test -- N_VCloneVectorArray, Proc %d missing required operations\n",
           myid);
    return(1);
  }

  /* clone array of vectors */
  start_time = get_time();
  vs = N_VCloneVectorArray(count, W);
  stop_time = get_time();

  /* check array of vectors */
  if (count <= 0 && vs != NULL) {
    printf(">>> FAILED test -- N_VCloneVectorArray, Proc %d \n", myid);
    printf("    count = %d, expected *vs = NULL \n\n",count);
    return(1);
  }

  /* check vectors in array */
  for(i=0; i<count; i++) {
    if (vs[i] == NULL) {
      printf(">>> FAILED test -- N_VCloneVectorArray, Proc %d \n", myid);
      printf("    Vector[%d] = NULL \n\n",i);
      N_VDestroyVectorArray(vs, count);
      return(1);
    }

    N_VConst(ONE,vs[i]);
    failure = check_ans(ONE, vs[i], local_length);
    if (failure) {
      printf(">>> FAILED test -- N_VCloneVectorArray, Proc %d \n", myid);
      printf("    Vector[%d] failed N_VConst check \n\n",i);
      N_VDestroyVectorArray(vs, count);
      return(1);
    }
  }

  N_VDestroyVectorArray(vs, count);

  if (myid == 0)
    printf("PASSED test -- N_VCloneVectorArray \n");

  /* find max time across all processes */
  maxt = max_time(W, stop_time - start_time);
  PRINT_TIME("N_VCloneVectorArray", maxt);

  return(0);
}

/* ----------------------------------------------------------------------
 * N_VCloneVectorArrayEmpty Test
 * --------------------------------------------------------------------*/
int Test_N_VCloneEmptyVectorArray(int count, N_Vector W, int myid)
{
  int      i;
  double   start_time, stop_time, maxt;
  N_Vector *vs;

  /* clone empty array */
  start_time = get_time();
  vs = N_VCloneEmptyVectorArray(count, W);
  stop_time = get_time();

  /* check array of vectors */
  if (count <= 0 && vs != NULL) {
    printf(">>> FAILED test -- N_VCloneEmptyVectorArray, Proc %d \n", myid);
    printf("    count = %d, expected *vs = NULL \n\n",count);
    return(1);
  }

  /* check vectors in array */
  for(i=0; i<count; i++) {
    if (vs[i] == NULL) {
      printf(">>> FAILED test -- N_VCloneEmptyVectorArray, Proc %d \n", myid);
      printf("    Vector[%d] = NULL \n\n",i);
      N_VDestroyVectorArray(vs, count);
      return(1);
    }

    if (has_data(vs[i])) {
      printf(">>> FAILED test -- N_VCloneEmptyVectorArray, Proc %d \n", myid);
      printf("    Vector[%d] data != NULL \n\n",i);
      N_VDestroyVectorArray(vs, count);
      return(1);
    }
  }

  N_VDestroyVectorArray(vs, count);

  if (myid == 0)
    printf("PASSED test -- N_VCloneEmptyVectorArray \n");

  /* find max time across all processes */
  maxt = max_time(W, stop_time - start_time);
  PRINT_TIME("N_VCloneEmptyVectorArray", maxt);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VCloneEmpty Test
 * --------------------------------------------------------------------*/
int Test_N_VCloneEmpty(N_Vector W, int myid)
{
  double   start_time, stop_time, maxt;
  N_Vector X;

  /* clone empty vector */
  start_time = get_time();
  X = N_VCloneEmpty(W);
  stop_time = get_time();

  /* check vector */
  if (X == NULL) {
    printf(">>> FAILED test -- N_VCloneEmpty, Proc %d \n", myid);
    printf("    After N_VCloneEmpty, X == NULL \n\n");
    return(1);
  }

  /* check vector data */
  if (has_data(X)) {
    printf(">>> FAILED test -- N_VCloneEmpty, Proc %d \n", myid);
    printf("    Vector data != NULL \n\n");
    N_VDestroy(X);
    return(1);
  }

  N_VDestroy(X);

  if (myid == 0)
    printf("PASSED test -- N_VCloneEmpty \n");

  /* find max time across all processes */
  maxt = max_time(W, stop_time - start_time);
  PRINT_TIME("N_VCloneEmpty", maxt);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VClone Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VClone(N_Vector W, sunindextype local_length, int myid)
{
  int      failure;
  double   start_time, stop_time, maxt;
  N_Vector X;

  /* check if the required operations are implemented */
  if (W->ops->nvconst == NULL) {
    printf(">>> FAILED test -- N_VClone, Proc %d missing required operations\n",
           myid);
    return(1);
  }

  /* clone vector */
  start_time = get_time();
  X = N_VClone(W);
  stop_time = get_time();

  /* check cloned vector */
  if (X == NULL) {
    printf(">>> FAILED test -- N_VClone, Proc %d \n", myid);
    printf("    After N_VClone, X == NULL \n\n");
    return(1);
  }

  /* check cloned vector data */
  if (!has_data(X)) {
    printf(">>> FAILED test -- N_VClone, Proc %d \n", myid);
    printf("    Vector data == NULL \n\n");
    N_VDestroy(X);
    return(1);
  }

  N_VConst(ONE,X);
  failure = check_ans(ONE, X, local_length);
  if (failure) {
    printf(">>> FAILED test -- N_VClone, Proc %d \n", myid);
    printf("    Failed N_VClone check \n\n");
    N_VDestroy(X);
    return(1);
  }
  N_VDestroy(X);

  if (myid == 0)
    printf("PASSED test -- N_VClone \n");

  /* find max time across all processes */
  maxt = max_time(W, stop_time - start_time);
  PRINT_TIME("N_VClone", maxt);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VGetArrayPointer Test
 *
 * For now commenting this out to surpress warning messages (pointer set,
 * but not used). Do we really need to time access to the vector
 * data pointer?
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VGetArrayPointer(N_Vector W, sunindextype local_length, int myid)
{
  int      failure = 0;
  double   start_time, stop_time, maxt;
  realtype *Wdata;

  /* check if the required operations are implemented */
  if (W->ops->nvconst == NULL) {
    printf(">>> FAILED test -- N_VGetArrayPointer, Proc %d missing required operations\n",
           myid);
    return(1);
  }

  /* get vector data, time it and set it to NULL */
  start_time = get_time();
  Wdata = N_VGetArrayPointer(W);
  stop_time = get_time();
  Wdata++; Wdata=NULL; /* Do something with pointer to surpress warning */

  /* check vector data */
  if (!has_data(W)) {
    printf(">>> FAILED test -- N_VGetArrayPointer, Proc %d \n", myid);
    printf("    Vector data == NULL \n\n");
    return(1);
  }

  N_VConst(NEG_HALF,W);
  failure = check_ans(NEG_HALF, W, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VGetArrayPointer, Proc %d \n", myid);
    printf("    Failed N_VConst check \n\n");
    return(1);
  }

  if (myid == 0)
    printf("PASSED test -- N_VGetArrayPointer \n");

  /* find max time across all processes */
  maxt = max_time(W, stop_time - start_time);
  PRINT_TIME("N_VGetArrayPointer", maxt);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VSetArrayPointer Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VSetArrayPointer(N_Vector W, sunindextype local_length, int myid)
{
  int          failure = 0;
  sunindextype i;
  double       start_time, stop_time, maxt;
  realtype     *Wdata;

  /* check if the required operations are implemented */
  if (W->ops->nvconst == NULL) {
    printf(">>> FAILED test -- N_VSetArrayPointer, Proc %d missing required operations\n",
           myid);
    return(1);
  }

  /* create vector data */
  Wdata = (realtype*) malloc(local_length * sizeof(realtype));
  for(i=0; i < local_length; i++){
    Wdata[i] = ONE;
  }

  /* attach data to vector */
  start_time = get_time();
  N_VSetArrayPointer(Wdata, W);
  stop_time = get_time();

  /* check vector data */
  N_VConst(NEG_HALF,W);
  for(i=0; i < local_length; i++){
    failure += SUNRCompare(Wdata[i], NEG_HALF);
  }

  if (failure) {
    printf(">>> FAILED test -- N_VSetArrayPointer, Proc %d \n", myid);
    printf("    Failed N_VConst check \n\n");
    free(Wdata);
    return(1);
  }

  free(Wdata);

  if (myid == 0)
    printf("PASSED test -- N_VSetArrayPointer \n");

  /* find max time across all processes */
  maxt = max_time(W, stop_time - start_time);
  PRINT_TIME("N_VSetArrayPointer", maxt);

  return(0);
}


/* ----------------------------------------------------------------------
 * Test_N_VGetLength Test
 *
 * NOTE: This routine depends on N_VConst and N_VDotProd.  It also
 * modifies the values inside the input vector W.
 * --------------------------------------------------------------------*/
int Test_N_VGetLength(N_Vector W, int myid)
{
  sunindextype Wlength, Wlength2;

  /* check if the required operations are implemented */
  if (W->ops->nvconst == NULL || W->ops->nvdotprod == NULL) {
    printf(">>> FAILED test -- N_VGetLength, Proc %d missing required operations\n",
           myid);
    return(1);
  }

  /* ask W for it's overall length */
  Wlength = N_VGetLength(W);

  /* use N_VConst and N_VDotProd to compute length */
  N_VConst(RCONST(1.0), W);
  Wlength2 = (sunindextype) N_VDotProd(W, W);
  sync_device(W);

  /* return error if lengths disagree */
  if (Wlength != Wlength2) {
    printf(">>> FAILED test -- N_VGetLength, Proc %d (%li != %li)\n", myid,
           (long int) Wlength, (long int) Wlength2);
    return(1);
  }
  if (myid == 0)
    printf("PASSED test -- N_VGetLength\n");
  return(0);
}


/* ----------------------------------------------------------------------
 * Test_N_VGetCommunicator Test (without MPI dependency)
 * --------------------------------------------------------------------*/
int Test_N_VGetCommunicator(N_Vector W, void *comm, int myid)
{
  void* wcomm;

  /* ask W for its communicator */
  wcomm = NULL;
  wcomm = N_VGetCommunicator(W);

  /* return with success if both are NULL */
  if ((wcomm == NULL) && (comm == NULL))  {
    printf("PASSED test -- N_VGetCommunicator\n");
    return(0);
  }

  /* return with failure if either is NULL */
  if (wcomm == NULL) {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (incorrectly reports NULL comm)\n", myid);
    return(1);
  }
  if (comm == NULL) {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (incorrectly reports non-NULL comm)\n", myid);
    return(1);
  }

  if (myid == 0)
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d has non-NULL comm with MPI disabled\n", myid);
  return(0);
}


/* ----------------------------------------------------------------------
 * N_VLinearSum Tests
 * --------------------------------------------------------------------*/
int Test_N_VLinearSum(N_Vector X, N_Vector Y, N_Vector Z,
                      sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;

  /*
   * Case 1a: y = x + y, (Vaxpy Case 1)
   */

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(NEG_TWO, Y);

  start_time = get_time();
  N_VLinearSum(ONE, X, ONE, Y, Y);
  sync_device(X);
  stop_time = get_time();

  /* Y should be vector of -1 */
  failure = check_ans(NEG_ONE, Y, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 1b: y = -x + y, (Vaxpy Case 2)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(TWO, Y);

  start_time = get_time();
  N_VLinearSum(NEG_ONE, X, ONE, Y, Y);
  sync_device(X);
  stop_time = get_time();

  /* Y should be vector of +1 */
  failure = check_ans(ONE, Y, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 1c: y = ax + y, (Vaxpy Case 3)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_TWO, Y);

  start_time = get_time();
  N_VLinearSum(HALF, X, ONE, Y, Y);
  sync_device(X);
  stop_time = get_time();

  /* Y should be vector of -1 */
  failure = check_ans(NEG_ONE, Y, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 1c, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 1c \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum Time", maxt);

  /*
   * Case 2a: x = x + y, (Vaxpy Case 1)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_ONE, Y);

  start_time = get_time();
  N_VLinearSum(ONE, X, ONE, Y, X);
  sync_device(X);
  stop_time = get_time();

  /* X should be vector of +1 */
  failure = check_ans(ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 2a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 2a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 2b: x = x - y, (Vaxpy Case 2)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(TWO, Y);

  start_time = get_time();
  N_VLinearSum(ONE, X, NEG_ONE, Y, X);
  sync_device(X);
  stop_time = get_time();

  /* X should be vector of -1 */
  failure = check_ans(NEG_ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 2b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 2b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 2c: x = x + by, (Vaxpy Case 3)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_HALF, Y);

  start_time = get_time();
  N_VLinearSum(ONE, X, TWO, Y, X);
  sync_device(X);
  stop_time = get_time();

  /* X should be vector of +1 */
  failure = check_ans(ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 2c, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 2c \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 3: z = x + y, (VSum)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(ONE, X, ONE, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 3, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 3 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 4a: z = x - y, (VDiff)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(ONE, X, NEG_ONE, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 4a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 4a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 4b: z = -x + y, (VDiff)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(NEG_ONE, X, ONE, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 4b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 4b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 5a: z = x + by, (VLin1)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(ONE, X, TWO, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 5a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 5a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 5b: z = ax + y, (VLin1)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(HALF, X);
  N_VConst(NEG_TWO, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(TWO, X, ONE, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 5b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 5b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 6a: z = -x + by, (VLin2)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(NEG_ONE, X, TWO, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 6a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 6a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum Time", maxt);

  /*
   * Case 6b: z = ax - y, (VLin2)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(HALF, X);
  N_VConst(TWO, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(TWO, X, NEG_ONE, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 6b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 6b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 7: z = a(x + y), (VScaleSum)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(TWO, X, TWO, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 7, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 7 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 8: z = a(x - y), (VScaleDiff)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(HALF, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(TWO, X, NEG_TWO, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 8, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 8 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  /*
   * Case 9: z = ax + by, All Other Cases
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(NEG_TWO, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VLinearSum(TWO, X, HALF, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 9, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSum Case 9 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearSum", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VConst Test
 * --------------------------------------------------------------------*/
int Test_N_VConst(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;

  /* fill vector data with zeros to prevent passing in the case where
     the input vector is a vector of ones */
  set_element_range(X, 0, local_length-1, ZERO);

  start_time = get_time();
  N_VConst(ONE,X);
  sync_device(X);
  stop_time = get_time();

  /* X should be vector of +1 */
  failure = check_ans(ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VConst, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VConst \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VConst", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VProd Test
 * --------------------------------------------------------------------*/
int Test_N_VProd(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VProd(X, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VProd, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VProd \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VProd", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VDiv Test
 * --------------------------------------------------------------------*/
int Test_N_VDiv(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(TWO, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VDiv(X, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1/2 */
  failure = check_ans(HALF, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VDiv, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDiv \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDiv", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VScale Tests
 * --------------------------------------------------------------------*/
int Test_N_VScale(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;

  /*
   * Case 1: x = cx, VScaleBy
   */

  /* fill vector data */
  N_VConst(HALF, X);

  start_time = get_time();
  N_VScale(TWO, X, X);
  sync_device(X);
  stop_time = get_time();

  /* X should be vector of +1 */
  failure = check_ans(ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VScale Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScale Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScale", maxt);

  /*
   * Case 2: z = x, VCopy
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_ONE, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VScale(ONE, X, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VScale Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScale Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScale", maxt);

  /*
   * Case 3: z = -x, VNeg
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_ONE, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VScale(NEG_ONE, X, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VScale Case 3, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScale Case 3 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScale", maxt);

  /*
   * Case 4: z = cx, All other cases
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VScale(TWO, X, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VScale Case 4, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScale Case 4 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScale", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VAbs Test
 * --------------------------------------------------------------------*/
int Test_N_VAbs(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;

  /* fill vector data */
  N_VConst(NEG_ONE, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VAbs(X,Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VAbs, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VAbs \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VAbs", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VInv Test
 * --------------------------------------------------------------------*/
int Test_N_VInv(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VInv(X,Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of +1/2 */
  failure = check_ans(HALF, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VInv, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VInv \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VInv", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VAddConst Test
 * --------------------------------------------------------------------*/
int Test_N_VAddConst(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VAddConst(X,NEG_TWO,Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VAddConst, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VAddConst \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VAddConst", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VDotProd Test
 * --------------------------------------------------------------------*/
int Test_N_VDotProd(N_Vector X, N_Vector Y, sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype global_length;
  realtype     ans;

  /* get global length */
  global_length = N_VGetLength(X);

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(HALF, Y);

  start_time = get_time();
  ans = N_VDotProd(X,Y);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal global vector length */
  failure = SUNRCompare(ans, (realtype) global_length);

  if (failure) {
    printf(">>> FAILED test -- N_VDotProd, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProd \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProd", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMaxNorm Test
 * --------------------------------------------------------------------*/
int Test_N_VMaxNorm(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans;

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  if (myid == 0)
    set_element(X, local_length-1, NEG_TWO);
  else
    set_element(X, local_length-1, ONE);

  start_time = get_time();
  ans = N_VMaxNorm(X);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal 2 */
  failure = (ans < ZERO) ? 1 : SUNRCompare(ans, TWO);

  if (failure) {
    printf(">>> FAILED test -- N_VMaxNorm, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMaxNorm \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VMaxNorm", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNorm Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNorm(N_Vector X, N_Vector W, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans;

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(HALF, W);

  start_time = get_time();
  ans = N_VWrmsNorm(X, W);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal 1/4 */
  failure = (ans < ZERO) ? 1 : SUNRCompare(ans, HALF*HALF);

  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNorm, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWrmsNorm \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWrmsNorm", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNormMask Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNormMask(N_Vector X, N_Vector W, N_Vector ID,
                         sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype global_length;
  realtype     ans;
  realtype     fac;

  /* get global length */
  global_length = N_VGetLength(X);

  /* factor used in checking solutions */
  fac = SUNRsqrt((realtype) (global_length - 1)/(global_length))*HALF*HALF;

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(HALF, W);

  /* use all elements except one */
  N_VConst(ONE, ID);
  if (myid == 0)
    set_element(ID, local_length-1, ZERO);

  start_time = get_time();
  ans = N_VWrmsNormMask(X, W, ID);
  sync_device(X);
  stop_time = get_time();

  /* check ans */
  failure = (ans < ZERO) ? 1 : SUNRCompare(ans, fac);

  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNormMask, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWrmsNormMask \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWrmsNormMask", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMin Test
 * --------------------------------------------------------------------*/
int Test_N_VMin(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans;

  /* fill vector data */
  N_VConst(TWO, X);
  if (myid == 0)
    set_element(X, local_length-1, HALF);
  else
    set_element(X, local_length-1, ONE);

  start_time = get_time();
  ans = N_VMin(X);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal 0.5 */
  failure = SUNRCompare(ans, HALF);

  if (failure) {
    printf(">>> FAILED test -- N_VMin Case 1, Proc %d \n", myid);
    printf("    min = %" FSYM ", expected %" FSYM "\n", ans, HALF);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMin \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VMin", maxt);

  /* fill vector data */
  N_VConst(TWO, X);
  if (myid == 0)
    set_element(X, local_length-1, NEG_TWO);
  else
    set_element(X, local_length-1, NEG_ONE);

  start_time = get_time();
  ans = N_VMin(X);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal -2 */
  failure = SUNRCompare(ans, NEG_TWO);

  if (failure) {
    printf(">>> FAILED test -- N_VMin, Proc %d \n", myid);
    printf("    min = %" FSYM ", expected %" FSYM "\n", ans, NEG_TWO);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMin \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VMin", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWL2Norm Test
 * --------------------------------------------------------------------*/
int Test_N_VWL2Norm(N_Vector X, N_Vector W, sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype global_length;
  realtype     ans;

  /* get global length */
  global_length = N_VGetLength(X);

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(HALF, W);

  start_time = get_time();
  ans = N_VWL2Norm(X, W);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal 1/4 * sqrt(global_length) */
  failure = (ans < ZERO) ? 1 : SUNRCompare(ans, HALF*HALF*SUNRsqrt((realtype) global_length));

  if (failure) {
    printf(">>> FAILED test -- N_VWL2Norm, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWL2Norm \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWL2Norm", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VL1Norm Test
 * --------------------------------------------------------------------*/
int Test_N_VL1Norm(N_Vector X, sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype global_length;
  realtype     ans;

  /* get global length */
  global_length = N_VGetLength(X);

  /* fill vector data */
  N_VConst(NEG_ONE, X);

  start_time = get_time();
  ans = N_VL1Norm(X);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal global_length */
  failure = (ans < ZERO) ? 1 : SUNRCompare(ans, (realtype) global_length);

  if (failure) {
    printf(">>> FAILED test -- N_VL1Norm, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VL1Norm \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VL1Norm", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VCompare
 * --------------------------------------------------------------------*/
int Test_N_VCompare(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype i;

  if (local_length < 3) {
    printf("Error Test_N_VCompare: Local vector length is %ld, length must be >= 3\n",
           (long int) local_length);
    return(1);
  }

  /* fill vector data */
  for(i=0; i < local_length; i++){
    set_element(Z, i, NEG_ONE);

    switch(i % 3) {

    case 0 :
      /* abs(X[i]) < c */
      set_element(X, i, ZERO);
      break;

    case 1 :
      /* abs(X[i]) = c */
      set_element(X, i, NEG_ONE);
      break;

    case 2 :
      /* abs(X[i]) > c */
      set_element(X, i, NEG_TWO);
      break;
    }
  }

  start_time = get_time();
  N_VCompare(ONE, X, Z);
  sync_device(X);
  stop_time = get_time();

  /* check return vector */
  for(i=0; i < local_length; i++){

    switch(i % 3) {

    case 0 :
      /* Z[i] == 0 */
      if (get_element(Z, i) != ZERO)
        failure = 1;
      break;

    case 1 :
      /* Z[i] == 1 */
      if (get_element(Z, i) != ONE)
        failure = 1;
      break;

    case 2 :
      /* Z[i] == 1 */
      if (get_element(Z, i) != ONE)
        failure = 1;
      break;
    }
  }

  if (failure) {
    printf(">>> FAILED test -- N_VCompare, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VCompare \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VCompare", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VInvTest
 * --------------------------------------------------------------------*/
int Test_N_VInvTest(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype i;
  booleantype  ans, exp;

  if (local_length < 2) {
    printf("Error Test_N_VInvTest: Local vector length is %ld, length must be >= 2\n",
           (long int) local_length);
    return(1);
  }

  /*
   * Case 1: All elements Nonzero, z[i] = 1/x[i], return True
   */

  /* fill vector data */
  N_VConst(HALF, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  ans = N_VInvTest(X, Z);
  sync_device(X);
  stop_time = get_time();

  /* we expect no zeros */
  exp = SUNTRUE;

  /* Z should be vector of +2 */
  failure = check_ans(TWO, Z, local_length);

  if (failure || (ans != exp)) {
    printf(">>> FAILED test -- N_VInvTest Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VInvTest Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VInvTest", maxt);

  /*
   * Case 2: Some elements Zero, z[i] = 1/x[i] for x[i] != 0, return False
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ZERO, Z);
  for(i=0; i < local_length; i++){
    if (i % 2) {
      set_element(X, i, HALF);
    } else {
      exp = SUNFALSE;
      set_element(X, i, ZERO);
    }
  }

  start_time = get_time();
  ans = N_VInvTest(X, Z);
  sync_device(X);
  stop_time = get_time();

  /* check return vector */
  for(i=0; i < local_length; i++){
    if (i % 2) {
      if (get_element(Z, i) != TWO)
        failure = 1;
    } else {
      if (get_element(Z, i) != ZERO)
        failure = 1;
    }
  }

  if (failure || (ans != exp)) {
    printf(">>> FAILED test -- N_VInvTest Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VInvTest Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VInvTest", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VConstrMask
 * --------------------------------------------------------------------*/
int Test_N_VConstrMask(N_Vector C, N_Vector X, N_Vector M,
                       sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype i;
  booleantype  test;

  if (local_length < 7) {
    printf("Error Test_N_VConstrMask: Local vector length is %ld, length must be >= 7\n",
           (long int) local_length);
    return(1);
  }

  /*
   * Case 1: Return True
   */

  /* fill vector data */
  for(i=0; i < local_length; i++){
    set_element(M, i, NEG_ONE);

    switch(i % 7) {
    case 0 :
      /* c = -2, test for < 0*/
      set_element(C, i, NEG_TWO);
      set_element(X, i, NEG_TWO);
      break;

    case 1 :
      /* c = -1, test for <= 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, NEG_ONE);
      break;

    case 2 :
      /* c = -1, test for == 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, ZERO);
      break;

    case 3 :
      /* c = 0, no test */
      set_element(C, i, ZERO);
      set_element(X, i, HALF);
      break;

    case 4 :
      /* c = 1, test for == 0*/
      set_element(C, i, ONE);
      set_element(X, i, ZERO);
      break;

    case 5 :
      /* c = 1, test for >= 0*/
      set_element(C, i, ONE);
      set_element(X, i, ONE);
      break;

    case 6:
      /* c = 2, test for > 0 */
      set_element(C, i, TWO);
      set_element(X, i, TWO);
      break;
    }
  }

  start_time = get_time();
  test = N_VConstrMask(C, X, M);
  sync_device(X);
  stop_time = get_time();

  /* M should be vector of 0 */
  failure = check_ans(ZERO, M, local_length);

  if (failure || !test) {
    printf(">>> FAILED test -- N_VConstrMask Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VConstrMask Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VConstrMask", maxt);

  /*
   * Case 2: Return False
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  for(i=0; i < local_length; i++){
    set_element(M, i, NEG_ONE);

    switch(i % 5) {
    case 0 :
      /* c = -2, test for < 0*/
      set_element(C, i, NEG_TWO);
      set_element(X, i, TWO);
      break;

    case 1 :
      /* c = -1, test for <= 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, ONE);
      break;

    case 2 :
      /* c = 0, no test */
      set_element(C, i, ZERO);
      set_element(X, i, HALF);
      break;

    case 3 :
      /* c = 1, test for >= 0*/
      set_element(C, i, ONE);
      set_element(X, i, NEG_ONE);
      break;

    case 4 :
      /* c = 2, test for > 0 */
      set_element(C, i, TWO);
      set_element(X, i, NEG_TWO);
      break;
    }
  }

  start_time = get_time();
  test = N_VConstrMask(C, X, M);
  sync_device(X);
  stop_time = get_time();

  /* check mask vector */
  for(i=0; i < local_length; i++){
    if ((i % 5) == 2){
      if (get_element(M, i) != ZERO)
        failure = 1;
    } else {
      if (get_element(M, i) != ONE)
        failure = 1;
    }
  }

  if (failure || test) {
    printf(">>> FAILED test -- N_VConstrMask Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VConstrMask Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VConstrMask", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMinQuotient Test
 * --------------------------------------------------------------------*/
int Test_N_VMinQuotient(N_Vector NUM, N_Vector DENOM,
                        sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans;

  /*
   * Case 1: Pass
   */

  /* fill vector data */
  N_VConst(TWO, NUM);
  N_VConst(TWO, DENOM);
  if (myid == 0)
    set_element(NUM, local_length-1, HALF);
  else
    set_element(NUM, local_length-1, ONE);

  start_time = get_time();
  ans = N_VMinQuotient(NUM, DENOM);
  sync_device(NUM);
  stop_time = get_time();

  /* ans should equal 1/4 */
  failure = SUNRCompare(ans, HALF*HALF);

  if (failure) {
    printf(">>> FAILED test -- N_VMinQuotient Case 1, Proc %d \n", myid);
    printf("    min = %" FSYM ", expected %" FSYM "\n", ans, HALF*HALF);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMinQuotient Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(NUM, stop_time - start_time);
  PRINT_TIME("N_VMinQuotient", maxt);

  /*
   * Case 2: Fail
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, NUM);
  N_VConst(ZERO, DENOM);

  start_time = get_time();
  ans = N_VMinQuotient(NUM, DENOM);
  sync_device(NUM);
  stop_time = get_time();

  /* ans should equal BIG_REAL */
  failure = SUNRCompare(ans, BIG_REAL);

  if (failure) {
    printf(">>> FAILED test -- N_VMinQuotient Case 2, Proc %d \n", myid);
    printf("    min = %" FSYM ", expected %" FSYM "\n", ans, BIG_REAL);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMinQuotient Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(NUM, stop_time - start_time);
  PRINT_TIME("N_VMinQuotient", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VLinearCombination Test
 * --------------------------------------------------------------------*/
int Test_N_VLinearCombination(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  N_Vector Y1, Y2, Y3;
  N_Vector V[3];
  realtype c[3];

  /* create vectors for testing */
  Y1 = N_VClone(X);
  Y2 = N_VClone(X);
  Y3 = N_VClone(X);

  /* set vectors in vector array */
  V[0] = Y1;
  V[1] = Y2;
  V[2] = Y3;

  /* initialize c values */
  c[0] = ZERO;
  c[1] = ZERO;
  c[2] = ZERO;

  /*
   * Case 1a: V[0] = a V[0], N_VScale
   */

  /* fill vector data */
  N_VConst(TWO, Y1);

  /* set scaling factors */
  c[0] = HALF;

  start_time = get_time();
  ierr = N_VLinearCombination(1, c, V, Y1);
  sync_device(X);
  stop_time = get_time();

  /* Y1 should be vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, Y1, local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombination Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombination Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearCombination", maxt);

  /*
   * Case 1b: X = a V[0], N_VScale
   */

  /* fill vector data and scaling factors */
  N_VConst(TWO, Y1);

  c[0] = HALF;

  N_VConst(ZERO, X);

  start_time = get_time();
  ierr = N_VLinearCombination(1, c, V, X);
  sync_device(X);
  stop_time = get_time();

  /* X should be vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, X, local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombination Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombination Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearCombination", maxt);

  /*
   * Case 2a: V[0] = a V[0] + b V[1], N_VLinearSum
   */

  /* fill vector data */
  N_VConst(NEG_TWO, Y1);
  N_VConst(ONE,     Y2);

  /* set scaling factors */
  c[0] = HALF;
  c[1] = TWO;

  start_time = get_time();
  ierr = N_VLinearCombination(2, c, V, Y1);
  sync_device(X);
  stop_time = get_time();

  /* Y1 should be vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, Y1, local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombination Case 2a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombination Case 2a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearCombination", maxt);

  /*
   * Case 2b: X = a V[0] + b V[1], N_VLinearSum
   */

  /* fill vector data and scaling factors */
  N_VConst(ONE,     Y1);
  N_VConst(NEG_TWO, Y2);

  c[0] = TWO;
  c[1] = HALF;

  N_VConst(ZERO, X);

  start_time = get_time();
  ierr = N_VLinearCombination(2, c, V, X);
  sync_device(X);
  stop_time = get_time();

  /* X should be vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, X, local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombination Case 2b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombination Case 2b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearCombination", maxt);


  /*
   * Case 3a: V[0] = V[0] + b V[1] + c V[2]
   */

  /* fill vector data */
  N_VConst(TWO,     Y1);
  N_VConst(NEG_TWO, Y2);
  N_VConst(NEG_ONE, Y3);

  /* set scaling factors */
  c[0] = ONE;
  c[1] = HALF;
  c[2] = NEG_TWO;

  start_time = get_time();
  ierr = N_VLinearCombination(3, c, V, Y1);
  sync_device(X);
  stop_time = get_time();

  /* Y1 should be vector of +3 */
  if (ierr == 0)
    failure = check_ans(TWO+ONE, Y1, local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombination Case 3a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombination Case 3a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearCombination", maxt);

  /*
   * Case 3b: V[0] = a V[0] + b V[1] + c V[2]
   */

  /* fill vector data */
  N_VConst(ONE,     Y1);
  N_VConst(NEG_TWO, Y2);
  N_VConst(NEG_ONE, Y3);

  /* set scaling factors */
  c[0] = TWO;
  c[1] = HALF;
  c[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VLinearCombination(3, c, V, Y1);
  sync_device(X);
  stop_time = get_time();

  /* Y1 should be vector of +2 */
  if (ierr == 0)
    failure = check_ans(TWO, Y1, local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombination Case 3b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombination Case 3b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearCombination", maxt);

  /*
   * Case 3c: X = a V[0] + b V[1] + c V[2]
   */

  /* fill vector data and set scaling factors */
  N_VConst(ONE,     Y1);
  N_VConst(NEG_TWO, Y2);
  N_VConst(NEG_ONE, Y3);

  c[0] = TWO;
  c[1] = HALF;
  c[2] = NEG_ONE;

  N_VConst(ZERO, X);

  start_time = get_time();
  ierr = N_VLinearCombination(3, c, V, X);
  sync_device(X);
  stop_time = get_time();

  /* X should be vector of +2 */
  if (ierr == 0)
    failure = check_ans(TWO, X, local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombination Case 3c, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombination Case 3c \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VLinearCombination", maxt);

  /* Free vectors */
  N_VDestroy(Y1);
  N_VDestroy(Y2);
  N_VDestroy(Y3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VScaleAddMulti Test
 * --------------------------------------------------------------------*/
int Test_N_VScaleAddMulti(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  realtype avals[3];
  N_Vector *V, *Z;

  /* create vectors for testing */
  Z = N_VCloneVectorArray(3, X);
  V = N_VCloneVectorArray(3, X);

  /* initialize a values */
  avals[0] = ZERO;
  avals[1] = ZERO;
  avals[2] = ZERO;

  /*
   * Case 1a: V[0] = a[0] x + V[0], N_VLinearSum
   */

  /* fill vector data */
  N_VConst(ONE,     X);
  N_VConst(NEG_ONE, V[0]);

  /* set scaling factors */
  avals[0] = TWO;

  start_time = get_time();
  ierr = N_VScaleAddMulti(1, avals, X, V, V);
  sync_device(X);
  stop_time = get_time();

  /* V[0] should be vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, V[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMulti Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMulti Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMulti", maxt);

  /*
   * Case 1b: Z[0] = a[0] x + V[0], N_VLinearSum
   */

  /* fill vector data and set scaling factors */
  N_VConst(ONE,     X);
  N_VConst(NEG_ONE, V[0]);

  avals[0] = TWO;

  N_VConst(ZERO, Z[0]);

  start_time = get_time();
  ierr = N_VScaleAddMulti(1, avals, X, V, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z[0] should be vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, Z[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMulti Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMulti Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMulti", maxt);

  /*
   * Case 2a: V[i] = a[i] x + V[i], N_VScaleAddMulti
   */

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(NEG_TWO, V[0]);
  N_VConst(TWO,     V[1]);
  N_VConst(NEG_ONE, V[2]);

  /* set scaling factors */
  avals[0] = ONE;
  avals[1] = NEG_TWO;
  avals[2] = TWO;

  start_time = get_time();
  ierr = N_VScaleAddMulti(3, avals, X, V, V);
  sync_device(X);
  stop_time = get_time();

  /* V[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, V[0], local_length);
    failure += check_ans(ZERO,    V[1], local_length);
    failure += check_ans(ONE,     V[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMulti Case 2a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMulti Case 2a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMulti", maxt);

  /*
   * Case 2b: Z[i] = a[i] x + V[i], N_VScaleAddMulti
   */

  /* fill vector data and set scaling factors */
  N_VConst(ONE, X);
  N_VConst(NEG_TWO, V[0]);
  N_VConst(TWO,     V[1]);
  N_VConst(NEG_ONE, V[2]);

  avals[0] = ONE;
  avals[1] = NEG_TWO;
  avals[2] = TWO;

  N_VConst(TWO, Z[0]);
  N_VConst(TWO, Z[1]);
  N_VConst(TWO, Z[2]);

  start_time = get_time();
  ierr = N_VScaleAddMulti(3, avals, X, V, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(ZERO,    Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  }
  else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMulti Case 2b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMulti Case 2b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMulti", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(Z, 3);
  N_VDestroyVectorArray(V, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VDotProdMulti Test
 * --------------------------------------------------------------------*/
int Test_N_VDotProdMulti(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  sunindextype  global_length;
  N_Vector     *V;
  realtype      dotprods[3];

  /* get global length */
  global_length = N_VGetLength(X);

  /* create vectors for testing */
  V = N_VCloneVectorArray(3, X);

  /*
   * Case 1: d[0] = z . V[0], N_VDotProd
   */

  /* fill vector data */
  N_VConst(TWO,  X);
  N_VConst(HALF, V[0]);

  start_time = get_time();
  ierr = N_VDotProdMulti(1, X, V, dotprods);
  sync_device(X);
  stop_time = get_time();

  /* dotprod[0] should equal the global vector length */
  if (ierr == 0)
    failure = SUNRCompare(dotprods[0], (realtype) global_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdMulti Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdMulti Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdMulti", maxt);

  /*
   * Case 2: d[i] = z . V[i], N_VDotProd
   */

  /* fill vector data */
  N_VConst(TWO,      X);
  N_VConst(NEG_HALF, V[0]);
  N_VConst(HALF,     V[1]);
  N_VConst(ONE,      V[2]);

  start_time = get_time();
  ierr = N_VDotProdMulti(3, X, V, dotprods);
  sync_device(X);
  stop_time = get_time();

  /* dotprod[i] should equal -1, +1, and 2 times the global vector length */
  if (ierr == 0) {
    failure  = SUNRCompare(dotprods[0], (realtype) -1*global_length);
    failure += SUNRCompare(dotprods[1], (realtype)    global_length);
    failure += SUNRCompare(dotprods[2], (realtype)  2*global_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdMulti Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdMulti Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdMulti", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(V, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VLinearSumVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VLinearSumVectorArray(N_Vector V, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  N_Vector *X, *Y, *Z;

  /* create vectors for testing */
  X = N_VCloneVectorArray(3, V);
  Y = N_VCloneVectorArray(3, V);
  Z = N_VCloneVectorArray(3, V);

  /*
   * Case 0: Z[0] = a X[0] + b Y[0], N_VLinearSum
   */

  /* fill vector data */
  N_VConst(NEG_HALF, X[0]);
  N_VConst(TWO,      Y[0]);
  N_VConst(TWO,      Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(1, TWO, X, HALF, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[0] should be a vector of 0 */
  if (ierr == 0)
    failure = check_ans(ZERO, Z[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 0, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 0 \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 1a: Y[i] = X[i] + Y[i], (VaxpyVectorArray Case 1)
   */

  /* fill vector data */
  N_VConst(NEG_TWO, X[0]);
  N_VConst(ONE,     Y[0]);

  N_VConst(TWO,     X[1]);
  N_VConst(NEG_TWO, Y[1]);

  N_VConst(TWO,     X[2]);
  N_VConst(NEG_ONE, Y[2]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, ONE, X, ONE, Y, Y);
  sync_device(V);
  stop_time = get_time();

  /* Y[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Y[0], local_length);
    failure += check_ans(ZERO,    Y[1], local_length);
    failure += check_ans(ONE,     Y[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 1b: y = -x + y, (VaxpyVectorArray Case 2)
   */

  /* fill vector data */
  N_VConst(TWO, X[0]);
  N_VConst(ONE, Y[0]);

  N_VConst(NEG_TWO, X[1]);
  N_VConst(NEG_TWO, Y[1]);

  N_VConst(NEG_TWO, X[2]);
  N_VConst(NEG_ONE, Y[2]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, NEG_ONE, X, ONE, Y, Y);
  sync_device(V);
  stop_time = get_time();

  /* Y[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Y[0], local_length);
    failure += check_ans(ZERO,    Y[1], local_length);
    failure += check_ans(ONE,     Y[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 1c: y = ax + y, (VaxpyVectorArray Case 3)
   */

  /* fill vector data */
  N_VConst(TWO,     X[0]);
  N_VConst(NEG_TWO, Y[0]);

  N_VConst(TWO,     X[1]);
  N_VConst(NEG_ONE, Y[1]);

  N_VConst(NEG_TWO, X[2]);
  N_VConst(TWO,     Y[2]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, HALF, X, ONE, Y, Y);
  sync_device(V);
  stop_time = get_time();

  /* Y[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Y[0], local_length);
    failure += check_ans(ZERO,    Y[1], local_length);
    failure += check_ans(ONE,     Y[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 1c, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 1c \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 2a: x = x + y, (VaxpyVectorArray Case 1)
   */

  /* fill vector data */
  N_VConst(NEG_TWO, X[0]);
  N_VConst(ONE,     Y[0]);

  N_VConst(TWO,     X[1]);
  N_VConst(NEG_TWO, Y[1]);

  N_VConst(TWO,     X[2]);
  N_VConst(NEG_ONE, Y[2]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, ONE, X, ONE, Y, X);
  sync_device(V);
  stop_time = get_time();

  /* X[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, X[0], local_length);
    failure += check_ans(ZERO,    X[1], local_length);
    failure += check_ans(ONE,     X[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 2a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 2a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 2b: x = x - y, (VaxpyVectorArray Case 2)
   */

  /* fill vector data */
  N_VConst(ONE, X[0]);
  N_VConst(TWO, Y[0]);

  N_VConst(NEG_TWO, X[1]);
  N_VConst(NEG_TWO, Y[1]);

  N_VConst(NEG_ONE, X[2]);
  N_VConst(NEG_TWO, Y[2]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, ONE, X, NEG_ONE, Y, X);
  sync_device(V);
  stop_time = get_time();

  /* X[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, X[0], local_length);
    failure += check_ans(ZERO,    X[1], local_length);
    failure += check_ans(ONE,     X[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 2b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 2b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 2c: x = x + by, (VaxpyVectorArray Case 3)
   */

  /* fill vector data */
  N_VConst(NEG_TWO, X[0]);
  N_VConst(TWO,     Y[0]);

  N_VConst(NEG_ONE, X[1]);
  N_VConst(TWO,     Y[1]);

  N_VConst(TWO,     X[2]);
  N_VConst(NEG_TWO, Y[2]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, ONE, X, HALF, Y, X);
  sync_device(V);
  stop_time = get_time();

  /* X[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, X[0], local_length);
    failure += check_ans(ZERO,    X[1], local_length);
    failure += check_ans(ONE,     X[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 2c, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 2c \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 3: z = x + y, (VSumVectorArray)
   */

  /* fill vector data */
  N_VConst(NEG_TWO, X[0]);
  N_VConst(ONE,     Y[0]);
  N_VConst(TWO,     Z[0]);

  N_VConst(NEG_ONE, X[1]);
  N_VConst(ONE,     Y[1]);
  N_VConst(TWO,     Z[0]);

  N_VConst(TWO,     X[2]);
  N_VConst(NEG_ONE, Y[2]);
  N_VConst(TWO,     Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, ONE, X, ONE, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(ZERO,    Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 3, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 3 \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 4a: z = x - y, (VDiffVectorArray)
   */

  /* fill vector data */
  N_VConst(NEG_TWO, X[0]);
  N_VConst(NEG_ONE, Y[0]);
  N_VConst(TWO,     Z[0]);

  N_VConst(NEG_ONE, X[1]);
  N_VConst(NEG_ONE, Y[1]);
  N_VConst(TWO,     Z[0]);

  N_VConst(TWO, X[2]);
  N_VConst(ONE, Y[2]);
  N_VConst(TWO, Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, ONE, X, NEG_ONE, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(ZERO,    Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 4a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 4a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 4b: z = -x + y, (VDiffVectorArray)
   */

  /* fill vector data */
  N_VConst(TWO, X[0]);
  N_VConst(ONE, Y[0]);
  N_VConst(TWO, Z[0]);

  N_VConst(NEG_ONE, X[1]);
  N_VConst(NEG_ONE, Y[1]);
  N_VConst(TWO,     Z[0]);

  N_VConst(NEG_TWO, X[2]);
  N_VConst(NEG_ONE, Y[2]);
  N_VConst(TWO,     Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, NEG_ONE, X, ONE, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(ZERO,    Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 4b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 4b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 5a: z = x + by, (VLin1VectorArray)
   */

  /* fill vector data */
  N_VConst(NEG_TWO, X[0]);
  N_VConst(TWO,     Y[0]);
  N_VConst(TWO,     Z[0]);

  N_VConst(ONE,     X[1]);
  N_VConst(NEG_TWO, Y[1]);
  N_VConst(TWO,     Z[0]);

  N_VConst(HALF, X[2]);
  N_VConst(ONE,  Y[2]);
  N_VConst(TWO,  Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, ONE, X, HALF, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(ZERO,    Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 5a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 5a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 5b: z = ax + y, (VLin1VectorArray)
   */

  /* fill vector data */
  N_VConst(NEG_TWO, X[0]);
  N_VConst(NEG_TWO, Y[0]);
  N_VConst(TWO,     Z[0]);

  N_VConst(ONE,  X[1]);
  N_VConst(HALF, Y[1]);
  N_VConst(TWO,  Z[0]);

  N_VConst(TWO, X[2]);
  N_VConst(TWO, Y[2]);
  N_VConst(TWO, Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, NEG_HALF, X, ONE, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(ZERO,    Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 5b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 5b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 6a: z = -x + by, (VLin2VectorArray)
   */

  /* fill vector data */
  N_VConst(HALF,    X[0]);
  N_VConst(NEG_ONE, Y[0]);
  N_VConst(TWO,     Z[0]);

  N_VConst(ONE, X[1]);
  N_VConst(TWO, Y[1]);
  N_VConst(TWO, Z[0]);

  N_VConst(NEG_TWO, X[2]);
  N_VConst(NEG_TWO, Y[2]);
  N_VConst(TWO,     Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, NEG_ONE, X, HALF, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(ZERO,    Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 6a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 6a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 6b: z = ax - y, (VLin2VectorArray)
   */

  /* fill vector data */
  N_VConst(HALF, X[0]);
  N_VConst(TWO,  Y[0]);
  N_VConst(TWO,  Z[0]);

  N_VConst(ONE, X[1]);
  N_VConst(TWO, Y[1]);
  N_VConst(TWO, Z[0]);

  N_VConst(NEG_HALF, X[2]);
  N_VConst(NEG_TWO,  Y[2]);
  N_VConst(TWO,      Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, TWO, X, NEG_ONE, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(ZERO,    Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 6b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 6b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 7: z = a(x + y), (VScaleSumVectorArray)
   */

  /* fill vector data */
  N_VConst(NEG_ONE, X[0]);
  N_VConst(HALF,    Y[0]);
  N_VConst(TWO,     Z[0]);

  N_VConst(ONE,  X[1]);
  N_VConst(HALF, Y[1]);
  N_VConst(TWO,  Z[0]);

  N_VConst(ONE,      X[2]);
  N_VConst(NEG_HALF, Y[2]);
  N_VConst(TWO,      Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, TWO, X, TWO, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 3, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(TWO+ONE, Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 7, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 7 \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 8: z = a(x - y), (VScaleDiffVectorArray)
   */

  /* fill vector data */
  N_VConst(HALF, X[0]);
  N_VConst(ONE,  Y[0]);
  N_VConst(TWO,  Z[0]);

  N_VConst(TWO,  X[1]);
  N_VConst(HALF, Y[1]);
  N_VConst(TWO,  Z[0]);

  N_VConst(NEG_HALF, X[2]);
  N_VConst(NEG_ONE,  Y[2]);
  N_VConst(TWO,      Z[0]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, TWO, X, NEG_TWO, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of -1, 3, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0], local_length);
    failure += check_ans(TWO+ONE, Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 8, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 8 \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /*
   * Case 9: z = ax + by, All Other Cases
   */

  /* fill vector data */
  N_VConst(NEG_HALF, X[0]);
  N_VConst(TWO,      Y[0]);

  N_VConst(ONE,     X[1]);
  N_VConst(NEG_TWO, Y[1]);

  N_VConst(HALF, X[2]);
  N_VConst(TWO,  Y[2]);

  start_time = get_time();
  ierr = N_VLinearSumVectorArray(3, TWO, X, HALF, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should be a vector of 0, +1, +2 */
  if (ierr == 0) {
    failure  = check_ans(ZERO, Z[0], local_length);
    failure += check_ans(ONE,  Z[1], local_length);
    failure += check_ans(TWO,  Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSumVectorArray Case 9, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearSumVectorArray Case 9 \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearSumVectorArray", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(X, 3);
  N_VDestroyVectorArray(Y, 3);
  N_VDestroyVectorArray(Z, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VScaleVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VScaleVectorArray(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  realtype c[3];
  N_Vector *Y, *Z;

  /* create vectors for testing */
  Y = N_VCloneVectorArray(3, X);
  Z = N_VCloneVectorArray(3, X);

  /*
   * Case 1a: Y[0] = c[0] Y[0], N_VScale
   */

  /* fill vector data */
  N_VConst(HALF, Y[0]);

  c[0] = TWO;

  start_time = get_time();
  ierr = N_VScaleVectorArray(1, c, Y, Y);
  sync_device(X);
  stop_time = get_time();

  /* Y[0] should be a vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, Y[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VScaleVectorArray Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleVectorArray Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScaleVectorArray", maxt);

  /*
   * Case 1b: Z[0] = c[0] Y[0], N_VScale
   */

  /* fill vector data */
  N_VConst(HALF, Y[0]);

  c[0] = TWO;

  start_time = get_time();
  ierr = N_VScaleVectorArray(1, c, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z[0] should be a vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, Z[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VScaleVectorArray Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleVectorArray Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScaleVectorArray", maxt);

  /*
   * Case 2a: Y[i] = c[i] Y[i]
   */

  /* fill vector data */
  N_VConst(HALF,    Y[0]);
  N_VConst(NEG_TWO, Y[1]);
  N_VConst(NEG_ONE, Y[2]);

  c[0] = TWO;
  c[1] = HALF;
  c[2] = NEG_TWO;

  start_time = get_time();
  ierr = N_VScaleVectorArray(3, c, Y, Y);
  sync_device(X);
  stop_time = get_time();

  /* Y[i] should be a vector of +1, -1, 2 */
  if (ierr == 0) {
    failure  = check_ans(ONE,     Y[0], local_length);
    failure += check_ans(NEG_ONE, Y[1], local_length);
    failure += check_ans(TWO,     Y[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleVectorArray Case 2a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleVectorArray Case 2a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScaleVectorArray", maxt);

  /*
   * Case 2b: Z[i] = c[i] Y[i]
   */

  /* fill vector data */
  N_VConst(HALF,    Y[0]);
  N_VConst(NEG_TWO, Y[1]);
  N_VConst(NEG_ONE, Y[2]);

  c[0] = TWO;
  c[1] = HALF;
  c[2] = NEG_TWO;

  start_time = get_time();
  ierr = N_VScaleVectorArray(3, c, Y, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z[i] should be a vector of +1, -1, 2 */
  if (ierr == 0) {
    failure  = check_ans(ONE,     Z[0], local_length);
    failure += check_ans(NEG_ONE, Z[1], local_length);
    failure += check_ans(TWO,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleVectorArray Case 2b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleVectorArray Case 2b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VScaleVectorArray", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(Y, 3);
  N_VDestroyVectorArray(Z, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VConstVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VConstVectorArray(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  N_Vector *Z;

  /* create vectors for testing */
  Z = N_VCloneVectorArray(3, X);

  /*
   * Case 1a: Z[0] = c, N_VConst
   */

  /* fill vector data */
  N_VConst(ZERO, Z[0]);

  start_time = get_time();
  ierr = N_VConstVectorArray(1, ONE, Z);
  sync_device(X);
  stop_time = get_time();

  /* Y[0] should be a vector of 1 */
  if (ierr == 0)
    failure = check_ans(ONE, Z[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VConstVectorArray Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VConstVectorArray Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VConstVectorArray", maxt);

  /*
   * Case 1b: Z[i] = c
   */

  /* fill vector data */
  N_VConst(ZERO, Z[0]);
  N_VConst(ZERO, Z[1]);
  N_VConst(ZERO, Z[2]);

  start_time = get_time();
  ierr = N_VConstVectorArray(3, ONE, Z);
  sync_device(X);
  stop_time = get_time();

  /* Y[i] should be a vector of 1 */
  if (ierr == 0) {
    failure  = check_ans(ONE, Z[0], local_length);
    failure += check_ans(ONE, Z[1], local_length);
    failure += check_ans(ONE, Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VConstVectorArray Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VConstVectorArray Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VConstVectorArray", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(Z, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNormVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNormVectorArray(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  realtype nrm[3];
  N_Vector *Z;
  N_Vector *W;

  /* create vectors for testing */
  Z = N_VCloneVectorArray(3, X);
  W = N_VCloneVectorArray(3, X);

  /*
   * Case 1a: nrm[0] = ||Z[0]||, N_VWrmsNorm
   */

  /* fill vector data */
  N_VConst(NEG_HALF, Z[0]);
  N_VConst(HALF,     W[0]);

  nrm[0] = NEG_ONE;
  nrm[1] = NEG_ONE;
  nrm[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VWrmsNormVectorArray(1, Z, W, nrm);
  sync_device(X);
  stop_time = get_time();

  /* nrm should equal 1/4 */
  if (ierr == 0)
    failure = (nrm[0] < ZERO) ? 1 : SUNRCompare(nrm[0], HALF*HALF);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNormVectorArray Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWrmsNormVectorArray Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWrmsNormVectorArray", maxt);

  /*
   * Case 1b: nrm[i] = ||Z[i]||
   */

  /* fill vector data */
  N_VConst(NEG_HALF, Z[0]);
  N_VConst(TWO*TWO,  Z[1]);
  N_VConst(HALF,     Z[2]);

  N_VConst(HALF,      W[0]);
  N_VConst(HALF*HALF, W[1]);
  N_VConst(ONE,       W[2]);

  nrm[0] = NEG_ONE;
  nrm[1] = NEG_ONE;
  nrm[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VWrmsNormVectorArray(3, Z, W, nrm);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal 1/4, 1, 1/2 */
  if (ierr == 0) {
    failure  = (nrm[0] < ZERO) ? 1 : SUNRCompare(nrm[0], HALF*HALF);
    failure += (nrm[1] < ZERO) ? 1 : SUNRCompare(nrm[1], ONE);
    failure += (nrm[2] < ZERO) ? 1 : SUNRCompare(nrm[2], HALF);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNormVectorArray Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWrmsNormVectorArray Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWrmsNormVectorArray", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(Z, 3);
  N_VDestroyVectorArray(W, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNormMaskVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNormMaskVectorArray(N_Vector X, sunindextype local_length,
                                    int myid)
{
  int          fails = 0, failure = 0, ierr = 0;
  double       start_time, stop_time, maxt;
  sunindextype global_length;

  realtype fac;
  realtype nrm[3];
  N_Vector *Z;
  N_Vector *W;

  /* get global length */
  global_length = N_VGetLength(X);

  /* factor used in checking solutions */
  fac = SUNRsqrt((realtype) (global_length - 1)/(global_length));

  /* create vectors for testing */
  Z = N_VCloneVectorArray(3, X);
  W = N_VCloneVectorArray(3, X);

  /*
   * Case 1: nrm[0] = ||Z[0]||
   */

  /* fill vector data */
  N_VConst(NEG_HALF, Z[0]);
  N_VConst(HALF,     W[0]);

  /* use all elements except one */
  N_VConst(ONE, X);
  if (myid == 0)
    set_element(X, local_length-1, ZERO);

  nrm[0] = NEG_ONE;
  nrm[1] = NEG_ONE;
  nrm[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VWrmsNormMaskVectorArray(1, Z, W, X, nrm);
  sync_device(X);
  stop_time = get_time();

  /* nrm should equal fac/4 */
  if (ierr == 0)
    failure = (nrm[0] < ZERO) ? 1 : SUNRCompare(nrm[0], fac*HALF*HALF);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNormMaskVectorArray Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWrmsNormMaskVectorArray Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWrmsNormVectorArray", maxt);

  /*
   * Case 2: nrm[i] = ||Z[i]||
   */

  /* fill vector data */
  N_VConst(NEG_HALF, Z[0]);
  N_VConst(TWO*TWO,  Z[1]);
  N_VConst(HALF,     Z[2]);

  N_VConst(HALF,      W[0]);
  N_VConst(HALF*HALF, W[1]);
  N_VConst(ONE,       W[2]);

  /* use all elements except one */
  N_VConst(ONE, X);
  if (myid == 0)
    set_element(X, local_length-1, ZERO);

  nrm[0] = NEG_ONE;
  nrm[1] = NEG_ONE;
  nrm[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VWrmsNormMaskVectorArray(3, Z, W, X, nrm);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal fac/4, fac, fac/2] */
  if (ierr == 0) {
    failure  = (nrm[0] < ZERO) ? 1 : SUNRCompare(nrm[0], fac*HALF*HALF);
    failure += (nrm[1] < ZERO) ? 1 : SUNRCompare(nrm[1], fac);
    failure += (nrm[2] < ZERO) ? 1 : SUNRCompare(nrm[2], fac*HALF);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNormMaskVectorArray Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWrmsNormMaskVectorArray Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWrmsNormVectorArray", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(Z, 3);
  N_VDestroyVectorArray(W, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VScaleAddMultiVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VScaleAddMultiVectorArray(N_Vector V, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  realtype  a[3];
  N_Vector* X;
  N_Vector* Y[3];
  N_Vector* Z[3];


  /* create vectors for testing */
  X = N_VCloneVectorArray(3, V);

  Y[0] = N_VCloneVectorArray(3, V);
  Y[1] = N_VCloneVectorArray(3, V);
  Y[2] = N_VCloneVectorArray(3, V);

  Z[0] = N_VCloneVectorArray(3, V);
  Z[1] = N_VCloneVectorArray(3, V);
  Z[2] = N_VCloneVectorArray(3, V);

  /*
   * Case 1a (nvec = 1, nsum = 1):
   * Z[0][0] = a[0] X[0] + Y[0][0], N_VLinearSum
   */

  /* fill scaling and vector data */
  a[0] = TWO;

  N_VConst(ONE,     X[0]);
  N_VConst(NEG_ONE, Y[0][0]);

  start_time = get_time();
  ierr = N_VScaleAddMultiVectorArray(1, 1, a, X, Y, Y);
  sync_device(V);
  stop_time = get_time();

  /* Y[0][0] should be vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, Y[0][0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMultiVectorArray Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMultiVectorArray Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMultiVectorArray", maxt);

  /*
   * Case 1b (nvec = 1, nsum = 1):
   * Z[0][0] = a[0] X[0] + Y[0][0], N_VLinearSum
   */

  /* fill scaling and vector data */
  a[0] = TWO;

  N_VConst(ONE,     X[0]);
  N_VConst(NEG_ONE, Y[0][0]);
  N_VConst(ZERO,    Z[0][0]);

  start_time = get_time();
  ierr = N_VScaleAddMultiVectorArray(1, 1, a, X, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[0][0] should be vector of +1 */
  if (ierr == 0)
    failure = check_ans(ONE, Z[0][0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMultiVectorArray Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMultiVectorArray Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMultiVectorArray", maxt);

  /*
   * Case 2a (nvec = 1, nsum > 1):
   * Y[j][0] = a[j] X[0] + Y[j][0], N_VScaleAddMulti
   */

  /* fill scaling and vector data */
  a[0] = ONE;
  a[1] = NEG_TWO;
  a[2] = TWO;

  N_VConst(ONE, X[0]);

  N_VConst(NEG_TWO, Y[0][0]);
  N_VConst(TWO,     Y[1][0]);
  N_VConst(NEG_ONE, Y[2][0]);

  start_time = get_time();
  ierr = N_VScaleAddMultiVectorArray(1, 3, a, X, Y, Y);
  sync_device(V);
  stop_time = get_time();

  /* Y[i][0] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Y[0][0], local_length);
    failure += check_ans(ZERO,    Y[1][0], local_length);
    failure += check_ans(ONE,     Y[2][0], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMultiVectorArray Case 2a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMultiVectorArray Case 2a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMultiVectorArray", maxt);

  /*
   * Case 2b (nvec = 1, nsum > 1):
   * Z[j][0] = a[j] X[0] + Y[j][0], N_VScaleAddMulti
   */

  /* fill scaling and vector data */
  a[0] = ONE;
  a[1] = NEG_TWO;
  a[2] = TWO;

  N_VConst(ONE, X[0]);

  N_VConst(NEG_TWO, Y[0][0]);
  N_VConst(TWO,     Y[1][0]);
  N_VConst(NEG_ONE, Y[2][0]);

  N_VConst(ZERO, Z[0][0]);
  N_VConst(ONE,  Z[1][0]);
  N_VConst(TWO,  Z[2][0]);

  start_time = get_time();
  ierr = N_VScaleAddMultiVectorArray(1, 3, a, X, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i][0] should be a vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0][0], local_length);
    failure += check_ans(ZERO,    Z[1][0], local_length);
    failure += check_ans(ONE,     Z[2][0], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMultiVectorArray Case 2b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMultiVectorArray Case 2b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMultiVectorArray", maxt);

  /*
   * Case 3a (nvec > 1, nsum = 1):
   * Y[0][i] = a[0] X[i] + Y[0][i], N_VLinearSumVectorArray
   */

  /* fill scaling and vector data */
  a[0] = TWO;

  N_VConst(HALF,    X[0]);
  N_VConst(NEG_ONE, X[1]);
  N_VConst(ONE,     X[2]);

  N_VConst(NEG_TWO, Y[0][0]);
  N_VConst(TWO,     Y[0][1]);
  N_VConst(NEG_ONE, Y[0][2]);

  start_time = get_time();
  ierr = N_VScaleAddMultiVectorArray(3, 1, a, X, Y, Y);
  sync_device(V);
  stop_time = get_time();

  /* Y[0][i] should be vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Y[0][0], local_length);
    failure += check_ans(ZERO,    Y[0][1], local_length);
    failure += check_ans(ONE,     Y[0][2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMultiVectorArray Case 3a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMultiVectorArray Case 3a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMultiVectorArray", maxt);

  /*
   * Case 3b (nvec > 1, nsum = 1):
   * Z[j][0] = a[j] X[0] + Y[j][0], N_VLinearSumVectorArray
   */

  /* fill scaling and vector data */
  a[0] = TWO;

  N_VConst(HALF,    X[0]);
  N_VConst(NEG_ONE, X[1]);
  N_VConst(ONE,     X[2]);

  N_VConst(NEG_TWO, Y[0][0]);
  N_VConst(TWO,     Y[0][1]);
  N_VConst(NEG_ONE, Y[0][2]);

  N_VConst(TWO, Z[0][0]);
  N_VConst(TWO, Z[0][1]);
  N_VConst(TWO, Z[0][2]);

  start_time = get_time();
  ierr = N_VScaleAddMultiVectorArray(3, 1, a, X, Y, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[0][i] should be vector of -1, 0, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE, Z[0][0], local_length);
    failure += check_ans(ZERO,    Z[0][1], local_length);
    failure += check_ans(ONE,     Z[0][2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMultiVectorArray Case 3b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMultiVectorArray Case 3b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMultiVectorArray", maxt);

  /*
   * Case 4a (nvec > 1, nsum > 1):
   * Y[j][i] = a[j] X[i] + Y[j][i], N_VScaleAddMultiVectorArray
   */

  /* fill scaling and vector data */
  a[0] = TWO;
  a[1] = ONE;
  a[2] = NEG_TWO;

  N_VConst(HALF,     X[0]);
  N_VConst(NEG_TWO,  Y[0][0]);
  N_VConst(NEG_HALF, Y[1][0]);
  N_VConst(TWO,      Y[2][0]);

  N_VConst(ONE,     X[1]);
  N_VConst(NEG_ONE, Y[0][1]);
  N_VConst(NEG_TWO, Y[1][1]);
  N_VConst(TWO,     Y[2][1]);

  N_VConst(NEG_TWO,     X[2]);
  N_VConst(TWO,         Y[0][2]);
  N_VConst(TWO*TWO,     Y[1][2]);
  N_VConst(NEG_TWO*TWO, Y[2][2]);

  start_time = get_time();
  ierr = N_VScaleAddMultiVectorArray(3, 3, a, X, Y, Y);
  sync_device(V);
  stop_time = get_time();

  if (ierr == 0) {
    /* Y[i][0] should be vector of -1, 0, +1 */
    failure  = check_ans(NEG_ONE, Y[0][0], local_length);
    failure += check_ans(ZERO,    Y[1][0], local_length);
    failure += check_ans(ONE,     Y[2][0], local_length);

    /* Y[i][1] should be vector of +1, -1, 0 */
    failure += check_ans(ONE,     Y[0][1], local_length);
    failure += check_ans(NEG_ONE, Y[1][1], local_length);
    failure += check_ans(ZERO,    Y[2][1], local_length);

    /* Y[i][2] should be vector of -2, 2, 0 */
    failure += check_ans(NEG_TWO, Y[0][2], local_length);
    failure += check_ans(TWO,     Y[1][2], local_length);
    failure += check_ans(ZERO,    Y[2][2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMultiVectorArray Case 4a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMultiVectorArray Case 4a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMultiVectorArray", maxt);

  /*
   * Case 4b (nvec > 1, nsum > 1):
   * Z[j][i] = a[j] X[i] + Y[j][i], N_VScaleAddMultiVectorArray
   */

  /* fill scaling and vector data */
  a[0] = TWO;
  a[1] = ONE;
  a[2] = NEG_TWO;

  N_VConst(HALF, X[0]);

  N_VConst(NEG_TWO,  Y[0][0]);
  N_VConst(NEG_HALF, Y[1][0]);
  N_VConst(TWO,      Y[2][0]);

  N_VConst(HALF, Z[0][0]);
  N_VConst(HALF, Z[1][0]);
  N_VConst(HALF, Z[2][0]);

  N_VConst(ONE, X[1]);

  N_VConst(NEG_ONE, Y[0][1]);
  N_VConst(NEG_TWO, Y[1][1]);
  N_VConst(TWO,     Y[2][1]);

  N_VConst(HALF, Z[0][1]);
  N_VConst(HALF, Z[1][1]);
  N_VConst(HALF, Z[2][1]);

  N_VConst(NEG_TWO, X[2]);

  N_VConst(TWO,         Y[0][2]);
  N_VConst(TWO*TWO,     Y[1][2]);
  N_VConst(NEG_TWO*TWO, Y[2][2]);

  N_VConst(HALF, Z[0][2]);
  N_VConst(HALF, Z[1][2]);
  N_VConst(HALF, Z[2][2]);

  start_time = get_time();
  ierr = N_VScaleAddMultiVectorArray(3, 3, a, X, Y, Z);
  sync_device(V);
  stop_time = get_time();

  if (ierr == 0) {
    /* Z[i][0] should be vector of -1, 0, +1 */
    failure  = check_ans(NEG_ONE, Z[0][0], local_length);
    failure += check_ans(ZERO,    Z[1][0], local_length);
    failure += check_ans(ONE,     Z[2][0], local_length);

    /* Z[i][1] should be vector of +1, -1, 0 */
    failure += check_ans(ONE,     Z[0][1], local_length);
    failure += check_ans(NEG_ONE, Z[1][1], local_length);
    failure += check_ans(ZERO,    Z[2][1], local_length);

    /* Z[i][2] should be vector of -2, 2, 0 */
    failure += check_ans(NEG_TWO, Z[0][2], local_length);
    failure += check_ans(TWO,     Z[1][2], local_length);
    failure += check_ans(ZERO,    Z[2][2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VScaleAddMultiVectorArray Case 4b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VScaleAddMultiVectorArray Case 4b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VScaleAddMultiVectorArray", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(X, 3);
  N_VDestroyVectorArray(Y[0], 3);
  N_VDestroyVectorArray(Y[1], 3);
  N_VDestroyVectorArray(Y[2], 3);
  N_VDestroyVectorArray(Z[0], 3);
  N_VDestroyVectorArray(Z[1], 3);
  N_VDestroyVectorArray(Z[2], 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VLinearCombinationVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VLinearCombinationVectorArray(N_Vector V, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  realtype c[3];
  N_Vector *Z;
  N_Vector* X[3];

  /* create vectors for testing */
  Z = N_VCloneVectorArray(3, V);

  X[0] = N_VCloneVectorArray(3, V);
  X[1] = N_VCloneVectorArray(3, V);
  X[2] = N_VCloneVectorArray(3, V);

  /*
   * Case 1a: (nvec = 1, nsum = 1), N_VScale
   * X[0][0] = c[0] X[0][0]
   */

  /* fill vector data and scaling factor */
  N_VConst(HALF, X[0][0]);
  c[0] = TWO;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(1, 1, c, X, X[0]);
  sync_device(V);
  stop_time = get_time();

  /* X[0][0] should equal +1 */
  if (ierr == 0)
    failure = check_ans(ONE, X[0][0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 1a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 1a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 1b: (nvec = 1, nsum = 1), N_VScale
   * Z[0] = c[0] X[0][0]
   */

  /* fill vector data and scaling factor */
  N_VConst(HALF, X[0][0]);
  N_VConst(ZERO, Z[0]);
  c[0] = TWO;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(1, 1, c, X, Z);
  sync_device(V);
  stop_time = get_time();

  /* X[0][0] should equal +1 */
  if (ierr == 0)
    failure = check_ans(ONE, Z[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 1b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 1b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 2a: (nvec = 1, nsum = 2), N_VLinearSum
   * X[0][0] = c[0] X[0][0] + c[1] X[1][0]
   */

  /* fill vector data and scaling factor */
  N_VConst(HALF,    X[0][0]);
  N_VConst(NEG_ONE, X[1][0]);

  c[0] = TWO;
  c[1] = NEG_ONE;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(1, 2, c, X, X[0]);
  sync_device(V);
  stop_time = get_time();

  /* X[0][0] should equal +2 */
  if (ierr == 0)
    failure = check_ans(TWO, X[0][0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 2a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 2a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 2b: (nvec = 1, nsum = 2), N_VLinearSum
   * Z[0] = c[0] X[0][0] + c[1] X[1][0]
   */

  /* fill vector data and scaling factor */
  N_VConst(HALF,    X[0][0]);
  N_VConst(NEG_ONE, X[1][0]);

  c[0] = TWO;
  c[1] = NEG_ONE;

  N_VConst(ZERO, Z[0]);

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(1, 2, c, X, Z);
  sync_device(V);
  stop_time = get_time();

  /* X[0][0] should equal +2 */
  if (ierr == 0)
    failure = check_ans(TWO, Z[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 2b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 2b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 3a: (nvec = 1, nsum > 2), N_VLinearCombination
   * X[0][0] = c[0] X[0][0] + c[1] X[1][0] + c[2] X[2][0]
   */

  /* fill vector data */
  N_VConst(ONE,     X[0][0]);
  N_VConst(NEG_TWO, X[1][0]);
  N_VConst(NEG_ONE, X[2][0]);

  /* set scaling factors */
  c[0] = TWO;
  c[1] = HALF;
  c[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(1, 3, c, X, X[0]);
  sync_device(V);
  stop_time = get_time();

  /* X[0][0] should equal +2 */
  if (ierr == 0)
    failure = check_ans(TWO, X[0][0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 3a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 3a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 3b: (nvec = 1, nsum > 2), N_VLinearCombination
   * Z[0] = c[0] X[0][0] + c[1] X[1][0] + c[2] X[2][0]
   */

  /* fill vector data */
  N_VConst(ONE,     X[0][0]);
  N_VConst(NEG_TWO, X[1][0]);
  N_VConst(NEG_ONE, X[2][0]);

  /* set scaling factors */
  c[0] = TWO;
  c[1] = HALF;
  c[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(1, 3, c, X, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[0] should equal +2 */
  if (ierr == 0)
    failure = check_ans(TWO, Z[0], local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 3b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 3b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 4a: (nvec > 1, nsum = 1), N_VScaleVectorArray
   * X[0][i] = c[0] X[0][i]
   */

  /* fill vector data and set scaling factors */
  N_VConst(NEG_TWO, X[0][0]);
  N_VConst(NEG_ONE, X[0][1]);
  N_VConst(TWO,     X[0][2]);

  c[0] = HALF;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(3, 1, c, X, X[0]);
  sync_device(V);
  stop_time = get_time();

  /* X[0][i] should equal to -1, -1/2, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE,  X[0][0], local_length);
    failure += check_ans(NEG_HALF, X[0][1], local_length);
    failure += check_ans(ONE,      X[0][2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 4a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 4a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 4b: (nvec > 1, nsum = 1), N_VScaleVectorArray
   * Z[i] = c[0] X[0][i]
   */

  /* fill vector data and set scaling factors */
  N_VConst(NEG_TWO, X[0][0]);
  N_VConst(NEG_ONE, X[0][1]);
  N_VConst(TWO,     X[0][2]);

  c[0] = HALF;

  N_VConst(ZERO, Z[0]);
  N_VConst(ZERO, Z[1]);
  N_VConst(ZERO, Z[2]);

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(3, 1, c, X, Z);
  sync_device(V);
  stop_time = get_time();

  /* X[0][i] should equal to -1, -1/2, +1 */
  if (ierr == 0) {
    failure  = check_ans(NEG_ONE,  Z[0], local_length);
    failure += check_ans(NEG_HALF, Z[1], local_length);
    failure += check_ans(ONE,      Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 4b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 4b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 5a: (nvec > 1, nsum = 2), N_VLinearSumVectorArray
   * X[0][i] = c[0] X[0][i] + c[1] X[1][i]
   */

  /* fill vector data and set scaling factors */
  N_VConst(NEG_TWO, X[0][0]);
  N_VConst(TWO,     X[1][0]);

  N_VConst(TWO,  X[0][1]);
  N_VConst(HALF, X[1][1]);

  N_VConst(ZERO, X[0][2]);
  N_VConst(HALF, X[1][2]);

  c[0] = HALF;
  c[1] = TWO;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(3, 2, c, X, X[0]);
  sync_device(V);
  stop_time = get_time();

  /* X[0][i] should equal to +3, +2, +1 */
  if (ierr == 0) {
    failure  = check_ans(ONE+TWO, X[0][0], local_length);
    failure += check_ans(TWO,     X[0][1], local_length);
    failure += check_ans(ONE,     X[0][2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 5a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 5a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 5b: (nvec > 1, nsum = 2), N_VLinearSumVectorArray
   * Z[0] = c[0] X[0][i] + c[1] X[1][i]
   */

  /* fill vector data and set scaling factors */
  N_VConst(NEG_TWO, X[0][0]);
  N_VConst(TWO,     X[1][0]);

  N_VConst(TWO,  X[0][1]);
  N_VConst(HALF, X[1][1]);

  N_VConst(ZERO, X[0][2]);
  N_VConst(HALF, X[1][2]);

  c[0] = HALF;
  c[1] = TWO;

  N_VConst(ZERO, Z[0]);
  N_VConst(ZERO, Z[1]);
  N_VConst(ZERO, Z[2]);

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(3, 2, c, X, Z);
  sync_device(V);
  stop_time = get_time();

  /* X[0][i] should equal to +3, +2, +1 */
  if (ierr == 0) {
    failure  = check_ans(ONE+TWO, Z[0], local_length);
    failure += check_ans(TWO,     Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 5b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 5b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 6a: (nvec > 1, nsum > 2)
   * X[0][i] += c[1] X[1][i] + c[2] X[2][i]
   */

  /* fill vector data and set scaling factors */
  N_VConst(TWO,     X[0][0]);
  N_VConst(NEG_TWO, X[1][0]);
  N_VConst(NEG_ONE, X[2][0]);

  N_VConst(ONE,     X[0][1]);
  N_VConst(TWO,     X[1][1]);
  N_VConst(ONE,     X[2][1]);

  N_VConst(NEG_ONE, X[0][2]);
  N_VConst(TWO,     X[1][2]);
  N_VConst(TWO,     X[2][2]);

  c[0] = ONE;
  c[1] = NEG_HALF;
  c[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(3, 3, c, X, X[0]);
  sync_device(V);
  stop_time = get_time();

  /* X[0][i] should equal to +4, -1, -4 */
  if (ierr == 0) {
    failure  = check_ans(TWO+TWO,  X[0][0], local_length);
    failure += check_ans(NEG_ONE,  X[0][1], local_length);
    failure += check_ans(-TWO-TWO, X[0][2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 6a, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 6a \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 6b: (nvec > 1, nsum > 2)
   * X[0][i] = c[0] X[0][i] + c[1] X[1][i] + c[2] X[2][i]
   */

  /* fill vector data and set scaling factors */
  N_VConst(ONE,     X[0][0]);
  N_VConst(NEG_TWO, X[1][0]);
  N_VConst(NEG_ONE, X[2][0]);

  N_VConst(NEG_ONE, X[0][1]);
  N_VConst(TWO,     X[1][1]);
  N_VConst(ONE,     X[2][1]);

  N_VConst(HALF,    X[0][2]);
  N_VConst(TWO,     X[1][2]);
  N_VConst(ONE,     X[2][2]);

  c[0] = TWO;
  c[1] = HALF;
  c[2] = NEG_ONE;

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(3, 3, c, X, X[0]);
  sync_device(V);
  stop_time = get_time();

  /* X[0][i] should equal to +2, -2, +1 */
  if (ierr == 0) {
    failure  = check_ans(TWO,     X[0][0], local_length);
    failure += check_ans(NEG_TWO, X[0][1], local_length);
    failure += check_ans(ONE,     X[0][2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 6b, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 6b \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /*
   * Case 6c: (nvec > 1, nsum > 2)
   * Z[i] = c[0] X[0][i] + c[1] X[1][i] + c[2] X[2][i]
   */

  /* fill vector data and set scaling factors */
  N_VConst(ONE,     X[0][0]);
  N_VConst(NEG_TWO, X[1][0]);
  N_VConst(NEG_ONE, X[2][0]);

  N_VConst(NEG_ONE, X[0][1]);
  N_VConst(TWO,     X[1][1]);
  N_VConst(ONE,     X[2][1]);

  N_VConst(HALF,    X[0][2]);
  N_VConst(TWO,     X[1][2]);
  N_VConst(ONE,     X[2][2]);

  c[0] = TWO;
  c[1] = HALF;
  c[2] = NEG_ONE;

  N_VConst(ZERO, Z[0]);
  N_VConst(ZERO, Z[1]);
  N_VConst(ZERO, Z[2]);

  start_time = get_time();
  ierr = N_VLinearCombinationVectorArray(3, 3, c, X, Z);
  sync_device(V);
  stop_time = get_time();

  /* Z[i] should equal to +2, -2, +1 */
  if (ierr == 0) {
    failure  = check_ans(TWO,     Z[0], local_length);
    failure += check_ans(NEG_TWO, Z[1], local_length);
    failure += check_ans(ONE,     Z[2], local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VLinearCombinationVectorArray Case 6c, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VLinearCombinationVectorArray Case 6c \n");
  }

  /* find max time across all processes */
  maxt = max_time(V, stop_time - start_time);
  PRINT_TIME("N_VLinearCombinationVectorArray", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(Z, 3);
  N_VDestroyVectorArray(X[0], 3);
  N_VDestroyVectorArray(X[1], 3);
  N_VDestroyVectorArray(X[2], 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VDotProdLocal test
 * --------------------------------------------------------------------*/
int Test_N_VDotProdLocal(N_Vector X, N_Vector Y, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans, rmyid, locleninv;

  /* fill vector data */
  rmyid = (realtype) myid;
  locleninv = ONE/local_length;
  set_element_range(X, 0, local_length-1, rmyid);
  set_element_range(Y, 0, local_length-1, locleninv);

  start_time = get_time();
  ans = N_VDotProdLocal(X,Y);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal rmyid */
  failure = SUNRCompareTol(ans, rmyid, SUNRsqrt(UNIT_ROUNDOFF));

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdLocal, Proc %d\n", myid);
    printf("ans = %" FSYM " expected = %" FSYM "\n", ans, rmyid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdLocal\n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdLocal", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMaxNormLocal test
 * --------------------------------------------------------------------*/
int Test_N_VMaxNormLocal(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans, myidp1;

  /* fill vector data */
  myidp1 = (realtype) (myid+1);
  N_VConst(NEG_HALF, X);
  set_element(X, local_length-1, -myidp1);

  start_time = get_time();
  ans = N_VMaxNormLocal(X);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal myidp1 */
  failure = (ans < ZERO) ? 1 : SUNRCompareTol(ans, myidp1, SUNRsqrt(UNIT_ROUNDOFF));

  if (failure) {
    printf(">>> FAILED test -- N_VMaxNormLocal, Proc %d\n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMaxNormLocal\n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VMaxNormLocal", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMinLocal test
 * --------------------------------------------------------------------*/
int Test_N_VMinLocal(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans, negmyid;

  /* fill vector data */
  negmyid = (realtype) -myid;
  N_VConst(TWO, X);
  set_element(X, local_length-1, negmyid);

  start_time = get_time();
  ans = N_VMinLocal(X);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal negmyid */
  failure = SUNRCompareTol(ans, negmyid, SUNRsqrt(UNIT_ROUNDOFF));

  if (failure) {
    printf(">>> FAILED test -- N_VMinLocal, Proc %d\n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMinLocal\n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VMinLocal", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VL1NormLocal test
 * --------------------------------------------------------------------*/
int Test_N_VL1NormLocal(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans, val;

  /* fill vector data */
  val = -(ONE*myid)/local_length;
  set_element_range(X, 0, local_length-1, val);

  start_time = get_time();
  ans = N_VL1NormLocal(X);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal myid */
  failure = (ans < ZERO) ? 1 : SUNRCompareTol(ans, (realtype) myid,
                                       SUNRsqrt(UNIT_ROUNDOFF));

  if (failure) {
    printf(">>> FAILED test -- N_VL1NormLocal, Proc %d\n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VL1NormLocal\n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VL1NormLocal", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWSqrSumLocal test
 * --------------------------------------------------------------------*/
int Test_N_VWSqrSumLocal(N_Vector X, N_Vector W, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans, xval, wval;

  /* fill vector data */
  xval = SUNRsqrt(myid);
  wval = ONE/SUNRsqrt(local_length);
  set_element_range(X, 0, local_length-1, xval);
  set_element_range(W, 0, local_length-1, wval);

  start_time = get_time();
  ans = N_VWSqrSumLocal(X, W);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal myid */
  failure = (ans < ZERO) ? 1 : SUNRCompareTol(ans, (realtype) myid,
                                       SUNRsqrt(UNIT_ROUNDOFF));

  if (failure) {
    printf(">>> FAILED test -- N_VWSqrSumLocal, Proc %d\n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWSqrSumLocal\n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWL2NormLocal", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWSqrSumMaskLocal test
 * --------------------------------------------------------------------*/
int Test_N_VWSqrSumMaskLocal(N_Vector X, N_Vector W, N_Vector ID,
                             sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans, xval, wval;

  /* fill vector data */
  xval = SUNRsqrt(myid);
  wval = ONE/SUNRsqrt(local_length-1);
  set_element_range(X, 0, local_length-1, xval);
  set_element_range(W, 0, local_length-1, wval);

  /* use all elements except one */
  N_VConst(ONE, ID);
  set_element(ID, local_length-1, ZERO);

  start_time = get_time();
  ans = N_VWSqrSumMaskLocal(X, W, ID);
  sync_device(X);
  stop_time = get_time();

  /* ans should equal myid */
  failure = (ans < ZERO) ? 1 : SUNRCompareTol(ans, (realtype) myid,
                                       SUNRsqrt(UNIT_ROUNDOFF));

  if (failure) {
    printf(">>> FAILED test -- N_VWSqrSumMaskLocal, Proc %d\n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VWSqrSumMaskLocal\n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VWSqrSumMaskLocal", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VInvTestLocal test
 * --------------------------------------------------------------------*/
int Test_N_VInvTestLocal(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  realtype     xval;
  sunindextype i;
  booleantype  test;

  if (local_length < 2) {
    printf("Error Test_N_VInvTestLocal: Local vector length is %ld, length must be >= 2\n",
           (long int) local_length);
    return(1);
  }

  /*
   * Case 1: All elements Nonzero, z[i] = 1/x[i], return True
   */

  /* fill vector data */
  xval = ONE/(myid+2);
  set_element_range(X, 0, local_length-1, xval);
  N_VConst(ZERO, Z);

  start_time = get_time();
  test = N_VInvTestLocal(X, Z);
  sync_device(X);
  stop_time = get_time();

  /* Z should be vector of myid+2 */
  failure = check_ans(myid+2, Z, local_length);

  if (failure || !test) {
    printf(">>> FAILED test -- N_VInvTestLocal Case 1, Proc %d\n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VInvTestLocal Case 1\n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VInvTestLocal", maxt);

  /*
   * Case 2: Some elements Zero, z[i] = 1/x[i] for x[i] != 0, return False
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ZERO, Z);
  for (i=0; i < local_length; i++) {
    if (i % 2)
      set_element(X, i, HALF);
    else
      set_element(X, i, ZERO);
  }

  start_time = get_time();
  test = N_VInvTestLocal(X, Z);
  sync_device(X);
  stop_time = get_time();

  /* check return vector */
  for (i=0; i < local_length; i++) {
    if (i % 2) {
      if (get_element(Z, i) != TWO)
        failure = 1;
    } else {
      if (get_element(Z, i) != ZERO)
        failure = 1;
    }
  }

  if (failure || test) {
    printf(">>> FAILED test -- N_VInvTestLocal Case 2, Proc %d\n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VInvTestLocal Case 2\n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VInvTestLocal", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VConstrMaskLocal test
 * --------------------------------------------------------------------*/
int Test_N_VConstrMaskLocal(N_Vector C, N_Vector X, N_Vector M,
                            sunindextype local_length, int myid)
{
  int          fails = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype i;
  booleantype  test;

  if (local_length < 7) {
    printf("Error Test_N_VConstrMaskLocal: Local vector length is %ld, length must be >= 7\n",
           (long int) local_length);
    return(1);
  }

  /*
   * Case 1: Return True
   */

  /* fill vector data */
  for (i=0; i < local_length; i++) {
    set_element(M, i, NEG_ONE);

    switch(i % 7) {
    case 0 :
      /* c = -2, test for < 0*/
      set_element(C, i, NEG_TWO);
      set_element(X, i, NEG_TWO);
      break;

    case 1 :
      /* c = -1, test for <= 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, NEG_ONE);
      break;

    case 2 :
      /* c = -1, test for == 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, ZERO);
      break;

    case 3 :
      /* c = 0, no test */
      set_element(C, i, ZERO);
      set_element(X, i, HALF);
      break;

    case 4 :
      /* c = 1, test for == 0*/
      set_element(C, i, ONE);
      set_element(X, i, ZERO);
      break;

    case 5 :
      /* c = 1, test for >= 0*/
      set_element(C, i, ONE);
      set_element(X, i, ONE);
      break;

    case 6:
      /* c = 2, test for > 0 */
      set_element(C, i, TWO);
      set_element(X, i, TWO);
      break;
    }
  }

  start_time = get_time();
  test = N_VConstrMaskLocal(C, X, M);
  sync_device(X);
  stop_time = get_time();

  /* M should be vector of 0 */
  failure = check_ans(ZERO, M, local_length);

  if (failure || !test) {
    printf(">>> FAILED test -- N_VConstrMaskLocal Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VConstrMaskLocal Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VConstrMaskLocal", maxt);

  /*
   * Case 2: Return False
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  for (i=0; i < local_length; i++) {
    set_element(M, i, NEG_ONE);

    switch(i % 5) {
    case 0 :
      /* c = -2, test for < 0*/
      set_element(C, i, NEG_TWO);
      set_element(X, i, TWO);
      break;

    case 1 :
      /* c = -1, test for <= 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, ONE);
      break;

    case 2 :
      /* c = 0, no test */
      set_element(C, i, ZERO);
      set_element(X, i, HALF);
      break;

    case 3 :
      /* c = 1, test for >= 0*/
      set_element(C, i, ONE);
      set_element(X, i, NEG_ONE);
      break;

    case 4 :
      /* c = 2, test for > 0 */
      set_element(C, i, TWO);
      set_element(X, i, NEG_TWO);
      break;
    }
  }

  start_time = get_time();
  test = N_VConstrMaskLocal(C, X, M);
  sync_device(X);
  stop_time = get_time();

  /* check mask vector */
  for (i=0; i < local_length; i++) {

    if (i % 5 == 2){
      if (get_element(M, i) != ZERO)
        failure = 1;
    } else {
      if (get_element(M, i) != ONE)
        failure = 1;
    }
  }

  if (failure || test) {
    printf(">>> FAILED test -- N_VConstrMaskLocal Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VConstrMaskLocal Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VConstrMaskLocal", maxt);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMinQuotientLocal test
 * --------------------------------------------------------------------*/
int Test_N_VMinQuotientLocal(N_Vector NUM, N_Vector DENOM,
                             sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time, maxt;
  realtype ans;

  /*
   * Case 1: Pass
   */

  /* fill vector data */
  set_element_range(NUM, 0, local_length-1, TWO*(myid+2));
  N_VConst(TWO, DENOM);
  set_element(NUM, local_length-1, TWO*myid);

  start_time = get_time();
  ans = N_VMinQuotientLocal(NUM, DENOM);
  sync_device(NUM);
  stop_time = get_time();

  /* ans should equal myid */
  failure = SUNRCompareTol(ans, (realtype) myid, SUNRsqrt(UNIT_ROUNDOFF));

  if (failure) {
    printf(">>> FAILED test -- N_VMinQuotientLocal Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMinQuotientLocal Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(NUM, stop_time - start_time);
  PRINT_TIME("N_VMinQuotientLocal", maxt);

  /*
   * Case 2: Fail
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, NUM);
  N_VConst(ZERO, DENOM);

  start_time = get_time();
  ans = N_VMinQuotientLocal(NUM, DENOM);
  sync_device(NUM);
  stop_time = get_time();

  /* ans should equal BIG_REAL */
  failure = SUNRCompareTol(ans, BIG_REAL, SUNRsqrt(UNIT_ROUNDOFF));

  if (failure) {
    printf(">>> FAILED test -- N_VMinQuotientLocal Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VMinQuotientLocal Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(NUM, stop_time - start_time);
  PRINT_TIME("N_VMinQuotientLocal", maxt);

  return(fails);
}



/* ----------------------------------------------------------------------
 * N_VDotProdMultiLocal Test
 * --------------------------------------------------------------------*/
int Test_N_VDotProdMultiLocal(N_Vector X, sunindextype local_length, int myid)
{
  int    fails = 0, failure = 0, ierr = 0;
  double start_time, stop_time, maxt;

  N_Vector *V;
  realtype dotprods[3];

  /* create vectors for testing */
  V = N_VCloneVectorArray(3, X);

  /*
   * Case 1: d[0] = z . V[0], N_VDotProd
   */

  /* fill vector data */
  N_VConst(TWO,  X);
  N_VConst(HALF, V[0]);

  start_time = get_time();
  ierr = N_VDotProdMultiLocal(1, X, V, dotprods);
  sync_device(X);
  stop_time = get_time();

  /* dotprod[0] should equal the local vector length */
  if (ierr == 0)
    failure = SUNRCompare(dotprods[0], (realtype) local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdMultiLocal Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdMultiLocal Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdMultiLocal", maxt);

  /*
   * Case 2: d[i] = z . V[i], N_VDotProd
   */

  /* fill vector data */
  N_VConst(TWO,      X);
  N_VConst(NEG_HALF, V[0]);
  N_VConst(HALF,     V[1]);
  N_VConst(ONE,      V[2]);

  start_time = get_time();
  ierr = N_VDotProdMultiLocal(3, X, V, dotprods);
  sync_device(X);
  stop_time = get_time();

  /* dotprod[i] should equal -1, +1, and 2 times the local vector length */
  if (ierr == 0) {
    failure  = SUNRCompare(dotprods[0], (realtype) -1 * local_length);
    failure += SUNRCompare(dotprods[1], (realtype)      local_length);
    failure += SUNRCompare(dotprods[2], (realtype)  2 * local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdMultiLocal Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdMultiLocal Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdMultiLocal", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(V, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VDotProdMultiAllReduce Test
 * --------------------------------------------------------------------*/
int Test_N_VDotProdMultiAllReduce(N_Vector X, sunindextype local_length,
                                  int myid)
{
  int      fails = 0, failure = 0, ierr = 0;
  double   start_time, stop_time, maxt;

  sunindextype  global_length;
  N_Vector     *V;
  realtype      dotprods[3];

  /* only test if the operation is implemented, local vectors (non-MPI) do not
     provide this function */
  if (!(X->ops->nvdotprodmultiallreduce)) return 0;

  /* get global length */
  global_length = N_VGetLength(X);

  /* create vectors for testing */
  V = N_VCloneVectorArray(3, X);

  /*
   * Case 1: d[0] = z . V[0], N_VDotProd
   */

  /* fill vector data */
  N_VConst(TWO,  X);
  N_VConst(HALF, V[0]);

  start_time = get_time();
  ierr = N_VDotProdMultiLocal(1, X, V, dotprods);
  sync_device(X);
  stop_time = get_time();

  /* dotprod[0] should equal the local vector length */
  if (ierr == 0)
    failure = SUNRCompare(dotprods[0], (realtype) local_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdMultiAllReduce Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdMultiAllReduce Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdMultiLocal", maxt);

  /* perform the global reduction */
  start_time = get_time();
  ierr = N_VDotProdMultiAllReduce(1, X, dotprods);
  sync_device(X);
  stop_time = get_time();

  /* dotprod[0] should equal the global vector length */
  if (ierr == 0)
    failure = SUNRCompare(dotprods[0], (realtype) global_length);
  else
    failure = 1;

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdMultiAllReduce Case 1, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdMultiAllReduce Case 1 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdMultiAllReduce", maxt);

  /*
   * Case 2: d[i] = z . V[i], N_VDotProd
   */

  /* fill vector data */
  N_VConst(TWO,      X);
  N_VConst(NEG_HALF, V[0]);
  N_VConst(HALF,     V[1]);
  N_VConst(ONE,      V[2]);

  start_time = get_time();
  ierr = N_VDotProdMultiLocal(3, X, V, dotprods);
  sync_device(X);
  stop_time = get_time();

  /* dotprod[i] should equal -1, +1, and 2 times the local vector length */
  if (ierr == 0) {
    failure  = SUNRCompare(dotprods[0], (realtype) -1 * local_length);
    failure += SUNRCompare(dotprods[1], (realtype)      local_length);
    failure += SUNRCompare(dotprods[2], (realtype)  2 * local_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdMultiLocal Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdMultiLocal Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdMultiLocal", maxt);

  /* perform the global reduction */
  start_time = get_time();
  ierr = N_VDotProdMultiAllReduce(3, X, dotprods);
  sync_device(X);
  stop_time = get_time();

  /* dotprod[i] should equal -1, +1, and 2 times the global vector length */
  if (ierr == 0) {
    failure  = SUNRCompare(dotprods[0], (realtype) -1 * global_length);
    failure += SUNRCompare(dotprods[1], (realtype)      global_length);
    failure += SUNRCompare(dotprods[2], (realtype)  2 * global_length);
  } else {
    failure = 1;
  }

  if (failure) {
    printf(">>> FAILED test -- N_VDotProdMultiAllReduce Case 2, Proc %d \n", myid);
    fails++;
  } else if (myid == 0) {
    printf("PASSED test -- N_VDotProdMultiAllReduce Case 2 \n");
  }

  /* find max time across all processes */
  maxt = max_time(X, stop_time - start_time);
  PRINT_TIME("N_VDotProdMultiAllReduce", maxt);

  /* Free vectors */
  N_VDestroyVectorArray(V, 3);

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VBufSize test
 * --------------------------------------------------------------------*/
int Test_N_VBufSize(N_Vector x, sunindextype local_length, int myid)
{
  int          flag = 0;
  double       start_time, stop_time, maxt;
  sunindextype size;

  /* get buffer size */
  start_time = get_time();
  flag = N_VBufSize(x, &size);
  sync_device(x);
  stop_time = get_time();

  /* check return value */
  if (flag != 0) {
    printf(">>> FAILED test -- N_VBufSize returned %d, Proc %d \n", flag, myid);
    return(1);
  }

  /* check buffer size */
  if (size != local_length * ((sunindextype)sizeof(realtype))) {
    printf(">>> FAILED test -- N_VBufSize, Proc %d \n", myid);
    return(1);
  }

  if (myid == 0) {
    printf("PASSED test -- N_VBufSize\n");
  }

  /* find max time across all processes */
  maxt = max_time(x, stop_time - start_time);
  PRINT_TIME("N_VBufSize", maxt);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VBufPack test
 * --------------------------------------------------------------------*/
int Test_N_VBufPack(N_Vector x, sunindextype local_length, int myid)
{
  int          flag = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype i, size;
  realtype     *buf;

  /* get buffer size */
  flag = N_VBufSize(x, &size);
  if (flag != 0) {
    printf(">>> FAILED test -- N_VBufSize returned %d, Proc %d \n", flag, myid);
    return(1);
  }

  /* create and initialize buffer */
  buf = NULL;
  buf = (realtype*) malloc((size_t)size);
  if (buf == NULL) {
    printf(">>> FAILED test -- malloc failed, Proc %d \n", myid);
    return(1);
  }

  for (i = 0; i < local_length; i++) {
    buf[i] = ZERO;
  }

  /* set vector data */
  N_VConst(ONE, x);

  /* fill buffer */
  start_time = get_time();
  flag = N_VBufPack(x, (void*)buf);
  sync_device(x);
  stop_time = get_time();

  if (flag != 0) {
    free(buf);
    printf(">>> FAILED test -- N_VBufPack returned %d, Proc %d \n", flag, myid);
    return(1);
  }

  /* check buffer values */
  for(i = 0; i < local_length; i++) {
    failure += SUNRCompare(buf[i], ONE);
  }

  if (failure) {
    free(buf);
    printf(">>> FAILED test -- N_VBufPack failed, Proc %d \n", myid);
    return(1);
  } else if (myid == 0) {
    printf("PASSED test -- N_VBufPack\n");
  }

  /* find max time across all processes */
  maxt = max_time(x, stop_time - start_time);
  PRINT_TIME("N_VBufPack", maxt);

  /* free buffer */
  free(buf);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VBufUnpack test
 * --------------------------------------------------------------------*/
int Test_N_VBufUnpack(N_Vector x, sunindextype local_length, int myid)
{
  int          flag = 0, failure = 0;
  double       start_time, stop_time, maxt;
  sunindextype i, size;
  realtype     *buf;

  /* get buffer size */
  flag = N_VBufSize(x, &size);
  if (flag != 0) {
    printf(">>> FAILED test -- N_VBufSize returned %d, Proc %d \n", flag, myid);
    return(1);
  }

  /* create and initialize buffer */
  buf = NULL;
  buf = (realtype*) malloc((size_t)size);
  if (buf == NULL) {
    printf(">>> FAILED test -- malloc failed, Proc %d \n", myid);
    return(1);
  }

  for(i = 0; i < local_length; i++) {
    buf[i] = ONE;
  }

  /* clear the input vector values */
  N_VConst(ZERO, x);

  /* fill vector data */
  start_time = get_time();
  flag = N_VBufUnpack(x, (void*)buf);
  sync_device(x);
  stop_time = get_time();

  if (flag != 0) {
    free(buf);
    printf(">>> FAILED test -- N_VBufUnPack returned %d, Proc %d \n", flag, myid);
    return(1);
  }

  /* x should be vector of ones */
  failure = check_ans(ONE, x, local_length);

  if (failure) {
    free(buf);
    printf(">>> FAILED test -- N_VBufUnpack failed, Proc %d \n", myid);
    return(1);
  } else if (myid == 0) {
    printf("PASSED test -- N_VBufUnpack\n");
  }

  /* find max time across all processes */
  maxt = max_time(x, stop_time - start_time);
  PRINT_TIME("N_VBufUnpack", maxt);

  /* free buffer */
  free(buf);

  return(0);
}


/* ======================================================================
 * Private functions
 * ====================================================================*/


void SetTiming(int onoff, int myid)
{
#if defined(SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime(CLOCK_MONOTONIC, &spec);
  base_time_tv_sec = spec.tv_sec;

  clock_getres(CLOCK_MONOTONIC, &spec);
  if (myid == 0)
    printf("Timer resolution: %ld ns = %g s\n",
           spec.tv_nsec, ((double)(spec.tv_nsec) / 1E9));
#endif

  /* only print from the root process */
  print_time = (myid == 0) ? onoff : 0;
}

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
static double get_time()
{
  double time;
#if defined(SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime(CLOCK_MONOTONIC, &spec);
  time = (double)(spec.tv_sec - base_time_tv_sec) + ((double)(spec.tv_nsec) / 1E9);
#else
  time = 0;
#endif
  return time;
}
