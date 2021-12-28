/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * These test functions are designed to evaluate the performance of
 * an NVECTOR module implementation. They do not check for accuracy.
 * To test the accuracy of an implementation see the test_nvector.c
 * file.
 * -----------------------------------------------------------------*/

/* Minimum POSIX version needed for struct timespec and clock_monotonic */
#if !defined(_POSIX_C_SOURCE) || (_POSIX_C_SOURCE < 199309L)
#define _POSIX_C_SOURCE 199309L
#endif

/* POSIX timers */
#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
#include <time.h>
#include <stddef.h>
#include <unistd.h>
#endif

#include <sundials/sundials_config.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include "test_nvector_performance.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>


/* private functions */
static double get_time();
static void time_stats(N_Vector X, double *times, int start, int ntimes,
                       double *avg, double *sdev, double *min, double *max);

int print_time = 0; /* flag for printing timing data */
int nextra = 1;     /* number of extra tests to perform and ignore in average */

#if defined(SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
time_t base_time_tv_sec = 0; /* Base time; makes time values returned
                                by get_time easier to read when
                                printed since they will be zero
                                based.
                              */
#endif

#define FMT1 "%33s %22.15e %22.15e %22.15e %22.15e\n"
#define PRINT_TIME1(test, time, sdev, min, max) \
  if(print_time) printf(FMT1, test, time, sdev, min, max)

#define FMT2 "%33s %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n"
#define PRINT_TIME2(test, time1, sdev1, min1, max1, time2, sdev2, min2, max2) \
  if(print_time) printf(FMT2, test, time1, sdev1, min1, max1, time2, sdev2, min2, max2)


/* -----------------------------------------------------------------------------
 * N_VLinearSum Tests
 * ---------------------------------------------------------------------------*/
int Test_N_VLinearSum(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  realtype a, b;
  int      i;
  N_Vector Y, Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Y = N_VClone(X);
  Z = N_VClone(X);

  /*
   * Case 1a: y = x + y, (Vaxpy Case 1)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector with random data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VLinearSum(ONE, X, ONE, Y, Y);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-1a", avgtime, sdevtime, mintime, maxtime);


  /*
   * Case 1b: y = -x + y, (Vaxpy Case 2)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VLinearSum(NEG_ONE, X, ONE, Y, Y);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-1b", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 1c: y = ax + y, (Vaxpy Case 3)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    a = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(a, X, ONE, Y, Y);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-1c", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 2a: x = x + y, (Vaxpy Case 1)
   */

  for (i=0; i < ntests+nextra; i++) {

    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VLinearSum(ONE, X, ONE, Y, X);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-2a", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 2b: x = x - y, (Vaxpy Case 2)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VLinearSum(ONE, X, NEG_ONE, Y, X);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-2b", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 2c: x = x + by, (Vaxpy Case 3)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    b = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(ONE, X, b, Y, X);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-2c", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 3: z = x + y, (VSum)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VLinearSum(ONE, X, ONE, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-3", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 4a: z = x - y, (VDiff)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VLinearSum(ONE, X, NEG_ONE, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-4a", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 4b: z = -x + y, (VDiff)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VLinearSum(NEG_ONE, X, ONE, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-4b", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 5a: z = x + by, (VLin1)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    b = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(ONE, X, b, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-5a", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 5b: z = ax + y, (VLin1)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    a = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(a, X, ONE, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-5b", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 6a: z = -x + by, (VLin2)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    b = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(NEG_ONE, X, b, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-6a", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 6b: z = ax - y, (VLin2)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    a = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(a, X, NEG_ONE, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-6b", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 7: z = a(x + y), (VScaleSum)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    a = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(a, X, a, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-7", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 8: z = a(x - y), (VScaleDiff)
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    a = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
    b = -1.0*a;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(a, X, b, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-8", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 9: z = ax + by, All Other Cases
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    a = 2.0*((realtype)rand() / (realtype)RAND_MAX) - 1.0;
    b = 2.0*((realtype)rand() / (realtype)RAND_MAX) - 1.0;

    ClearCache();
    start_time = get_time();
    N_VLinearSum(a, X, b, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VLinearSum-9", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Y);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VConst Test
 * --------------------------------------------------------------------*/
int Test_N_VConst(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  realtype c;
  int      i;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  for (i=0; i < ntests+nextra; i++) {
    c = 2.0*((realtype)rand() / (realtype)RAND_MAX) - 1.0;

    ClearCache();
    start_time = get_time();
    N_VConst(c, X);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VConst", avgtime, sdevtime, mintime, maxtime);

  free(times);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VProd Test
 * --------------------------------------------------------------------*/
int Test_N_VProd(N_Vector X, sunindextype local_length, int ntests)
{
  double start_time, stop_time;
  double *times;
  double avgtime, sdevtime, mintime, maxtime;
  int    i;
  N_Vector Y, Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Y = N_VClone(X);
  Z = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VProd(X, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VProd", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Y);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VDiv Test
 * --------------------------------------------------------------------*/
int Test_N_VDiv(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;
  N_Vector Y, Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Y = N_VClone(X);
  Z = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, ONE, TEN);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VDiv(X, Y, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VDiv", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Y);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VScale Tests
 * --------------------------------------------------------------------*/
int Test_N_VScale(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  realtype c;
  int      i;
  N_Vector Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Z = N_VClone(X);

  /*
   * Case 1: x = cx, VScaleBy
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    c = 2.0*((realtype)rand() / (realtype)RAND_MAX) - 1.0;

    ClearCache();
    start_time = get_time();
    N_VScale(c, X, X);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VScale-1", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 2: z = x, VCopy
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VScale(ONE, X, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VScale-2", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 3: z = -x, VNeg
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VScale(NEG_ONE, X, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VScale-3", avgtime, sdevtime, mintime, maxtime);

  /*
   * Case 4: z = cx, All other cases
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    c = 2.0*((realtype)rand() / (realtype)RAND_MAX) - 1.0;

    ClearCache();
    start_time = get_time();
    N_VScale(c, X, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VScale-4", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VAbs Test
 * --------------------------------------------------------------------*/
int Test_N_VAbs(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;
  N_Vector Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Z = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VAbs(X, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VAbs", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VInv Test
 * --------------------------------------------------------------------*/
int Test_N_VInv(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;
  N_Vector Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Z = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, ONE, TEN);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VInv(X, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VInv", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VAddConst Test
 * --------------------------------------------------------------------*/
int Test_N_VAddConst(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  realtype c;
  int      i;
  N_Vector Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Z = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    c = 2.0*((realtype)rand() / (realtype)RAND_MAX) - 1.0;

    ClearCache();
    start_time = get_time();
    N_VAddConst(X, c, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VAddConst", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VDotProd Test
 * --------------------------------------------------------------------*/
int Test_N_VDotProd(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;
  N_Vector Y;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Y = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VDotProd(X, Y);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VDotProd", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Y);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VMaxNorm Test
 * --------------------------------------------------------------------*/
int Test_N_VMaxNorm(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VMaxNorm(X);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VMaxNorm", avgtime, sdevtime, mintime, maxtime);

  free(times);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNorm Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNorm(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;
  N_Vector W;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  W = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(W, local_length, ONE, TWO);

    ClearCache();
    start_time = get_time();
    N_VWrmsNorm(X, W);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VWrmsNorm", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(W);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNormMask Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNormMask(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;
  N_Vector W, ID;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  W  = N_VClone(X);
  ID = N_VClone(X);

  /*
   * Case 2: use no elements, ID = 0
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(W, local_length, ONE, TWO);
    N_VRandZeroOne(ID, local_length);

    ClearCache();
    start_time = get_time();
    N_VWrmsNormMask(X, W, ID);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VWrmsNormMask", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(W);
  N_VDestroy(ID);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VMin Test
 * --------------------------------------------------------------------*/
int Test_N_VMin(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VMin(X);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VMin", avgtime, sdevtime, mintime, maxtime);

  free(times);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VWL2Norm Test
 * --------------------------------------------------------------------*/
int Test_N_VWL2Norm(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;
  N_Vector W;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  W = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(W, local_length, ONE, TWO);

    ClearCache();
    start_time = get_time();
    N_VWL2Norm(X, W);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VWL2Norm", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(W);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VL1Norm Test
 * --------------------------------------------------------------------*/
int Test_N_VL1Norm(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VL1Norm(X);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VL1Norm", avgtime, sdevtime, mintime, maxtime);

  free(times);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VCompare
 * --------------------------------------------------------------------*/
int Test_N_VCompare(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  realtype c;
  int      i;
  N_Vector Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Z = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);
    c = ((realtype)rand() / (realtype)RAND_MAX);

    ClearCache();
    start_time = get_time();
    N_VCompare(c, X, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VCompare", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VInvTest
 * --------------------------------------------------------------------*/
int Test_N_VInvTest(N_Vector X, sunindextype local_length, int ntests)
{
  double      start_time, stop_time;
  double      *times;
  double      avgtime, sdevtime, mintime, maxtime;
  int         i;
  N_Vector    Z;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Z = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VConst(ZERO, Z);

    ClearCache();
    start_time = get_time();
    N_VInvTest(X, Z);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VInvTest", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Z);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VConstrMask
 * --------------------------------------------------------------------*/
int Test_N_VConstrMask(N_Vector X, sunindextype local_length, int ntests)
{
  double      start_time, stop_time;
  double      *times;
  double      avgtime, sdevtime, mintime, maxtime;
  int         i;
  N_Vector    C, M;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  C = N_VClone(X);
  M = N_VClone(X);

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRandConstraints(C, local_length);
    N_VConst(ZERO, M);

    ClearCache();
    start_time = get_time();
    N_VConstrMask(C, X, M);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VConstrMask", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(C);
  N_VDestroy(M);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VMinQuotient Test
 * --------------------------------------------------------------------*/
int Test_N_VMinQuotient(N_Vector X, sunindextype local_length, int ntests)
{
  double   start_time, stop_time;
  double   *times;
  double   avgtime, sdevtime, mintime, maxtime;
  int      i;
  N_Vector Y;

  times = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors */
  Y = N_VClone(X);

  /*
   * Case 1: Pass
   */

  for (i=0; i < ntests+nextra; i++) {
    /* fill vector data */
    N_VRand(X, local_length, NEG_ONE, ONE);
    N_VRand(Y, local_length, NEG_ONE, ONE);

    ClearCache();
    start_time = get_time();
    N_VMinQuotient(X, Y);
    sync_device(X);
    stop_time = get_time();

    times[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, times, nextra, ntests, &avgtime, &sdevtime, &mintime, &maxtime);
  PRINT_TIME1("N_VMinQuotient", avgtime, sdevtime, mintime, maxtime);

  /* Free vectors */
  free(times);
  N_VDestroy(Y);

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VLinearCombination_Serial Test
 * --------------------------------------------------------------------*/
int Test_N_VLinearCombination(N_Vector X, sunindextype local_length, int nvecs, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, ier;
  realtype *c;
  N_Vector *Y;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors and array of scaling factors */
  c = (realtype*) malloc(nvecs*sizeof(realtype));
  Y = N_VCloneVectorArray(nvecs, X);

  /*
   * Case 1: Y[0] = sum c[i] Y[i], c[0] = 1
   */

  /* fill vector data */
  c[0] = ONE;
  N_VRand(Y[0], local_length, NEG_ONE, ONE);

  for (j=1; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    if (nvecs == 1) {
      ClearCache();
      start_time = get_time();
      N_VScale(ONE, Y[0], Y[0]);
      sync_device(X);
      stop_time = get_time();
    } else {
      ClearCache();
      start_time = get_time();
      for (j=1; j < nvecs; j++)
        N_VLinearSum(ONE, Y[0], c[j], Y[j], Y[0]);
      sync_device(X);
      stop_time = get_time();
    }

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  c[0] = ONE;
  N_VRand(Y[0], local_length, NEG_ONE, ONE);

  for (j=1; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VLinearCombination(nvecs, c, Y, Y[0]);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VLinearCombination-1", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);


  /*
   * Case 2: Y[0] = sum c[i] Y[i], c[0] != 1
   */

  /* fill vector data */
  c[0] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  N_VRand(Y[0], local_length, NEG_ONE, ONE);

  for (j=1; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    if (nvecs == 1) {
      ClearCache();
      start_time = get_time();
      N_VScale(c[0], Y[0], Y[0]);
      sync_device(X);
      stop_time = get_time();
    } else {
      ClearCache();
      start_time = get_time();
      N_VScale(c[0], Y[0], Y[0]);
      for (j=1; j < nvecs; j++)
        N_VLinearSum(ONE, Y[0], c[j], Y[j], Y[0]);
      sync_device(X);
      stop_time = get_time();
    }

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  c[0] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  N_VRand(Y[0], local_length, NEG_ONE, ONE);

  for (j=1; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VLinearCombination(nvecs, c, Y, Y[0]);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VLinearCombination-2", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);


  /*
   * Case 3: X = sum c[i] Y[i]
   */

  /* fill vector data */
  for (j=0; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    if (nvecs == 1) {
      ClearCache();
      start_time = get_time();
      N_VScale(c[0], Y[0], Y[0]);
      sync_device(X);
      stop_time = get_time();
    } else {
      ClearCache();
      start_time = get_time();
      N_VScale(c[0], Y[0], X);
      for (j=1; j < nvecs; j++)
      N_VLinearSum(ONE, Y[0], c[j], Y[j], X);
      sync_device(X);
      stop_time = get_time();
    }

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  for (j=0; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VLinearCombination(nvecs, c, Y, X);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VLinearCombination-3", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  free(c);
  N_VDestroyVectorArray(Y, nvecs);

  return(ier);
}


/* ----------------------------------------------------------------------
 * N_VScaleaddmulti Test
 * --------------------------------------------------------------------*/
int Test_N_VScaleAddMulti(N_Vector X, sunindextype local_length, int nvecs, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, ier;
  realtype *c;
  N_Vector *Y, *Z;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors and array of scaling factors */
  c = (realtype*) malloc(nvecs*sizeof(realtype));
  Y = N_VCloneVectorArray(nvecs, X);
  Z = N_VCloneVectorArray(nvecs, X);

  /*
   * Case 1: Y[i] = c[i] x + Y[i], N_VScaleAddMulti
   */

  /* fill vector data */
  N_VRand(X, local_length, NEG_ONE, ONE);
  for (j=0; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      N_VLinearSum(c[j], X, ONE, Y[j], Y[j]);
    sync_device(X);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  N_VRand(X, local_length, NEG_ONE, ONE);
  for (j=0; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VScaleAddMulti(nvecs, c, X, Y, Y);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VScaleAddMulti-1", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);


  /*
   * Case 2: Z[i] = c[i] x + Y[i], N_VScaleAddMulti
   */

  /* fill vector data */
  N_VRand(X, local_length, NEG_ONE, ONE);
  for (j=0; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      N_VLinearSum(c[j], X, ONE, Y[j], Z[j]);
    sync_device(X);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  N_VRand(X, local_length, NEG_ONE, ONE);
  for (j=0; j < nvecs; j++) {
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VScaleAddMulti(nvecs, c, X, Y, Z);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VScaleAddMulti-2", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  free(c);
  N_VDestroyVectorArray(Y, nvecs);
  N_VDestroyVectorArray(Z, nvecs);

  return(ier);
}


/* ----------------------------------------------------------------------
 * N_VDotProdMulti Test
 * --------------------------------------------------------------------*/
int Test_N_VDotProdMulti(N_Vector X, sunindextype local_length, int nvecs, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, ier;
  realtype *c;
  N_Vector *Y;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create additional nvectors and array of dot products */
  c = (realtype*) malloc(nvecs*sizeof(realtype));
  Y = N_VCloneVectorArray(nvecs, X);

  /* fill vector data */
  N_VRand(X, local_length, NEG_ONE, ONE);
  for (j=0; j < nvecs; j++)
    N_VRand(Y[j], local_length, NEG_ONE, ONE);

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      c[j] = N_VDotProd(X, Y[j]);
    sync_device(X);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  N_VRand(X, local_length, NEG_ONE, ONE);
  for (j=0; j < nvecs; j++)
    N_VRand(Y[j], local_length, NEG_ONE, ONE);

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VDotProdMulti(nvecs, X, Y, c);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VDotProdMulti", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  free(c);
  N_VDestroyVectorArray(Y, nvecs);

  return(ier);
}


/* ----------------------------------------------------------------------
 * N_VLinearSumVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VLinearSumVectorArray(N_Vector V, sunindextype local_length,
                                 int nvecs, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, ier;
  realtype a, b;
  N_Vector *X, *Y, *Z;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create nvectors arrays */
  X = N_VCloneVectorArray(nvecs, V);
  Y = N_VCloneVectorArray(nvecs, V);
  Z = N_VCloneVectorArray(nvecs, V);

  /*
   * Case 1: Z[i] = a X[i] + b Y[i]
   */

  /* fill vector data */
  a = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  b = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  for (j=0; j < nvecs; j++) {
    N_VRand(X[j], local_length, NEG_ONE, ONE);
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      N_VLinearSum(a, X[j], b, Y[j], Z[j]);
    sync_device(V);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  a = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  b = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  for (j=0; j < nvecs; j++) {
    N_VRand(X[j], local_length, NEG_ONE, ONE);
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VLinearSumVectorArray(nvecs, a, X, b, Y, Z);
    sync_device(V);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(V, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(V, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VLinearSumVectorArray", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  N_VDestroyVectorArray(X, nvecs);
  N_VDestroyVectorArray(Y, nvecs);
  N_VDestroyVectorArray(Z, nvecs);

  return(ier);
}


/* ----------------------------------------------------------------------
 * N_VScaleVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VScaleVectorArray(N_Vector X, sunindextype local_length,
                             int nvecs, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, ier;
  realtype *c;
  N_Vector *Y, *Z;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create nvectors arrays and scaling factor array */
  c = (realtype*) malloc(nvecs*sizeof(realtype));
  Y = N_VCloneVectorArray(nvecs, X);
  Z = N_VCloneVectorArray(nvecs, X);

  /*
   * Case 1: Y[j] = c[j] Y[j]
   */

  /* fill vector data */
  for (j=0; j < nvecs; j++) {
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      N_VScale(c[j], Y[j], Y[j]);
    sync_device(X);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  for (j=0; j < nvecs; j++) {
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VScaleVectorArray(nvecs, c, Y, Y);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VScaleVectorArray-1", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);


  /*
   * Case 2: Z[j] = c[j] Y[j]
   */

  /* fill vector data */
  for (j=0; j < nvecs; j++) {
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      N_VScale(c[j], Y[j], Z[j]);
    sync_device(X);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  for (j=0; j < nvecs; j++) {
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
    N_VRand(Y[j], local_length, NEG_ONE, ONE);
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VScaleVectorArray(nvecs, c, Y, Z);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VScaleVectorArray-2", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  free(c);
  N_VDestroyVectorArray(Y, nvecs);
  N_VDestroyVectorArray(Z, nvecs);

  return(ier);
 }


/* ----------------------------------------------------------------------
 * N_VConstVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VConstVectorArray(N_Vector X, sunindextype local_length,
                             int nvecs, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, ier;
  realtype c;
  N_Vector *Y;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create nvectors array */
  Y = N_VCloneVectorArray(nvecs, X);

  /*
   * Case: Y[j] = c
   */

  /* fill vector data */
  c = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  for (j=0; j < nvecs; j++)
    N_VRand(Y[j], local_length, NEG_ONE, ONE);

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      N_VConst(c, Y[j]);
    sync_device(X);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  c = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;
  for (j=0; j < nvecs; j++)
    N_VRand(Y[j], local_length, NEG_ONE, ONE);

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VConstVectorArray(nvecs, c, Y);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VConstVectorArray", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  N_VDestroyVectorArray(Y, nvecs);

  return(ier);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNormVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNormVectorArray(N_Vector X, sunindextype local_length,
                                int nvecs, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, ier;
  realtype *c;
  N_Vector *Z, *W;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create nvectors arrays and array for norms */
  c = (realtype*) malloc(nvecs*sizeof(realtype));
  Z = N_VCloneVectorArray(nvecs, X);
  W = N_VCloneVectorArray(nvecs, X);

  /*
   * Case: nrm[i] = ||Z[i]||
   */

  /* fill vector data */
  for (j=0; j < nvecs; j++) {
    N_VRand(Z[j], local_length, NEG_ONE, ONE);
    N_VRand(W[j], local_length, ONE, TWO);
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      c[j] = N_VWrmsNorm(Z[j], W[j]);
    sync_device(X);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  for (j=0; j < nvecs; j++) {
    N_VRand(Z[j], local_length, NEG_ONE, ONE);
    N_VRand(W[j], local_length, ONE, TWO);
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VWrmsNormVectorArray(nvecs, Z, W, c);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VWrmsNormVectorArray", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  free(c);
  N_VDestroyVectorArray(Z, nvecs);
  N_VDestroyVectorArray(W, nvecs);

  return(ier);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNormVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNormMaskVectorArray(N_Vector X, sunindextype local_length,
                                    int nvecs, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, ier;
  realtype *c;
  N_Vector *Z, *W, ID;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create nvectors arrays and array for norms */
  c  = (realtype*) malloc(nvecs*sizeof(realtype));
  Z  = N_VCloneVectorArray(nvecs, X);
  W  = N_VCloneVectorArray(nvecs, X);
  ID = N_VClone(X);

  /*
   * Case: nrm[i] = ||Z[i]||
   */

  /* fill vector data */
  N_VRandZeroOne(ID, local_length);
  for (j=0; j < nvecs; j++) {
    N_VRand(Z[j], local_length, NEG_ONE, ONE);
    N_VRand(W[j], local_length, ONE, TWO);
  }

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (j=0; j < nvecs; j++)
      c[j] = N_VWrmsNormMask(Z[j], W[j], ID);
    sync_device(X);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill vector data */
  N_VRandZeroOne(ID, local_length);
  for (j=0; j < nvecs; j++) {
    N_VRand(Z[j], local_length, NEG_ONE, ONE);
    N_VRand(W[j], local_length, ONE, TWO);
  }

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VWrmsNormMaskVectorArray(nvecs, Z, W, ID, c);
    sync_device(X);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(X, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(X, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VWrmsNormMaskVectorArray", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  free(c);
  N_VDestroyVectorArray(Z, nvecs);
  N_VDestroyVectorArray(W, nvecs);
  N_VDestroy(ID);

  return(ier);
}


/* ----------------------------------------------------------------------
 * N_VScaleAddMultiVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VScaleAddMultiVectorArray(N_Vector V, sunindextype local_length,
                                     int nvecs, int nsums, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, k, ier;
  realtype *c;
  N_Vector *X, **Y, **Z;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create nvectors arrays and array for norms */
  c = (realtype*) malloc(nsums*sizeof(realtype));
  X = N_VCloneVectorArray(nvecs, V);

  Y = (N_Vector**) malloc(nsums*sizeof(N_Vector*));
  Z = (N_Vector**) malloc(nsums*sizeof(N_Vector*));
  for (j=0; j < nsums; j++) {
    Y[j] = N_VCloneVectorArray(nvecs, V);
    Z[j] = N_VCloneVectorArray(nvecs, V);
  }

  /*
   * Case 1: Y[j][k] = c[j] X[k] + Y[j][k]
   */

  /* fill data */
  for (j=0; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (k=0; k < nvecs; k++)
    N_VRand(X[k], local_length, NEG_ONE, ONE);

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(Y[j][k], local_length, NEG_ONE, ONE);

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (k=0; k < nvecs; k++)
      for (j=0; j < nsums; j++)
        N_VLinearSum(c[j], X[k], ONE, Y[j][k], Y[j][k]);
    sync_device(V);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill data */
  for (j=0; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (k=0; k < nvecs; k++)
    N_VRand(X[k], local_length, NEG_ONE, ONE);

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(Y[j][k], local_length, NEG_ONE, ONE);

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VScaleAddMultiVectorArray(nvecs, nsums, c, X, Y, Y);
    sync_device(V);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(V, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(V, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VScaleAddMultiVectorArray-1", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);


  /*
   * Case 2: Z[j][k] = c[j] X[k] + Y[j][k]
   */

  /* fill data */
  for (j=0; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (k=0; k < nvecs; k++)
    N_VRand(X[k], local_length, NEG_ONE, ONE);

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(Y[j][k], local_length, NEG_ONE, ONE);

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    for (k=0; k < nvecs; k++)
      for (j=0; j < nsums; j++)
        N_VLinearSum(c[j], X[k], ONE, Y[j][k], Z[j][k]);
    sync_device(V);
    stop_time = get_time();

    utimes[i] = stop_time - start_time;
  }

  /* fill data */
  for (j=0; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (k=0; k < nvecs; k++)
    N_VRand(X[k], local_length, NEG_ONE, ONE);

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(Y[j][k], local_length, NEG_ONE, ONE);

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VScaleAddMultiVectorArray(nvecs, nsums, c, X, Y, Z);
    sync_device(V);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(V, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(V, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VScaleAddMultiVectorArray-2", favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  free(c);
  N_VDestroyVectorArray(X, nvecs);

  for (j=0; j < nsums; j++) {
    N_VDestroyVectorArray(Y[j], nvecs);
    N_VDestroyVectorArray(Z[j], nvecs);
  }
  free(Y);
  free(Z);

  return(ier);
}


/* ----------------------------------------------------------------------
 * N_VLinearCombinationVectorArray Test
 * --------------------------------------------------------------------*/
int Test_N_VLinearCombinationVectorArray(N_Vector V, sunindextype local_length,
                                         int nvecs, int nsums, int ntests)
{
  double   start_time, stop_time;
  double   favgtime, fsdevtime, fmintime, fmaxtime;
  double   uavgtime, usdevtime, umintime, umaxtime;
  double   *ftimes, *utimes;
  int      i, j, k, ier;
  realtype *c;
  N_Vector **X, *Z;

  /* allocate timing arrays */
  ftimes = (double*) malloc((ntests+nextra)*sizeof(double));
  utimes = (double*) malloc((ntests+nextra)*sizeof(double));

  /* create nvectors arrays and array for norms */
  c = (realtype*) malloc(nsums*sizeof(realtype));
  Z = N_VCloneVectorArray(nvecs, V);

  X = (N_Vector**) malloc(nsums*sizeof(N_Vector*));
  for (j=0; j < nsums; j++)
    X[j] = N_VCloneVectorArray(nvecs, V);

  /*
   * Case 1: X[0][k] = sum c[j] X[j][k], c[0] = 1
   */

  /* fill data */
  c[0] = ONE;
  for (j=1; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(X[j][k], local_length, NEG_ONE, ONE);

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    if (nsums == 1) {
      ClearCache();
      start_time = get_time();
      for (k=0; k < nvecs; k++)
        N_VScale(ONE, X[0][k], X[0][k]);
      sync_device(V);
      stop_time = get_time();
    } else {
      ClearCache();
      start_time = get_time();
      for (k=0; k < nvecs; k++)
        for (j=1; j < nsums; j++)
          N_VLinearSum(ONE, X[0][k], c[j], X[j][k], X[0][k]);
      sync_device(V);
      stop_time = get_time();
    }

    utimes[i] = stop_time - start_time;
  }

  /* fill data */
  c[0] = ONE;
  for (j=1; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(X[j][k], local_length, NEG_ONE, ONE);

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VLinearCombinationVectorArray(nvecs, nsums, c, X, X[0]);
    sync_device(V);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(V, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(V, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VLinearCombinationVectorArray-1",
              favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);


  /*
   * Case 2: X[0][k] = sum c[j] X[j][k], c[0] != 1
   */

  /* fill data */
  for (j=0; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(X[j][k], local_length, NEG_ONE, ONE);

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    if (nsums == 1) {
      ClearCache();
      start_time = get_time();
      for (k=0; k < nvecs; k++)
        N_VScale(c[0], X[0][k], X[0][k]);
      sync_device(V);
      stop_time = get_time();
    } else {
      ClearCache();
      start_time = get_time();
      for (k=0; k < nvecs; k++) {
        N_VScale(c[0], X[0][k], X[0][k]);
        for (j=1; j < nsums; j++) {
          N_VLinearSum(ONE, X[0][k], c[j], X[j][k], X[0][k]);
        }
      }
      sync_device(V);
      stop_time = get_time();
    }

    utimes[i] = stop_time - start_time;
  }

  /* fill data */
  for (j=0; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(X[j][k], local_length, NEG_ONE, ONE);

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VLinearCombinationVectorArray(nvecs, nsums, c, X, X[0]);
    sync_device(V);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(V, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(V, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VLinearCombinationVectorArray-2",
              favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);


  /*
   * Case 3: Z[j][k] = sum c[j] X[j][k]
   */

  /* fill data */
  for (j=0; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(X[j][k], local_length, NEG_ONE, ONE);

  /* unfused operation */
  for (i=0; i < ntests+nextra; i++) {

    if (nsums == 1) {
      ClearCache();
      start_time = get_time();
      for (k=0; k < nvecs; k++)
        N_VScale(c[0], X[0][k], Z[k]);
      sync_device(V);
      stop_time = get_time();
    } else {
      ClearCache();
      start_time = get_time();
      for (k=0; k < nvecs; k++) {
        N_VScale(c[0], X[0][k], Z[k]);
        for (j=1; j < nsums; j++) {
          N_VLinearSum(ONE, Z[k], c[j], X[j][k], Z[k]);
        }
      }
      sync_device(V);
      stop_time = get_time();
    }

    utimes[i] = stop_time - start_time;
  }

  /* fill data */
  for (j=0; j < nsums; j++)
    c[j] = ((realtype)rand() / (realtype)RAND_MAX) + 1.0;

  for (j=0; j < nsums; j++)
    for (k=0; k < nvecs; k++)
      N_VRand(X[j][k], local_length, NEG_ONE, ONE);

  /* fused operation */
  for (i=0; i < ntests+nextra; i++) {

    ClearCache();
    start_time = get_time();
    ier = N_VLinearCombinationVectorArray(nvecs, nsums, c, X, Z);
    sync_device(V);
    stop_time = get_time();

    ftimes[i] = stop_time - start_time;
  }

  /* get average time ignoring the first nextra tests */
  time_stats(V, ftimes, nextra, ntests, &favgtime, &fsdevtime, &fmintime, &fmaxtime);
  time_stats(V, utimes, nextra, ntests, &uavgtime, &usdevtime, &umintime, &umaxtime);
  PRINT_TIME2("N_VLinearCombinationVectorArray-3",
              favgtime, fsdevtime, fmintime, fmaxtime,
              uavgtime, usdevtime, umintime, umaxtime);

  /* Free vectors */
  free(ftimes);
  free(utimes);
  free(c);
  N_VDestroyVectorArray(Z, nvecs);

  for (j=0; j < nsums; j++)
    N_VDestroyVectorArray(X[j], nvecs);
  free(X);

  return(ier);
}

/* ======================================================================
 * Exported utility functions
 * ====================================================================*/

/* ----------------------------------------------------------------------
 * Print table headers for test output
 * --------------------------------------------------------------------*/
void PrintTableHeader(int type)
{
  switch(type) {

  case 1:
    printf("\n%33s %22s %22s %22s %22s\n",
           "Operation","Avg","Std Dev","Min","Max");
    break;

  case 2:
    printf("\n%33s %22s %22s %22s %22s %22s %22s %22s %22s\n",
           "Operation",
           "Avg Fused","Std Dev Fused","Min Fused","Max Fused",
           "Avg Unfused","Std Dev Unfused","Min Unfused","Max Unfused");
    break;
  }
}

void SetTiming(int onoff, int myid)
{
#if defined(SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime( CLOCK_MONOTONIC_RAW, &spec );
  base_time_tv_sec = spec.tv_sec;

  clock_getres( CLOCK_MONOTONIC_RAW, &spec);
  if (myid == 0)
    printf("Timer resolution: %ld ns = %g s\n", spec.tv_nsec, ((double)(spec.tv_nsec) / 1E9));
#endif

  /* only print from the root process */
  if (myid == 0)
    print_time = onoff;
  else
    print_time = 0;
}


/* ----------------------------------------------------------------------
 * Fill a realtype array with random numbers between lower and upper
 * using a linear congruential generator suggested in the C99 standard
 * --------------------------------------------------------------------*/
void rand_realtype(realtype *data, sunindextype len, realtype lower, realtype upper)
{
  int i, rand;
  realtype range;
  static int rand_max = 0x7fffffff; /* 2^32 - 1 */

  /* fill array with random data between lower and upper */
  range = upper - lower;
  rand  = (int) time(NULL);
  for (i=0; i < len; i++) {
    rand = (1103515245*rand+12345) & rand_max;
    data[i] = range*((realtype)rand / (realtype)rand_max) + lower;
  }

  return;
}


/* ----------------------------------------------------------------------
 * Fill a realtype array with random series of 0 or 1
 * --------------------------------------------------------------------*/
void rand_realtype_zero_one(realtype *data, sunindextype len)
{
  int i, rand;
  static int rand_max = 0x7fffffff; /* 2^32 - 1 */

  /* fill vector randomly with 0 or 1 */
  rand = (int) time(NULL);
  for (i=0; i < len; i++) {
    rand = (1103515245*rand+12345) & rand_max;
    data[i] = (realtype) (rand % 2);
  }

  return;
}


/* ----------------------------------------------------------------------
 * Fill a realtype array with random values for constraint testing
 * --------------------------------------------------------------------*/
void rand_realtype_constraints(realtype *data, sunindextype len)
{
  int i, rand;
  static int rand_max = 0x7fffffff; /* 2^32 - 1 */

  /* randomly fill vector with the values -2, -1, 0, 1, 2 */
  rand = (int) time(NULL);
  for (i=0; i < len; i++) {
    rand = (1103515245*rand+12345) & rand_max;
    data[i] = (realtype) (rand % 5 - 2);
  }

  return;
}


/* ======================================================================
 * Private functions
 * ====================================================================*/

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
static double get_time()
{
  double time;
#if defined(SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime( CLOCK_MONOTONIC_RAW, &spec );
  time = (double)(spec.tv_sec - base_time_tv_sec) + ((double)(spec.tv_nsec) / 1E9);
#else
  time = 0;
#endif
  return time;
}


/* ----------------------------------------------------------------------
 * compute average, standard deviation, max, and min
 * --------------------------------------------------------------------*/
static void time_stats(N_Vector X, double *times, int nextra, int ntests,
                       double *avg, double *sdev, double *min, double *max)
{
  int i, ntotal;

  /* total number of times collected */
  ntotal = nextra+ntests;

  /* if running in parallel collect data from all processes */
  collect_times(X, times, ntotal);

  /* compute timing stats */
  *avg = 0.0;
  *min = times[nextra];
  *max = times[nextra];

  for (i=nextra; i<ntotal; i++) {
    *avg += times[i];
    if (times[i] < *min) *min = times[i];
    if (times[i] > *max) *max = times[i];
  }
  *avg /= ntests;

  *sdev = 0.0;
  if (ntests > 1) {
    for (i=nextra; i<ntotal; i++)
      *sdev += (times[i] - *avg) * (times[i] - *avg);
    *sdev = sqrt(*sdev/(ntests-1));
  }
}
