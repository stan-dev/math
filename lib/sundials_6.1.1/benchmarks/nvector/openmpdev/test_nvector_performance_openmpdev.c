/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the testing routine to evaluate the performance of the
 * OpenMPDEV NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_openmpdev.h>
#include <sundials/sundials_math.h>
#include "test_nvector_performance.h"

/* private functions */
static int InitializeClearCache(int cachesize);
static int FinalizeClearCache();

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  N_Vector X;          /* test vector       */
  sunindextype veclen; /* vector length     */

  int print_timing;    /* output timings    */
  int ntests;          /* number of tests   */
  int nvecs;           /* number of tests   */
  int nsums;           /* number of sums    */
  int cachesize;       /* size of cache (MB) */
  int flag;            /* return flag       */

  printf("\nStart Tests\n");
  printf("Vector Name: OpenMPDEV\n");

  /* check input and set vector length */
  if (argc < 7){
    printf("ERROR: SIX (6) arguments required: ");
    printf("<vector length> <number of vectors> <number of sums> <number of tests> ");
    printf("<cachesize> <print timing>\n");
    return(-1);
  }

  veclen = atol(argv[1]);
  if (veclen <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    return(-1);
  }

  nvecs = atol(argv[2]);
  if (nvecs <= 0) {
    printf("ERROR: number of vectors must be a positive integer \n");
    return(-1);
  }

  nsums = atol(argv[3]);
  if (nsums <= 0) {
    printf("ERROR: number of sums must be a positive integer \n");
    return(-1);
  }

  ntests = atol(argv[4]);
  if (ntests <= 0) {
    printf("ERROR: number of tests must be a positive integer \n");
    return(-1);
  }

  cachesize = atol(argv[5]);
  if (cachesize < 0) {
    printf("ERROR: cache size (MB) must be a non-negative integer \n");
    return(-1);
  }
  InitializeClearCache(cachesize);

  print_timing = atoi(argv[6]);
  SetTiming(print_timing, 0);


  printf("\nRunning with: \n");
  printf("  vector length         %ld \n", (long int) veclen);
  printf("  max number of vectors %d  \n", nvecs);
  printf("  max number of sums    %d  \n", nsums);
  printf("  number of tests       %d  \n", ntests);
  printf("  timing on/off         %d  \n", print_timing);

  /* Create vectors */
  X = N_VNew_OpenMPDEV(veclen);

  /* run tests */
  if (print_timing) printf("\n\n standard operations:\n");
  if (print_timing) PrintTableHeader(1);
  flag = Test_N_VLinearSum(X, veclen, ntests);
  flag = Test_N_VConst(X, veclen, ntests);
  flag = Test_N_VProd(X, veclen, ntests);
  flag = Test_N_VDiv(X, veclen, ntests);
  flag = Test_N_VScale(X, veclen, ntests);
  flag = Test_N_VAbs(X, veclen, ntests);
  flag = Test_N_VInv(X, veclen, ntests);
  flag = Test_N_VAddConst(X, veclen, ntests);
  flag = Test_N_VDotProd(X, veclen, ntests);
  flag = Test_N_VMaxNorm(X, veclen, ntests);
  flag = Test_N_VWrmsNorm(X, veclen, ntests);
  flag = Test_N_VWrmsNormMask(X, veclen, ntests);
  flag = Test_N_VMin(X, veclen, ntests);
  flag = Test_N_VWL2Norm(X, veclen, ntests);
  flag = Test_N_VL1Norm(X, veclen, ntests);
  flag = Test_N_VCompare(X, veclen, ntests);
  flag = Test_N_VInvTest(X, veclen, ntests);
  flag = Test_N_VConstrMask(X, veclen, ntests);
  flag = Test_N_VMinQuotient(X, veclen, ntests);

  if (print_timing) printf("\n\n fused operations 1: nvecs= %d\n", nvecs);
  if (print_timing) PrintTableHeader(2);
  flag = Test_N_VLinearCombination(X, veclen, nvecs, ntests);
  flag = Test_N_VScaleAddMulti(X, veclen, nvecs, ntests);
  flag = Test_N_VDotProdMulti(X, veclen, nvecs, ntests);
  flag = Test_N_VLinearSumVectorArray(X, veclen, nvecs, ntests);
  flag = Test_N_VScaleVectorArray(X, veclen, nvecs, ntests);
  flag = Test_N_VConstVectorArray(X, veclen, nvecs, ntests);
  flag = Test_N_VWrmsNormVectorArray(X, veclen, nvecs, ntests);
  flag = Test_N_VWrmsNormMaskVectorArray(X, veclen, nvecs, ntests);

  if (print_timing) printf("\n\n fused operations 2: nvecs= %d nsums= %d\n", nvecs, nsums);
  if (print_timing) PrintTableHeader(2);
  flag = Test_N_VScaleAddMultiVectorArray(X, veclen, nvecs, nsums, ntests);
  flag = Test_N_VLinearCombinationVectorArray(X, veclen, nvecs, nsums, ntests);

  /* Free vectors */
  N_VDestroy(X);

  FinalizeClearCache();

  printf("\nFinished Tests\n");

  return(flag);
}


/* ----------------------------------------------------------------------
 * Functions required by testing routines to fill vector data
 * --------------------------------------------------------------------*/

/* random data between lower and upper */
void N_VRand(N_Vector Xvec, sunindextype Xlen, realtype lower, realtype upper)
{
  realtype *Xdata;

  Xdata = N_VGetHostArrayPointer_OpenMPDEV(Xvec);
  rand_realtype(Xdata, Xlen, lower, upper);
  N_VCopyToDevice_OpenMPDEV(Xvec);
}

/* series of 0 and 1 */
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen)
{
  realtype *Xdata;

  Xdata = N_VGetHostArrayPointer_OpenMPDEV(Xvec);
  rand_realtype_zero_one(Xdata, Xlen);
  N_VCopyToDevice_OpenMPDEV(Xvec);
}

/* random values for constraint array */
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen)
{
  realtype *Xdata;

  Xdata = N_VGetHostArrayPointer_OpenMPDEV(Xvec);
  rand_realtype_constraints(Xdata, Xlen);
  N_VCopyToDevice_OpenMPDEV(Xvec);
}


/* ----------------------------------------------------------------------
 * Functions required for MPI or GPU testing
 * --------------------------------------------------------------------*/

void collect_times(N_Vector X, double *times, int ntimes)
{
  /* not running with MPI, just return */
  return;
}

void sync_device(N_Vector x)
{
  /* not running on GPU, just return */
  return;
}


/* ----------------------------------------------------------------------
 * Functions required for clearing cache
 * --------------------------------------------------------------------*/

static int InitializeClearCache(int cachesize)
{
  return(0);
}

static int FinalizeClearCache()
{
  return(0);
}

void ClearCache()
{
  return;
}
