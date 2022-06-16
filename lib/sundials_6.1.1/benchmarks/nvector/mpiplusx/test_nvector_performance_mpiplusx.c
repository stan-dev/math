/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * MPIPlusX NVECTOR module implementation. It uses the serial vector
 * as the local vector.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_mpiplusx.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_nvector_performance.h"

#include <mpi.h>

/* private functions */
static int InitializeClearCache(int cachesize);
static int FinalizeClearCache();

/* private data for clearing cache */
static sunindextype N;  /* data length */
static realtype* data;  /* host data   */


/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  SUNContext   ctx = NULL;  /* SUNDIALS context */
  N_Vector     X   = NULL;  /* test vector      */
  sunindextype veclen;      /* vector length    */

  int print_timing;    /* output timings     */
  int ntests;          /* number of tests    */
  int nvecs;           /* number of tests    */
  int nsums;           /* number of sums     */
  int cachesize;       /* size of cache (MB) */
  int flag;            /* return flag        */

  MPI_Comm comm;          /* MPI Communicator   */
  int      nprocs, myid;  /* Num procs, proc id */

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  if (myid == 0) {
    printf("\nStart Tests\n");
    printf("Vector Name: MPIPlusX\n");
  }

  /* check input and set vector length */
  if (argc < 7){
    printf("ERROR: SIX (6) arguments required: ");
    printf("<vector length> <number of vectors> <number of sums> <number of tests> ");
    printf("<cache size (MB)> <print timing>\n");
    return(-1);
  }

  veclen = (sunindextype) atol(argv[1]);
  if (veclen <= 0) {
    printf("ERROR: local vector length must be a positive integer \n");
    return(-1);
  }

  nvecs = (int) atol(argv[2]);
  if (nvecs <= 0) {
    printf("ERROR: number of vectors must be a positive integer \n");
    return(-1);
  }

  nsums = (int) atol(argv[3]);
  if (nsums <= 0) {
    printf("ERROR: number of sums must be a positive integer \n");
    return(-1);
  }

  ntests = (int) atol(argv[4]);
  if (ntests <= 0) {
    printf("ERROR: number of tests must be a positive integer \n");
    return(-1);
  }

  cachesize = (int) atol(argv[5]);
  if (cachesize < 0) {
    printf("ERROR: cache size (MB) must be a non-negative integer \n");
    return(-1);
  }
  InitializeClearCache(cachesize);

  print_timing = atoi(argv[6]);
  SetTiming(print_timing, myid);

  if (myid == 0) {
    printf("\nRunning with: \n");
    printf("  local vector length   %ld \n", (long int) veclen);
    printf("  max number of vectors %d  \n", nvecs);
    printf("  max number of sums    %d  \n", nsums);
    printf("  number of tests       %d  \n", ntests);
    printf("  timing on/off         %d  \n", print_timing);
    printf("  number of MPI procs   %d  \n", nprocs);
  }

  flag = SUNContext_Create(&comm, &ctx);
  if (flag) return flag;

  /* Create vectors */
  X = N_VMake_MPIPlusX(comm, N_VNew_Serial(veclen, ctx), ctx);

  /* run tests */
  if (myid == 0 && print_timing) {
    printf("\n\n standard operations:\n");
    PrintTableHeader(1);
  }
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

  if (myid == 0 && print_timing) {
    printf("\n\n fused operations 1: nvecs= %d\n", nvecs);
    PrintTableHeader(2);
  }
  flag = Test_N_VLinearCombination(X, veclen, nvecs, ntests);
  flag = Test_N_VScaleAddMulti(X, veclen, nvecs, ntests);
  flag = Test_N_VDotProdMulti(X, veclen, nvecs, ntests);
  flag = Test_N_VLinearSumVectorArray(X, veclen, nvecs, ntests);
  flag = Test_N_VScaleVectorArray(X, veclen, nvecs, ntests);
  flag = Test_N_VConstVectorArray(X, veclen, nvecs, ntests);
  flag = Test_N_VWrmsNormVectorArray(X, veclen, nvecs, ntests);
  flag = Test_N_VWrmsNormMaskVectorArray(X, veclen, nvecs, ntests);

  if (myid == 0 && print_timing) {
    printf("\n\n fused operations 2: nvecs= %d nsums= %d\n", nvecs, nsums);
    PrintTableHeader(2);
  }
  flag = Test_N_VScaleAddMultiVectorArray(X, veclen, nvecs, nsums, ntests);
  flag = Test_N_VLinearCombinationVectorArray(X, veclen, nvecs, nsums, ntests);

  /* Free vectors */
  N_VDestroy(N_VGetLocalVector_MPIPlusX(X));
  N_VDestroy(X);

  FinalizeClearCache();

  flag = SUNContext_Free(&ctx);
  if (flag) return flag;

  if (myid == 0)
    printf("\nFinished Tests\n");

  MPI_Finalize();

  return(flag);
}


/* ----------------------------------------------------------------------
 * Functions required by testing routines to fill vector data
 * --------------------------------------------------------------------*/

/* random data between lower and upper */
void N_VRand(N_Vector Xvec, sunindextype Xlen, realtype lower, realtype upper)
{
  realtype *Xdata;

  Xdata = N_VGetArrayPointer(Xvec);
  rand_realtype(Xdata, Xlen, lower, upper);
}

/* series of 0 and 1 */
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen)
{
  realtype *Xdata;

  Xdata = N_VGetArrayPointer(Xvec);
  rand_realtype_zero_one(Xdata, Xlen);
}

/* random values for constraint array */
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen)
{
  realtype *Xdata;

  Xdata = N_VGetArrayPointer(Xvec);
  rand_realtype_constraints(Xdata, Xlen);
}


/* ----------------------------------------------------------------------
 * Functions required for MPI or GPU testing
 * --------------------------------------------------------------------*/

void collect_times(N_Vector X, double *times, int ntimes)
{
  MPI_Comm comm;
  int myid;

  comm = *((MPI_Comm*) N_VGetCommunicator(X));
  MPI_Comm_rank(comm, &myid);

  if (myid == 0)
    MPI_Reduce(MPI_IN_PLACE, times, ntimes, MPI_DOUBLE, MPI_MAX, 0, comm);
  else
    MPI_Reduce(times, times, ntimes, MPI_DOUBLE, MPI_MAX, 0, comm);

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
  size_t nbytes;  /* cache size in bytes */

  /* determine size of vector to clear cache, N = ceil(2 * nbytes/realtype) */
  nbytes = (size_t) (2 * cachesize * 1024 * 1024);
  N = (sunindextype) ((nbytes + sizeof(realtype) - 1)/sizeof(realtype));

  /* allocate data and fill random values */
  data = (realtype*) malloc(N*sizeof(realtype));
  rand_realtype(data, N, RCONST(-1.0), RCONST(1.0));

  return(0);
}

static int FinalizeClearCache()
{
  free(data);
  return(0);
}

void ClearCache()
{
  realtype     sum;
  sunindextype i;

  sum = RCONST(0.0);
  for (i=0; i<N; i++)
    sum += data[i];

  return;
}
