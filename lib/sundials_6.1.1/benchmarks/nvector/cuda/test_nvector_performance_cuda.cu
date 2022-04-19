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
 * This is the testing routine to check the performance of the
 * NVECTOR CUDA module implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_cuda.h>
#include <sundials/sundials_math.h>
#include "test_nvector_performance.h"

/* private functions */
static int InitializeClearCache(int cachesize);
static int FinalizeClearCache();

/* private data for clearing cache */
static sunindextype N;    /* data length */
static realtype* h_data;  /* host data   */
static realtype* h_sum;   /* host sum    */
static realtype* d_data;  /* device data */
static realtype* d_sum;   /* device sum  */
static int blocksPerGrid;

/* cuda reduction kernel to clearing cache between tests */
__global__
void ClearCacheKernel(sunindextype N, realtype* data, realtype* out)
{
  __shared__ realtype shared[256];

  int sharedidx = blockIdx.x;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  realtype tmp = 0;
  while (tid < N) {
    tmp += data[tid];
    tid += blockDim.x * gridDim.x;
  }
  shared[sharedidx] = tmp;
  __syncthreads();

  /* assues blockDim is a power of 2 */
  int i = blockDim.x/2;
  while (i != 0) {
    if (sharedidx < i)
      shared[sharedidx] += shared[sharedidx + i];
    __syncthreads();
    i /= 2;
  }

  if (sharedidx == 0)
    out[sharedidx] = shared[0];
}

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

  printf("\nStart Tests\n");
  printf("Vector Name: Cuda\n");

  /* check input and set vector length */
  if (argc < 7){
    printf("ERROR: SIX (6) arguments required: ");
    printf("<vector length> <number of vectors> <number of sums> <number of tests> ");
    printf("<cache size (MB)> <print timing>\n");
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

  flag = SUNContext_Create(NULL, &ctx);
  if (flag) return flag;

  /* Create vectors */
  X = N_VNew_Cuda(veclen, ctx);

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

  flag = SUNContext_Free(&ctx);
  if (flag) return flag;

  printf("\nFinished Tests\n");

  return(flag);
}


/* ----------------------------------------------------------------------
 * Functions required by testing routines to fill vector data
 * --------------------------------------------------------------------*/

/* random data between lower and upper */
void N_VRand(N_Vector Xvec, sunindextype Xlen, realtype lower, realtype upper)
{
  rand_realtype(N_VGetHostArrayPointer_Cuda(Xvec), Xlen, lower, upper);
  N_VCopyToDevice_Cuda(Xvec);
}

/* series of 0 and 1 */
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen)
{
  rand_realtype_zero_one(N_VGetHostArrayPointer_Cuda(Xvec), Xlen);
  N_VCopyToDevice_Cuda(Xvec);
}

/* random values for constraint array */
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen)
{
  rand_realtype_constraints(N_VGetHostArrayPointer_Cuda(Xvec), Xlen);
  N_VCopyToDevice_Cuda(Xvec);
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
  cudaDeviceSynchronize();
  return;
}


/* ----------------------------------------------------------------------
 * Functions required for clearing cache
 * --------------------------------------------------------------------*/

static int InitializeClearCache(int cachesize)
{
  cudaError_t err;     /* cuda error flag     */
  size_t      nbytes;  /* cache size in bytes */

  /* determine size of vector to clear cache, N = ceil(2 * nbytes/realtype) */
  nbytes = (size_t) (2 * cachesize * 1024 * 1024);
  N = (sunindextype) ((nbytes + sizeof(realtype) - 1)/sizeof(realtype));

  /* allocate host data */
  blocksPerGrid = SUNMIN(32,(N+255)/256);

  h_data = (realtype*) malloc(N*sizeof(realtype));
  h_sum  = (realtype*) malloc(blocksPerGrid*sizeof(realtype));

  /* allocate device data */
  err = cudaMalloc((void**) &d_data, N*sizeof(realtype));
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to allocate device vector (error code %d )!\n",err);
    return(-1);
  }

  err = cudaMalloc((void**) &d_sum, blocksPerGrid*sizeof(realtype));
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to allocate device vector (error code %d )!\n",err);
    return(-1);
  }

  /* fill host vector with random data and copy to device */
  rand_realtype(h_data, N, RCONST(-1.0), RCONST(1.0));

  err = cudaMemcpy(d_data, h_data, N*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to copy data from host to device (error code %d )!\n",err);
    return(-1);
  }

  return(0);
}

static int FinalizeClearCache()
{
  cudaError_t err;  /* cuda error flag */

  free(h_data);
  free(h_sum);

  err = cudaFree(d_data);
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to free device data (error code %d )!\n",err);
    return(-1);
  }

  err = cudaFree(d_sum);
  if (err != cudaSuccess) {
    fprintf(stderr,"Failed to free device data (error code %d )!\n",err);
    return(-1);
  }

  return(0);
}

void ClearCache()
{
  /* call cuda kernel to clear the cache */
  ClearCacheKernel<<<SUNMIN(32,(N+255)/256), 256>>>(N, d_data, d_sum);
  cudaMemcpy(h_sum, d_sum, blocksPerGrid*sizeof(realtype), cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  return;
}
