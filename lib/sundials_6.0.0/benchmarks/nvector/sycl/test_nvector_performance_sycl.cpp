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
 * This is the testing routine to check the performance of the
 * NVECTOR SYCL module implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_sycl.h>
#include <sundials/sundials_math.h>
#include "test_nvector_performance.h"

// private functions
static int InitializeClearCache(int cachesize);
static int FinalizeClearCache();

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  sundials::Context sunctx;

  // return flag
  int flag;

  printf("\nStart Tests\n");
  printf("Vector Name: Sycl\n");

  // check input and set vector length
  if (argc < 6){
    printf("ERROR: SIX (6) arguments required: ");
    printf("<vector length> <number of vectors> <number of sums> <number of tests> ");
    printf("<cachesize> <print timing>\n");
    return(-1);
  }

  // vector length
  sunindextype veclen = atol(argv[1]);
  if (veclen <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    return(-1);
  }

  // number of vectors in fused op tests
  int nvecs = atol(argv[2]);
  if (nvecs < 1) {
    printf("WARNING: Fused operation tests disabled\n");
  }

  // number of sums in fused op tests
  int nsums = atol(argv[3]);
  if (nsums < 1) {
    printf("WARNING: Some fused operation tests disabled\n");
  }

  // number of tests
  int ntests = atol(argv[4]);
  if (ntests <= 0) {
    printf("ERROR: number of tests must be a positive integer \n");
    return(-1);
  }

  // cache size
  int cachesize = atol(argv[5]);
  if (cachesize < 0) {
    printf("ERROR: cache size (MB) must be a non-negative integer \n");
    return(-1);
  }
  InitializeClearCache(cachesize);

  // output test timings
  int print_timing = atoi(argv[6]);
  SetTiming(print_timing, 0);

  printf("\nRunning with: \n");
  printf("  vector length         %ld \n", (long int) veclen);
  printf("  max number of vectors %d  \n", nvecs);
  printf("  max number of sums    %d  \n", nsums);
  printf("  number of tests       %d  \n", ntests);
  printf("  timing on/off         %d  \n", print_timing);

  /* Create an in-order GPU queue */
  ::sycl::gpu_selector selector;
  ::sycl::queue myQueue(selector,
                        ::sycl::property_list{::sycl::property::queue::in_order{}});

  ::sycl::device dev = myQueue.get_device();
  std::cout << "Running on "
            << (dev.get_info<::sycl::info::device::name>())
            << std::endl;
  std::cout << " is host? "
            << (dev.is_host() ? "Yes" : "No")
            << std::endl;
  std::cout << " is cpu? "
            << (dev.is_cpu() ? "Yes" : "No")
            << std::endl;
  std::cout << " is gpu? "
            << (dev.is_gpu() ? "Yes" : "No")
            << std::endl;
  std::cout << " is accelerator? "
            << (dev.is_accelerator() ? "Yes" : "No")
            << std::endl;
  std::cout << " is the queue in order? "
            << (myQueue.is_in_order() ? "Yes" : "No")
            << std::endl;
  std::cout << " supports usm host allocations? "
            << (dev.get_info<::sycl::info::device::usm_host_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " supports usm device allocations? "
            << (dev.get_info<::sycl::info::device::usm_device_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " suports usm shared allocations? "
            << (dev.get_info<::sycl::info::device::usm_shared_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " max work group size: "
            << dev.get_info<::sycl::info::device::max_work_group_size>()
            << std::endl;
  std::cout << " max global memory size (bytes): "
            << dev.get_info<::sycl::info::device::global_mem_size>()
            << std::endl;
  std::cout << " max local memory size (bytes): "
            << dev.get_info<::sycl::info::device::local_mem_size>()
            << std::endl;
  std::cout << std::endl;

  // Create vectors
  N_Vector X = N_VNew_Sycl(veclen, &myQueue, sunctx);

  // run tests
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

  if (nvecs > 0)
  {
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

    if (nsums > 0)
    {
      if (print_timing) printf("\n\n fused operations 2: nvecs= %d nsums= %d\n", nvecs, nsums);
      if (print_timing) PrintTableHeader(2);
      flag = Test_N_VScaleAddMultiVectorArray(X, veclen, nvecs, nsums, ntests);
      flag = Test_N_VLinearCombinationVectorArray(X, veclen, nvecs, nsums, ntests);
    }
  }

  // Sync with device
  sync_device(X);

  // Free vectors
  N_VDestroy(X);

  FinalizeClearCache();

  printf("\nFinished Tests\n");

  return(flag);
}


/* ----------------------------------------------------------------------
 * Functions required by testing routines to fill vector data
 * --------------------------------------------------------------------*/

// random data between lower and upper
void N_VRand(N_Vector Xvec, sunindextype Xlen, realtype lower, realtype upper)
{
  rand_realtype(N_VGetHostArrayPointer_Sycl(Xvec), Xlen, lower, upper);
  N_VCopyToDevice_Sycl(Xvec);
}

// series of 0 and 1
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen)
{
  rand_realtype_zero_one(N_VGetHostArrayPointer_Sycl(Xvec), Xlen);
  N_VCopyToDevice_Sycl(Xvec);
}

// random values for constraint array
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen)
{
  rand_realtype_constraints(N_VGetHostArrayPointer_Sycl(Xvec), Xlen);
  N_VCopyToDevice_Sycl(Xvec);
}


/* ----------------------------------------------------------------------
 * Functions required for MPI or GPU testing
 * --------------------------------------------------------------------*/

void collect_times(N_Vector X, double *times, int ntimes)
{
  // not running with MPI, just return
  return;
}

void sync_device(N_Vector x)
{
  ((N_VectorContent_Sycl)(x->content))->queue->wait();
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
