/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Daniel McGreer @ LLNL
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
 * This is the testing routine to check the NVECTOR Raja module
 * implementation.
 * -----------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>

#include <sundials/sundials_types.h>
#include <nvector/nvector_raja.h>
#include <sundials/sundials_math.h>

#include <RAJA/RAJA.hpp>

#if defined(SUNDIALS_RAJA_BACKENDS_CUDA) || defined(SUNDIALS_RAJA_BACKENDS_HIP)
#include "custom_memory_helper_gpu.h"
#elif defined(SUNDIALS_RAJA_BACKENDS_SYCL)
#include "custom_memory_helper_sycl.h"
#else
#error "Unknown RAJA backend"
#endif

#include "test_nvector.h"

/* Managed or unmanaged memory options */
enum mem_type { UNMANAGED, MANAGED, SUNMEMORY };

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails = 0;         /* counter for test failures */
  int          retval;            /* function return value     */
  sunindextype length;            /* vector length             */
  N_Vector     U, V, X, Y, Z;     /* test vectors              */
  int          print_timing;      /* turn timing on/off        */
  int          memtype;

  Test_Init(NULL);

  /* check input and set vector length */
  if (argc < 3){
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    Test_Abort(-1);
  }

  length = (sunindextype) atol(argv[1]);
  if (length <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    Test_Abort(-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, 0);

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  camp::resources::Resource* sycl_res = new camp::resources::Resource{camp::resources::Sycl()};
  ::RAJA::sycl::detail::setQueue(sycl_res);
  sycl::queue* myQueue = ::RAJA::sycl::detail::getQueue();
#endif

  /* test with both memory variants */
  for (memtype=UNMANAGED; memtype<=SUNMEMORY; ++memtype) {
    SUNMemoryHelper mem_helper = NULL;

    printf("=====> Beginning setup\n\n");

    if (memtype==UNMANAGED) {
      printf("Testing RAJA N_Vector\n");
    } else if (memtype==MANAGED) {
      printf("Testing RAJA N_Vector with managed memory\n");
    } else if (memtype==SUNMEMORY) {
      printf("Testing RAJA N_Vector with custom allocator\n");
      mem_helper = MyMemoryHelper(sunctx);
    }
    printf("Vector length: %ld \n", (long int) length);
    /* Create new vectors */
    if (memtype == UNMANAGED)       X = N_VNew_Raja(length, sunctx);
    else if (memtype == MANAGED)    X = N_VNewManaged_Raja(length, sunctx);
    else if (memtype == SUNMEMORY)  X = N_VNewWithMemHelp_Raja(length, SUNFALSE,
                                                               mem_helper, sunctx);
    if (X == NULL) {
      if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
      printf("FAIL: Unable to create a new vector \n\n");
      Test_Abort(1);
    }

    /* Fill vector with uniform random data in [-1,1] */
    realtype* xdata = N_VGetHostArrayPointer_Raja(X);
    for (sunindextype j=0; j<length; j++)
      xdata[j] = ((realtype) rand() / (realtype) RAND_MAX)*2-1;
    N_VCopyToDevice_Raja(X);

    /* Clone additional vectors for testing */
    Y = N_VClone(X);
    if (Y == NULL) {
      N_VDestroy(X);
      if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
      printf("FAIL: Unable to create a new vector \n\n");
      Test_Abort(1);
    }

    Z = N_VClone(X);
    if (Z == NULL) {
      N_VDestroy(X);
      N_VDestroy(Y);
      if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
      printf("FAIL: Unable to create a new vector \n\n");
      Test_Abort(1);
    }

    /* Fill vectors with uniform random data in [-1,1] */
    realtype* ydata = N_VGetHostArrayPointer_Raja(Y);
    realtype* zdata = N_VGetHostArrayPointer_Raja(Z);
    for (sunindextype j=0; j<length; j++) {
      ydata[j] = ((realtype) rand() / (realtype) RAND_MAX)*2-1;
      zdata[j] = ((realtype) rand() / (realtype) RAND_MAX)*2-1;
    }
    N_VCopyToDevice_Raja(Y);
    N_VCopyToDevice_Raja(Z);

    printf("\n=====> Setup complete\n");
    printf("=====> Beginning tests\n\n");

    /* Standard vector operation tests */
    printf("\nTesting standard vector operations:\n\n");

    /* Check vector ID */
    fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_RAJA, 0);

    /* Check vector length */
    fails += Test_N_VGetLength(X, 0);

    /* Check vector communicator */
    fails += Test_N_VGetCommunicator(X, NULL, 0);

    /* Test clone functions */
    fails += Test_N_VCloneEmpty(X, 0);
    fails += Test_N_VClone(X, length, 0);
    fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
    fails += Test_N_VCloneVectorArray(5, X, length, 0);

    /* Test vector math kernels */
    fails += Test_N_VConst(X, length, 0);
    fails += Test_N_VLinearSum(X, Y, Z, length, 0);
    fails += Test_N_VProd(X, Y, Z, length, 0);
    fails += Test_N_VDiv(X, Y, Z, length, 0);
    fails += Test_N_VScale(X, Z, length, 0);
    fails += Test_N_VAbs(X, Z, length, 0);
    fails += Test_N_VInv(X, Z, length, 0);
    fails += Test_N_VAddConst(X, Z, length, 0);
    fails += Test_N_VDotProd(X, Y, length, 0);
    fails += Test_N_VMaxNorm(X, length, 0);
    fails += Test_N_VWrmsNorm(X, Y, length, 0);
    fails += Test_N_VWrmsNormMask(X, Y, Z, length, 0);
    fails += Test_N_VMin(X, length, 0);
    fails += Test_N_VWL2Norm(X, Y, length, 0);
    fails += Test_N_VL1Norm(X, length, 0);
    fails += Test_N_VCompare(X, Z, length, 0);
    fails += Test_N_VInvTest(X, Z, length, 0);
    fails += Test_N_VConstrMask(X, Y, Z, length, 0);
    fails += Test_N_VMinQuotient(X, Y, length, 0);

    /* Fused and vector array operations tests (disabled) */
    printf("\nTesting fused and vector array operations (disabled):\n\n");

    /* create vector and disable all fused and vector array operations */
    U = N_VClone(X);
    retval = N_VEnableFusedOps_Raja(U, SUNFALSE);
    if (U == NULL || retval != 0) {
      N_VDestroy(X);
      N_VDestroy(Y);
      N_VDestroy(Z);
      if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
      printf("FAIL: Unable to create a new vector \n\n");
      Test_Abort(1);
    }

    /* fused operations */
    fails += Test_N_VLinearCombination(U, length, 0);
    fails += Test_N_VScaleAddMulti(U, length, 0);
    fails += Test_N_VDotProdMulti(U, length, 0);

    /* vector array operations */
    fails += Test_N_VLinearSumVectorArray(U, length, 0);
    fails += Test_N_VScaleVectorArray(U, length, 0);
    fails += Test_N_VConstVectorArray(U, length, 0);
    fails += Test_N_VWrmsNormVectorArray(U, length, 0);
    fails += Test_N_VWrmsNormMaskVectorArray(U, length, 0);
    fails += Test_N_VScaleAddMultiVectorArray(U, length, 0);
    fails += Test_N_VLinearCombinationVectorArray(U, length, 0);

    /* Fused and vector array operations tests (enabled) */
    printf("\nTesting fused and vector array operations (enabled):\n\n");

    /* create vector and enable all fused and vector array operations */
    V = N_VClone(X);
    retval = N_VEnableFusedOps_Raja(V, SUNTRUE);
    if (V == NULL || retval != 0) {
      N_VDestroy(X);
      N_VDestroy(Y);
      N_VDestroy(Z);
      N_VDestroy(U);
      if (mem_helper) SUNMemoryHelper_Destroy(mem_helper);
      printf("FAIL: Unable to create a new vector \n\n");
      Test_Abort(1);
    }

    /* fused operations */
    fails += Test_N_VLinearCombination(V, length, 0);
    fails += Test_N_VScaleAddMulti(V, length, 0);
    fails += Test_N_VDotProdMulti(V, length, 0);

    /* vector array operations */
    fails += Test_N_VLinearSumVectorArray(V, length, 0);
    fails += Test_N_VScaleVectorArray(V, length, 0);
    fails += Test_N_VConstVectorArray(V, length, 0);
    fails += Test_N_VWrmsNormVectorArray(V, length, 0);
    fails += Test_N_VWrmsNormMaskVectorArray(V, length, 0);
    fails += Test_N_VScaleAddMultiVectorArray(V, length, 0);
    fails += Test_N_VLinearCombinationVectorArray(V, length, 0);

    /* local reduction operations */
    printf("\nTesting local reduction operations:\n\n");

    fails += Test_N_VDotProdLocal(X, Y, length, 0);
    fails += Test_N_VMaxNormLocal(X, length, 0);
    fails += Test_N_VMinLocal(X, length, 0);
    fails += Test_N_VL1NormLocal(X, length, 0);
    fails += Test_N_VWSqrSumLocal(X, Y, length, 0);
    fails += Test_N_VWSqrSumMaskLocal(X, Y, Z, length, 0);
    fails += Test_N_VInvTestLocal(X, Z, length, 0);
    fails += Test_N_VConstrMaskLocal(X, Y, Z, length, 0);
    fails += Test_N_VMinQuotientLocal(X, Y, length, 0);

    /* XBraid interface operations */
    printf("\nTesting XBraid interface operations:\n\n");

    fails += Test_N_VBufSize(X, length, 0);
    fails += Test_N_VBufPack(X, length, 0);
    fails += Test_N_VBufUnpack(X, length, 0);

    printf("\n=====> Beginning teardown\n");

    /* Free vectors */
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
    N_VDestroy(V);
    if(mem_helper) SUNMemoryHelper_Destroy(mem_helper);

    /* Synchronize */
    sync_device(NULL);

    printf("=====> Teardown complete\n\n");
  }

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests \n\n", fails);
  } else {
    printf("SUCCESS: NVector module passed all tests \n\n");
  }

  sync_device(NULL);
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  cudaDeviceReset();
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  hipDeviceReset();
#endif
  Test_Finalize();
  return(fails);
}


/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int          failure = 0;
  sunindextype i;
  realtype     *Xdata;

  N_VCopyFromDevice_Raja(X);
  Xdata = N_VGetHostArrayPointer_Raja(X);

  /* check vector data */
  for (i = 0; i < local_length; i++) {
    failure += SUNRCompare(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  /* check if vector data is non-null */
  if ((N_VGetHostArrayPointer_Raja(X) == NULL) &&
      (N_VGetDeviceArrayPointer_Raja(X) == NULL))
    return SUNFALSE;
  return SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  /* set i-th element of data array */
  set_element_range(X, i, i, val);
}

void set_element_range(N_Vector X, sunindextype is, sunindextype ie,
                       realtype val)
{
  sunindextype i;
  realtype*    xd;

  /* set elements [is,ie] of the data array */
  N_VCopyFromDevice_Raja(X);
  xd = N_VGetHostArrayPointer_Raja(X);
  for(i = is; i <= ie; i++) xd[i] = val;
  N_VCopyToDevice_Raja(X);
}

realtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  N_VCopyFromDevice_Raja(X);
  return (N_VGetHostArrayPointer_Raja(X))[i];
}

double max_time(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return(time);
}

void sync_device(N_Vector x)
{
  /* sync with GPU */
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  cudaDeviceSynchronize();
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  hipDeviceSynchronize();
#elif defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  sycl::queue* myQueue = ::RAJA::sycl::detail::getQueue();
  myQueue->wait_and_throw();
#endif

  return;
}
