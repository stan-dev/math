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
 * This is the testing routine to check the POSIX Threads (Pthreads)
 * NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_pthreads.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails = 0;         /* counter for test failures */
  int          retval;            /* function return value     */
  sunindextype length;            /* vector length             */
  N_Vector     U, V, W, X, Y, Z;  /* test vectors              */
  int          print_timing;      /* turn timing on/off        */
  int          nthreads;          /* number of POSIX threads   */

  Test_Init(NULL);

  /* check input and set vector length */
  if (argc < 4){
    printf("ERROR: THREE (3) Inputs required: vector length, number of threads, print timing \n");
    Test_Abort(-1);
  }

  length = (sunindextype) atol(argv[1]);
  if (length <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    Test_Abort(-1);
  }

  nthreads = atoi(argv[2]);
  if (nthreads < 1) {
    printf("ERROR: number of threads must be at least 1 \n");
    Test_Abort(-1);
  }

  print_timing = atoi(argv[3]);
  SetTiming(print_timing, 0);

  printf("Testing the Pthreads N_Vector \n");
  printf("Vector length %ld \n", (long int) length);
  printf("Number of threads %d \n\n", nthreads);

  /* Create new vectors */
  W = N_VNewEmpty_Pthreads(length, nthreads, sunctx);
  if (W == NULL) {
    printf("FAIL: Unable to create a new empty vector \n\n");
    Test_Abort(1);
  }

  X = N_VNew_Pthreads(length, nthreads, sunctx);
  if (X == NULL) {
    N_VDestroy(W);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Abort(1);
  }

  /* Check vector ID */
  fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_PTHREADS, 0);

  /* Check vector length */
  fails += Test_N_VGetLength(X, 0);

  /* Check vector communicator */
  fails += Test_N_VGetCommunicator(X, NULL, 0);

  /* Test clone functions */
  fails += Test_N_VCloneEmpty(X, 0);
  fails += Test_N_VClone(X, length, 0);
  fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
  fails += Test_N_VCloneVectorArray(5, X, length, 0);

  /* Test setting/getting array data */
  fails += Test_N_VSetArrayPointer(W, length, 0);
  fails += Test_N_VGetArrayPointer(X, length, 0);

  /* Clone additional vectors for testing */
  Y = N_VClone(X);
  if (Y == NULL) {
    N_VDestroy(W);
    N_VDestroy(X);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Abort(1);
  }

  Z = N_VClone(X);
  if (Z == NULL) {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    printf("FAIL: Unable to create a new vector \n\n");
    Test_Abort(1);
  }

  /* Standard vector operation tests */
  printf("\nTesting standard vector operations:\n\n");

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
  U = N_VNew_Pthreads(length, nthreads, sunctx);
  retval = N_VEnableFusedOps_Pthreads(U, SUNFALSE);
  if (U == NULL || retval != 0) {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
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
  V = N_VNew_Pthreads(length, nthreads, sunctx);
  retval = N_VEnableFusedOps_Pthreads(V, SUNTRUE);
  if (V == NULL || retval != 0) {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
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

  /* local fused reduction operations */
  printf("\nTesting local fused reduction operations:\n\n");
  fails += Test_N_VDotProdMultiLocal(V, length, 0);

  /* XBraid interface operations */
  printf("\nTesting XBraid interface operations:\n\n");

  fails += Test_N_VBufSize(X, length, 0);
  fails += Test_N_VBufPack(X, length, 0);
  fails += Test_N_VBufUnpack(X, length, 0);

  /* Free vectors */
  N_VDestroy(W);
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);
  N_VDestroy(U);
  N_VDestroy(V);

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests \n\n", fails);
  } else {
    printf("SUCCESS: NVector module passed all tests \n\n");
  }

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

  Xdata = N_VGetArrayPointer(X);

  /* check vector data */
  for (i = 0; i < local_length; i++) {
    failure += SUNRCompare(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  /* check if data array is non-null */
  return (N_VGetArrayPointer(X) == NULL) ? SUNFALSE : SUNTRUE;
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

  /* set elements [is,ie] of the data array */
  realtype* xd = N_VGetArrayPointer(X);
  for(i = is; i <= ie; i++) xd[i] = val;
}

realtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  return NV_Ith_PT(X,i);
}

double max_time(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return(time);
}

void sync_device(N_Vector x)
{
  /* not running on GPU, just return */
  return;
}
