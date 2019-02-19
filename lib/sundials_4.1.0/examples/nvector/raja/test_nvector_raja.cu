/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the testing routine to check the NVECTOR Raja module
 * implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_raja.h>
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
  N_Vector     U, V, X, Y, Z;     /* test vectors              */
  int          print_timing;      /* turn timing on/off        */


  /* check input and set vector length */
  if (argc < 3){
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    return(-1);
  }

  length = atol(argv[1]);
  if (length <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    return(-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, 0);

  printf("Testing RAJA N_Vector \n");
  printf("Vector length %ld \n\n", (long int) length);

  /* Create new vectors */
  X = N_VNew_Raja(length);
  if (X == NULL) {
    printf("FAIL: Unable to create a new vector \n\n");
    return(1);
  }

  /* Check vector ID */
  fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_RAJA, 0);

  /* Test clone functions */
  fails += Test_N_VCloneEmpty(X, 0);
  fails += Test_N_VClone(X, length, 0);
  fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
  fails += Test_N_VCloneVectorArray(5, X, length, 0);

  /* Clone additional vectors for testing */
  Y = N_VClone(X);
  if (Y == NULL) {
    N_VDestroy(X);
    printf("FAIL: Unable to create a new vector \n\n");
    return(1);
  }

  Z = N_VClone(X);
  if (Z == NULL) {
    N_VDestroy(X);
    N_VDestroy(Y);
    printf("FAIL: Unable to create a new vector \n\n");
    return(1);
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
  fails += Test_N_VDotProd(X, Y, length, length, 0);
  fails += Test_N_VMaxNorm(X, length, 0);
  fails += Test_N_VWrmsNorm(X, Y, length, 0);
  fails += Test_N_VWrmsNormMask(X, Y, Z, length, length, 0);
  fails += Test_N_VMin(X, length, 0);
  fails += Test_N_VWL2Norm(X, Y, length, length, 0);
  fails += Test_N_VL1Norm(X, length, length, 0);
  fails += Test_N_VCompare(X, Z, length, 0);
  fails += Test_N_VInvTest(X, Z, length, 0);
  fails += Test_N_VConstrMask(X, Y, Z, length, 0);
  fails += Test_N_VMinQuotient(X, Y, length, 0);

  /* Fused and vector array operations tests (disabled) */
  printf("\nTesting fused and vector array operations (disabled):\n\n");

  /* create vector and disable all fused and vector array operations */
  U = N_VNew_Raja(length);
  retval = N_VEnableFusedOps_Raja(U, SUNFALSE);
  if (U == NULL || retval != 0) {
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    printf("FAIL: Unable to create a new vector \n\n");
    return(1);
  }

  /* fused operations */
  fails += Test_N_VLinearCombination(U, length, 0);
  fails += Test_N_VScaleAddMulti(U, length, 0);
  fails += Test_N_VDotProdMulti(U, length, length, 0);

  /* vector array operations */
  fails += Test_N_VLinearSumVectorArray(U, length, 0);
  fails += Test_N_VScaleVectorArray(U, length, 0);
  fails += Test_N_VConstVectorArray(U, length, 0);
  fails += Test_N_VWrmsNormVectorArray(U, length, 0);
  fails += Test_N_VWrmsNormMaskVectorArray(U, length, length, 0);
  fails += Test_N_VScaleAddMultiVectorArray(U, length, 0);
  fails += Test_N_VLinearCombinationVectorArray(U, length, 0);

  /* Fused and vector array operations tests (enabled) */
  printf("\nTesting fused and vector array operations (enabled):\n\n");

  /* create vector and enable all fused and vector array operations */
  V = N_VNew_Raja(length);
  retval = N_VEnableFusedOps_Raja(V, SUNTRUE);
  if (V == NULL || retval != 0) {
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
    printf("FAIL: Unable to create a new vector \n\n");
    return(1);
  }

  /* fused operations */
  fails += Test_N_VLinearCombination(V, length, 0);
  fails += Test_N_VScaleAddMulti(V, length, 0);
  fails += Test_N_VDotProdMulti(V, length, length, 0);

  /* vector array operations */
  fails += Test_N_VLinearSumVectorArray(V, length, 0);
  fails += Test_N_VScaleVectorArray(V, length, 0);
  fails += Test_N_VConstVectorArray(V, length, 0);
  fails += Test_N_VWrmsNormVectorArray(V, length, 0);
  fails += Test_N_VWrmsNormMaskVectorArray(V, length, length, 0);
  fails += Test_N_VScaleAddMultiVectorArray(V, length, 0);
  fails += Test_N_VLinearCombinationVectorArray(V, length, 0);

  /* Free vectors */
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
    failure += FNEQ(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  /* check if vector content is non-null */
  return (X->content == NULL ? SUNFALSE : SUNTRUE);
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  /* set i-th element of data array */
  N_VCopyFromDevice_Raja(X);
  (N_VGetHostArrayPointer_Raja(X))[i] = val;
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

void sync_device()
{
  /* sync with GPU */
  cudaDeviceSynchronize();
  return;
}
