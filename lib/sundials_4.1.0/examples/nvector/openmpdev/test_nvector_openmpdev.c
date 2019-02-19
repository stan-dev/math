/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the testing routine to check the OpenMP 4.5 NVECTOR
 * module implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_openmpdev.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <omp.h>

/* OpenMPDEV vector specific tests */
int Test_N_VMake_OpenMPDEV(N_Vector X, sunindextype length, int myid);

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

  /* check input and set vector length */
  if (argc < 3){
    printf("ERROR: TWO (2) Inputs required: vector length and print timing \n");
    return(-1);
  }

  length = atol(argv[1]);
  if (length <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    return(-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, 0);

  printf("Testing the OpenMP DEV N_Vector \n");
  printf("Vector length %ld \n", (long int) length);
  printf("\n omp_get_default_device = %d \n", omp_get_default_device());
  printf("\n omp_get_num_devices    = %d \n", omp_get_num_devices());
  printf("\n omp_get_initial_device = %d \n", omp_get_initial_device());
  printf("\n omp_is_initial_device  = %d \n", omp_is_initial_device());

  /* Create new vectors */
  W = N_VNewEmpty_OpenMPDEV(length);
  if (W == NULL) {
    printf("FAIL: Unable to create a new empty vector \n\n");
    return(1);
  }

  X = N_VNew_OpenMPDEV(length);
  if (X == NULL) {
    N_VDestroy(W);
    printf("FAIL: Unable to create a new vector \n\n");
    return(1);
  }

  /* Check vector ID */
  fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_OPENMPDEV, 0);

  /* Test clone functions */
  fails += Test_N_VCloneEmpty(X, 0);
  fails += Test_N_VClone(X, length, 0);
  fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
  fails += Test_N_VCloneVectorArray(5, X, length, 0);

  /* Clone additional vectors for testing */
  Y = N_VClone(X);
  if (Y == NULL) {
    N_VDestroy(W);
    N_VDestroy(X);
    printf("FAIL: Unable to create a new vector \n\n");
    return(1);
  }

  Z = N_VClone(X);
  if (Z == NULL) {
    N_VDestroy(W);
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
  U = N_VNew_OpenMPDEV(length);
  retval = N_VEnableFusedOps_OpenMPDEV(U, SUNFALSE);
  if (U == NULL || retval != 0) {
    N_VDestroy(W);
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
  V = N_VNew_OpenMPDEV(length);
  retval = N_VEnableFusedOps_OpenMPDEV(V, SUNTRUE);
  if (V == NULL || retval != 0) {
    N_VDestroy(W);
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
  N_VDestroy(U);
  N_VDestroy(V);
  N_VDestroy(W);
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests \n\n", fails);
  } else {
    printf("SUCCESS: NVector module passed all tests \n\n");
  }

  return(fails);
}

/* ----------------------------------------------------------------------
 * OpenMPDEV specific tests
 * --------------------------------------------------------------------*/

/* --------------------------------------------------------------------
 * Test for the CUDA N_Vector N_VMake_OpenMPDEV function. Requires N_VConst
 * to check data.
 */
int Test_N_VMake_OpenMPDEV(N_Vector X, sunindextype length, int myid)
{
  int failure = 0;
  realtype *h_data, *d_data;
  N_Vector Y;

  N_VConst(NEG_HALF, X);
  N_VCopyFromDevice_OpenMPDEV(X);

  h_data = N_VGetHostArrayPointer_OpenMPDEV(X);
  d_data = N_VGetDeviceArrayPointer_OpenMPDEV(X);

  /* Case 1: h_data and d_data are not null */
  Y = N_VMake_OpenMPDEV(length, h_data, d_data);
  if (Y == NULL) {
    printf(">>> FAILED test -- N_VMake_OpenMPDEV, Proc %d \n", myid);
    printf("    Vector is NULL \n \n");
    return(1);
  }

  if (N_VGetHostArrayPointer_OpenMPDEV(Y) == NULL) {
    printf(">>> FAILED test -- N_VMake_OpenMPDEV, Proc %d \n", myid);
    printf("    Vector host data == NULL \n \n");
    N_VDestroy(Y);
    return(1);
  }
  
  if (N_VGetDeviceArrayPointer_OpenMPDEV(Y) == NULL) {
    printf(">>> FAILED test -- N_VMake_OpenMPDEV, Proc %d \n", myid);
    printf("    Vector device data -= NULL \n \n");
    N_VDestroy(Y);
    return(1);
  }
  
  failure += check_ans(NEG_HALF, Y, length);
 
  if (failure) {
    printf(">>> FAILED test -- N_VMake_OpenMPDEV Case 1, Proc %d \n", myid);
    printf("    Failed N_VConst check \n \n");
    N_VDestroy(Y);
    return(1);
  }
  
  if (myid == 0) {
    printf("PASSED test -- N_VMake_OpenMPDEV Case 1 \n");
  }

  N_VDestroy(Y);

  /* Case 2: data is null */
  Y = N_VMake_OpenMPDEV(length, NULL, NULL);
  if (Y != NULL) {
    printf(">>> FAILED test -- N_VMake_OpenMPDEV Case 2, Proc %d \n", myid);
    printf("    Vector is not NULL \n \n");
    return(1);
  }

  if (myid == 0) {
    printf("PASSED test -- N_VMake_OpenMPDEV Case 2 \n");
  }
  
  N_VDestroy(Y);

  return(failure);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int          failure = 0;
  sunindextype i;
  realtype     *Xdata;

  N_VCopyFromDevice_OpenMPDEV(X);
  Xdata = N_VGetHostArrayPointer_OpenMPDEV(X);

  /* check vector data */
  for (i = 0; i < local_length; i++) {
    failure += FNEQ(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  realtype *Xdata = N_VGetHostArrayPointer_OpenMPDEV(X);
  if (Xdata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  realtype *xdev;
  int dev;
  xdev = N_VGetDeviceArrayPointer_OpenMPDEV(X);
  dev = omp_get_default_device();

  #pragma omp target map(to:val) is_device_ptr(xdev) device(dev)
  {
    xdev[i] = val;
  }
}

realtype get_element(N_Vector X, sunindextype i)
{
  realtype *data;

  N_VCopyFromDevice_OpenMPDEV(X);
  data = N_VGetHostArrayPointer_OpenMPDEV(X);

  return data[i];
}

double max_time(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return(time);
}

void sync_device()
{
  /* not running on DEV, just return */
  return;
}
