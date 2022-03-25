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
 * This is the testing routine to check the NVECTOR MPIPlusX
 * module implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_mpiplusx.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails = 0;         /* counter for local test failures  */
  int          globfails = 0;     /* counter for global test failures */
  int          retval;            /* function return value            */
  sunindextype loclen;            /* overall local vector length      */
  sunindextype local_length;      /* overall local vector length      */
  sunindextype global_length;     /* overall global vector length     */
  N_Vector     Xlocal;            /* local vector                     */
  N_Vector     U, V, W, X, Y, Z;  /* test vectors                     */
  int          print_timing;      /* turn timing on/off               */
  MPI_Comm     comm;              /* shared MPI Communicator          */
  int          nprocs, myid;      /* Number of procs, proc id         */

  /* Get processor number and total number of processes */
  retval = MPI_Init(&argc, &argv);
  if (retval != MPI_SUCCESS)  return(1);

  comm = MPI_COMM_WORLD;
  Test_Init(&comm);

  retval = MPI_Comm_size(comm, &nprocs);
  if (retval != MPI_SUCCESS)  Test_AbortMPI(&comm, -1);
  retval = MPI_Comm_rank(comm, &myid);
  if (retval != MPI_SUCCESS)  Test_AbortMPI(&comm, -1);

  /* check inputs */
  if (argc < 3) {
    if (myid == 0)
      printf("ERROR: TWO (2) Inputs required: vector local length, print timing \n");
    Test_AbortMPI(&comm, -1);
  }

  local_length = (sunindextype) atol(argv[1]);
  if (local_length < 1) {
    if (myid == 0)
      printf("ERROR: local vector length must be a positive integer \n");
    Test_AbortMPI(&comm, -1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, myid);

  /* overall global length */
  global_length = nprocs*local_length;

  if (myid == 0) {
    printf("Testing the MPIPlusX N_Vector with X being a serial N_Vector\n");
    printf("Vector local length %ld\n", (long int) local_length);
    printf("MPIPlusX vector global length %ld\n", (long int) global_length);
    printf("MPI processes %d\n", nprocs);
  }

  /* Create subvectors */
  Xlocal = N_VNew_Serial(local_length, sunctx);
  if (Xlocal == NULL) {
    printf("FAIL: Unable to create a new serial vector, Proc %d\n\n", myid);
    Test_AbortMPI(&comm, 1);
  }

  /* Create a new MPIPlusX */
  X = N_VMake_MPIPlusX(comm, Xlocal, sunctx);
  if (X == NULL) {
    N_VDestroy(Xlocal);
    printf("FAIL: Unable to create a new MPIPlusX vector, Proc %d\n\n", myid);
    Test_AbortMPI(&comm, 1);
  }

  /* Check vector ID */
  if (Test_N_VGetVectorID(X, SUNDIALS_NVEC_MPIPLUSX, myid)) {
    printf(">>> FAILED test -- N_VGetVectorID, Proc %d\n\n", myid);
    fails += 1;
  }

  /* Check vector length */
  if (Test_N_VGetLength(X, myid)) {
    printf(">>> FAILED test -- N_VGetLength, Proc %d\n\n", myid);
    fails += 1;
  }

  /* Check vector communicator */
  if (Test_N_VGetCommunicatorMPI(X, &comm, myid)) {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d\n\n", myid);
    fails += 1;
  }

  /* Test local vector accessors */
  U = N_VGetLocalVector_MPIPlusX(X);
  if (N_VGetLength(U) != local_length) {
    printf(">>> FAILED test -- N_VGetLocalVector_MPIPlusX, Proc %d\n\n", myid);
    fails += 1;
  }

  loclen = N_VGetLocalLength_MPIPlusX(X);
  if (N_VGetLength(U) != loclen) {
    printf(">>> FAILED test -- N_VGetLocalLength_MPIPlusX, Proc %d\n\n", myid);
    fails += 1;
  }

  /* Clone additional vectors for testing */
  W = N_VClone(X);
  if (W == NULL) {
    N_VDestroy(X);
    N_VDestroy(Xlocal);
    printf("FAIL: Unable to create a new vector, Proc %d\n\n", myid);
    Test_AbortMPI(&comm, 1);
  }

  /* Clone additional vectors for testing */
  Y = N_VClone(X);
  if (Y == NULL) {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Xlocal);
    printf("FAIL: Unable to create a new vector, Proc %d\n\n", myid);
    Test_AbortMPI(&comm, 1);
  }

  Z = N_VClone(X);
  if (Z == NULL) {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Xlocal);
    printf("FAIL: Unable to create a new vector, Proc %d\n\n", myid);
    Test_AbortMPI(&comm, 1);
  }

  /* Standard vector operation tests */
  if (myid == 0) printf("\nTesting standard vector operations:\n\n");

  fails += Test_N_VConst(X, local_length, myid);
  fails += Test_N_VLinearSum(X, Y, Z, local_length, myid);
  fails += Test_N_VProd(X, Y, Z, local_length, myid);
  fails += Test_N_VDiv(X, Y, Z, local_length, myid);
  fails += Test_N_VScale(X, Z, local_length, myid);
  fails += Test_N_VAbs(X, Z, local_length, myid);
  fails += Test_N_VInv(X, Z, local_length, myid);
  fails += Test_N_VAddConst(X, Z, local_length, myid);
  fails += Test_N_VDotProd(X, Y, local_length, myid);
  fails += Test_N_VMaxNorm(X, local_length, myid);
  fails += Test_N_VWrmsNorm(X, Y, local_length, myid);
  fails += Test_N_VWrmsNormMask(X, Y, Z, local_length, myid);
  fails += Test_N_VMin(X, local_length, myid);
  fails += Test_N_VWL2Norm(X, Y, local_length, myid);
  fails += Test_N_VL1Norm(X, local_length, myid);
  fails += Test_N_VCompare(X, Z, local_length, myid);
  fails += Test_N_VInvTest(X, Z, local_length, myid);
  fails += Test_N_VConstrMask(X, Y, Z, local_length, myid);
  fails += Test_N_VMinQuotient(X, Y, local_length, myid);

  /* Fused and vector array operations tests (disabled) */
  if (myid == 0) printf("\nTesting fused and vector array operations (disabled):\n\n");

  /* create vector and disable all fused and vector array operations */
  U = N_VClone(X);
  retval = N_VEnableFusedOps_Serial(N_VGetLocalVector_MPIPlusX(U), SUNFALSE);
  if (U == NULL || retval != 0) {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(Xlocal);
    printf("FAIL: Unable to create a new vector, Proc %d\n\n", myid);
    Test_AbortMPI(&comm, 1);
  }

  /* fused operations */
  fails += Test_N_VLinearCombination(U, local_length, myid);
  fails += Test_N_VScaleAddMulti(U, local_length, myid);
  fails += Test_N_VDotProdMulti(U, local_length, myid);

  /* vector array operations */
  fails += Test_N_VLinearSumVectorArray(U, local_length, myid);
  fails += Test_N_VScaleVectorArray(U, local_length, myid);
  fails += Test_N_VConstVectorArray(U, local_length, myid);
  fails += Test_N_VWrmsNormVectorArray(U, local_length, myid);
  fails += Test_N_VWrmsNormMaskVectorArray(U, local_length, myid);
  fails += Test_N_VScaleAddMultiVectorArray(U, local_length, myid);
  fails += Test_N_VLinearCombinationVectorArray(U, local_length, myid);

  /* Fused and vector array operations tests (enabled) */
  if (myid == 0) printf("\nTesting fused and vector array operations (enabled):\n\n");

  /* create vector and enable all fused and vector array operations */
  V = N_VClone(X);
  retval = N_VEnableFusedOps_Serial(N_VGetLocalVector_MPIPlusX(V), SUNTRUE);
  if (V == NULL || retval != 0) {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
    N_VDestroy(Xlocal);
    printf("FAIL: Unable to create a new vector, Proc %d\n\n", myid);
    Test_AbortMPI(&comm, 1);
  }

  /* fused operations */
  fails += Test_N_VLinearCombination(V, local_length, myid);
  fails += Test_N_VScaleAddMulti(V, local_length, myid);
  fails += Test_N_VDotProdMulti(V, local_length, myid);

  /* vector array operations */
  fails += Test_N_VLinearSumVectorArray(V, local_length, myid);
  fails += Test_N_VScaleVectorArray(V, local_length, myid);
  fails += Test_N_VConstVectorArray(V, local_length, myid);
  fails += Test_N_VWrmsNormVectorArray(V, local_length, myid);
  fails += Test_N_VWrmsNormMaskVectorArray(V, local_length, myid);
  fails += Test_N_VScaleAddMultiVectorArray(V, local_length, myid);
  fails += Test_N_VLinearCombinationVectorArray(V, local_length, myid);

  /* local reduction operations */
  if (myid == 0) printf("\nTesting local reduction operations:\n\n");

  fails += Test_N_VDotProdLocal(X, Y, local_length, myid);
  fails += Test_N_VMaxNormLocal(X, local_length, myid);
  fails += Test_N_VMinLocal(X, local_length, myid);
  fails += Test_N_VL1NormLocal(X, local_length, myid);
  fails += Test_N_VWSqrSumLocal(X, Y, local_length, myid);
  fails += Test_N_VWSqrSumMaskLocal(X, Y, Z, local_length, myid);
  fails += Test_N_VInvTestLocal(X, Z, local_length, myid);
  fails += Test_N_VConstrMaskLocal(X, Y, Z, local_length, myid);
  fails += Test_N_VMinQuotientLocal(X, Y, local_length, myid);

  /* local fused reduction operations */
  if (myid == 0) printf("\nTesting local fused reduction operations:\n\n");
  fails += Test_N_VDotProdMultiLocal(V, local_length, myid);
  fails += Test_N_VDotProdMultiAllReduce(V, local_length, myid);

  /* XBraid interface operations */
  if (myid == 0) printf("\nTesting XBraid interface operations:\n\n");

  fails += Test_N_VBufSize(X, local_length, myid);
  fails += Test_N_VBufPack(X, local_length, myid);
  fails += Test_N_VBufUnpack(X, local_length, myid);

  /* Free vectors */
  N_VDestroy(W);
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);
  N_VDestroy(U);
  N_VDestroy(V);
  N_VDestroy(Xlocal);

  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests, Proc %d \n\n", fails, myid);
  } else {
    if (myid == 0)
      printf("SUCCESS: NVector module passed all tests \n\n");
  }

  /* check if any other process failed */
  (void) MPI_Allreduce(&fails, &globfails, 1, MPI_INT, MPI_MAX, comm);

  MPI_Finalize();

  return(globfails);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int          failure = 0;
  sunindextype i;
  realtype     *x0;
  sunindextype x0len;

  x0len = N_VGetLocalLength_MPIPlusX(X);
  x0 = N_VGetArrayPointer(X);

  /* ensure that local_length = x0len + x1len */
  if (local_length != x0len)
    return(1);

  /* check vector data */
  for (i=0; i<x0len; i++)  failure += SUNRCompare(x0[i], ans);

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  /* should not be called in these tests */
  return SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  realtype *data;

  data = N_VGetArrayPointer(X);

  /* set i-th element of data array */
  data[i] = val;
}

void set_element_range(N_Vector X, sunindextype is, sunindextype ie, realtype val)
{
  sunindextype x0len, i;
  realtype *data;

  data = N_VGetArrayPointer(X);
  x0len = N_VGetLocalLength_MPIPlusX(X);

  /* set i-th element of data array */
  for (i=is; i<x0len; i++)  data[i] = val;
}

realtype get_element(N_Vector X, sunindextype i)
{
  realtype *data;

  data = N_VGetArrayPointer(X);

  /* get the i-th element of data array */
  return data[i];
}

double max_time(N_Vector X, double time)
{
  double maxt;
  MPI_Comm *comm;

  /* get max time across all MPI ranks */
  comm = (MPI_Comm *) N_VGetCommunicator(X);
  (void) MPI_Reduce(&time, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, *comm);
  return(maxt);
}

void sync_device(N_Vector x)
{
  /* not running on GPU, just return */
  return;
}
