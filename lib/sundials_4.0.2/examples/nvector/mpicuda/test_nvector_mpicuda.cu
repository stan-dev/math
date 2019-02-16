/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, and Cody J. Balos @ LLNL
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
 * This is the testing routine to check the NVECTOR CUDA module
 * implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/cuda/Vector.hpp>
#include <nvector/nvector_mpicuda.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>

/* CUDA vector specific tests */
int Test_N_VMake_Cuda(N_Vector X, sunindextype local_length, sunindextype global_length, MPI_Comm comm, int myid);
int Test_N_VMakeManaged_Cuda(N_Vector X, sunindextype local_length, sunindextype global_length, MPI_Comm comm, int myid);

/* CUDA vector can use unmanaged or managed memory */
enum mem_type { UNMANAGED, MANAGED };

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails = 0;         /* counter for test failures */
  int          globfails = 0;     /* counter for test failures */
  int          retval;            /* function return value     */
  sunindextype local_length;      /* local vector length       */
  sunindextype global_length;     /* global vector length      */
  N_Vector     U, V, X, Y, Z;     /* test vectors              */
  int          print_timing;      /* turn timing on/off        */
  MPI_Comm     comm;              /* MPI Communicator          */
  int          nprocs, myid;      /* Number of procs, proc id  */
  int          i;

  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* check inputs */
  if (argc < 3) {
    if (myid == 0)
      printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    MPI_Abort(comm, -1);
  }

  local_length = atol(argv[1]);
  if (local_length < 1) {
    if (myid == 0)
      printf("ERROR: local vector length must be a positive integer \n");
    MPI_Abort(comm, -1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, myid);

  /* global length */
  global_length = nprocs*local_length;

  for (i=UNMANAGED; i<=MANAGED; ++i) {
    if (myid == 0) {
      if (i==UNMANAGED) {
        printf("Testing CUDA N_Vector \n");
      } else {
        printf("\nTesting CUDA N_Vector with managed memory \n");
      }
      printf("Vector global length %ld \n", (long int) global_length);
      printf("MPI processes %d \n", nprocs);
    }

    /* Create new vectors */
    X = (i==UNMANAGED) ? N_VNew_Cuda(comm, local_length, global_length) : 
                         N_VNewManaged_Cuda(comm, local_length, global_length);
    if (X == NULL) {
      if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
      MPI_Abort(comm, 1);
    }

    /* Check vector ID */
    fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_CUDA, myid);

    /* Test clone functions */
    fails += Test_N_VCloneEmpty(X, myid);
    fails += Test_N_VClone(X, local_length, myid);
    fails += Test_N_VCloneEmptyVectorArray(5, X, myid);
    fails += Test_N_VCloneVectorArray(5, X, local_length, myid);

    /* Clone additional vectors for testing */
    Y = N_VClone(X);
    if (Y == NULL) {
      N_VDestroy(X);
      if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
      MPI_Abort(comm, 1);
    }

    Z = N_VClone(X);
    if (Z == NULL) {
      N_VDestroy(X);
      N_VDestroy(Y);
      if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
      MPI_Abort(comm, 1);
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
    fails += Test_N_VDotProd(X, Y, local_length, global_length, myid);
    fails += Test_N_VMaxNorm(X, local_length, myid);
    fails += Test_N_VWrmsNorm(X, Y, local_length, myid);
    fails += Test_N_VWrmsNormMask(X, Y, Z, local_length, global_length, myid);
    fails += Test_N_VMin(X, local_length, myid);
    fails += Test_N_VWL2Norm(X, Y, local_length, global_length, myid);
    fails += Test_N_VL1Norm(X, local_length, global_length, myid);
    fails += Test_N_VCompare(X, Z, local_length, myid);
    fails += Test_N_VInvTest(X, Z, local_length, myid);
    fails += Test_N_VConstrMask(X, Y, Z, local_length, myid);
    fails += Test_N_VMinQuotient(X, Y, local_length, myid);

    /* Fused and vector array operations tests (disabled) */
    if (myid == 0) printf("\nTesting fused and vector array operations (disabled):\n\n");

    /* create vector and disable all fused and vector array operations */
    U = (i==UNMANAGED) ? N_VNew_Cuda(comm, local_length, global_length) : 
                         N_VNewManaged_Cuda(comm, local_length, global_length);
    retval = N_VEnableFusedOps_Cuda(U, SUNFALSE);
    if (U == NULL || retval != 0) {
      N_VDestroy(X);
      N_VDestroy(Y);
      N_VDestroy(Z);
      if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
      MPI_Abort(comm, 1);
    }

    /* fused operations */
    fails += Test_N_VLinearCombination(U, local_length, myid);
    fails += Test_N_VScaleAddMulti(U, local_length, myid);
    fails += Test_N_VDotProdMulti(U, local_length, global_length, myid);

    /* vector array operations */
    fails += Test_N_VLinearSumVectorArray(U, local_length, myid);
    fails += Test_N_VScaleVectorArray(U, local_length, myid);
    fails += Test_N_VConstVectorArray(U, local_length, myid);
    fails += Test_N_VWrmsNormVectorArray(U, local_length, myid);
    fails += Test_N_VWrmsNormMaskVectorArray(U, local_length, global_length, myid);
    fails += Test_N_VScaleAddMultiVectorArray(U, local_length, myid);
    fails += Test_N_VLinearCombinationVectorArray(U, local_length, myid);

    /* Fused and vector array operations tests (enabled) */
    if (myid == 0) printf("\nTesting fused and vector array operations (enabled):\n\n");

    /* create vector and enable all fused and vector array operations */
    V = (i==UNMANAGED) ? N_VNew_Cuda(comm, local_length, global_length) : 
                         N_VNewManaged_Cuda(comm, local_length, global_length);
    retval = N_VEnableFusedOps_Cuda(V, SUNTRUE);
    if (V == NULL || retval != 0) {
      N_VDestroy(X);
      N_VDestroy(Y);
      N_VDestroy(Z);
      N_VDestroy(U);
      if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
      MPI_Abort(comm, 1);
    }

    /* fused operations */
    fails += Test_N_VLinearCombination(V, local_length, myid);
    fails += Test_N_VScaleAddMulti(V, local_length, myid);
    fails += Test_N_VDotProdMulti(V, local_length, global_length, myid);

    /* vector array operations */
    fails += Test_N_VLinearSumVectorArray(V, local_length, myid);
    fails += Test_N_VScaleVectorArray(V, local_length, myid);
    fails += Test_N_VConstVectorArray(V, local_length, myid);
    fails += Test_N_VWrmsNormVectorArray(V, local_length, myid);
    fails += Test_N_VWrmsNormMaskVectorArray(V, local_length, global_length, myid);
    fails += Test_N_VScaleAddMultiVectorArray(V, local_length, myid);
    fails += Test_N_VLinearCombinationVectorArray(V, local_length, myid);

    /* CUDA specific tests */
    if (i==UNMANAGED) {
      fails += Test_N_VMake_Cuda(X, local_length, global_length, comm, myid);
    } else {
      fails += Test_N_VMakeManaged_Cuda(X, local_length, global_length, comm, myid);
    }

    /* Free vectors */
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
    N_VDestroy(V);
  }

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
 * CUDA specific tests
 * --------------------------------------------------------------------*/

/* --------------------------------------------------------------------
 * Test for the CUDA N_Vector N_VMake_Cuda function. Requires N_VConst
 * to check data.
 */
int Test_N_VMake_Cuda(N_Vector X, sunindextype local_length, sunindextype global_length, MPI_Comm comm, int myid)
{
  int failure = 0;
  realtype *h_data, *d_data;
  N_Vector Y;

  N_VConst(NEG_HALF, X);
  N_VCopyFromDevice_Cuda(X);

  h_data = N_VGetHostArrayPointer_Cuda(X);
  d_data = N_VGetDeviceArrayPointer_Cuda(X);

  /* Case 1: h_data and d_data are not null */
  Y = N_VMake_Cuda(comm, local_length, global_length, h_data, d_data);
  if (Y == NULL) {
    printf(">>> FAILED test -- N_VMake_Cuda, Proc %d \n", myid);
    printf("    Vector is NULL \n \n");
    return(1);
  }

  if (N_VGetHostArrayPointer_Cuda(Y) == NULL) {
    printf(">>> FAILED test -- N_VMake_Cuda, Proc %d \n", myid);
    printf("    Vector host data == NULL \n \n");
    N_VDestroy(Y);
    return(1);
  }
  
  if (N_VGetDeviceArrayPointer_Cuda(Y) == NULL) {
    printf(">>> FAILED test -- N_VMake_Cuda, Proc %d \n", myid);
    printf("    Vector device data -= NULL \n \n");
    N_VDestroy(Y);
    return(1);
  }
  
  failure += check_ans(NEG_HALF, Y, local_length);
 
  if (failure) {
    printf(">>> FAILED test -- N_VMake_Cuda Case 1, Proc %d \n", myid);
    printf("    Failed N_VConst check \n \n");
    N_VDestroy(Y);
    return(1);
  }
  
  if (myid == 0) {
    printf("PASSED test -- N_VMake_Cuda Case 1 \n");
  }

  N_VDestroy(Y);

  /* Case 2: data is null */
  Y = N_VMake_Cuda(comm, local_length, global_length, NULL, NULL);
  if (Y != NULL) {
    printf(">>> FAILED test -- N_VMake_Cuda Case 2, Proc %d \n", myid);
    printf("    Vector is not NULL \n \n");
    return(1);
  }

  if (myid == 0) {
    printf("PASSED test -- N_VMake_Cuda Case 2 \n");
  }
  
  N_VDestroy(Y);

  return(failure);
}

/* --------------------------------------------------------------------
 * Test for the CUDA N_Vector N_VMakeManaged_Cuda function. Requires
 * N_VConst to check data. X must be using managed memory.
 */
int Test_N_VMakeManaged_Cuda(N_Vector X, sunindextype local_length, sunindextype global_length, MPI_Comm comm, int myid)
{
  int failure = 0;
  realtype *vdata;
  N_Vector Y;

  if(!N_VIsManagedMemory_Cuda(X)) {
    printf(">>> FAILED test -- N_VIsManagedMemory_Cuda, Proc %d \n", myid);
    return(1);
  }

  N_VConst(NEG_HALF, X);

  vdata = N_VGetHostArrayPointer_Cuda(X);
  
  /* Case 1: data is not null */
  Y = N_VMakeManaged_Cuda(comm, local_length, global_length, vdata);
  if (Y == NULL) {
    printf(">>> FAILED test -- N_VMakeManaged_Cuda, Proc %d \n", myid);
    printf("    Vector is NULL \n \n");
    return(1);
  }

  failure += check_ans(NEG_HALF, Y, local_length);
 
  /* Case 2: data is null */
  Y = N_VMakeManaged_Cuda(comm, local_length, global_length, NULL);
  if (Y != NULL) {
    printf(">>> FAILED test -- N_VMakeManaged_Cuda Case 2, Proc %d \n", myid);
    printf("    Vector is not NULL \n \n");
    return(1);
  }

  if (myid == 0) {
    printf("PASSED test -- N_VMakeManaged_Cuda Case 2 \n");
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

  N_VCopyFromDevice_Cuda(X);
  Xdata = N_VGetHostArrayPointer_Cuda(X);

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
  N_VCopyFromDevice_Cuda(X);
  (N_VGetHostArrayPointer_Cuda(X))[i] = val;
  N_VCopyToDevice_Cuda(X);
}

realtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  N_VCopyFromDevice_Cuda(X);
  return (N_VGetHostArrayPointer_Cuda(X))[i];
}

double max_time(N_Vector X, double time)
{
  MPI_Comm comm;
  double maxt;
 
  comm = N_VGetMPIComm_Cuda(X);
  
  /* get max time across all MPI ranks */
  (void) MPI_Reduce(&time, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  return(maxt);
}

void sync_device()
{
  /* sync with GPU */
  cudaDeviceSynchronize();
  return;
}
