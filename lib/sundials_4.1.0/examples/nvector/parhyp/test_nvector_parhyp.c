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
 * This is the testing routine to check the NVECTOR Parallel module
 * implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>

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

  HYPRE_Int       *partitioning;  /* Vector Partitioning               */
  HYPRE_ParVector Xhyp;           /* Instantiate hypre parallel vector */

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

  if (myid == 0) {
    printf("Testing the hypre ParVector N_Vector wrapper \n");
    printf("Vector global length %ld \n", (long int) global_length);
    printf("MPI processes %d \n\n", nprocs);
  }

  /* set partitioning */
  if(HYPRE_AssumedPartitionCheck()) {
    partitioning    = (HYPRE_Int*) malloc(2*sizeof(HYPRE_Int));
    partitioning[0] = myid*local_length;
    partitioning[1] = (myid+1)*local_length;
  } else {
    partitioning = (HYPRE_Int*) malloc((nprocs+1)*sizeof(HYPRE_Int));
    if (local_length < 1) {
      printf("Using global partition.\n");
      printf("I don't do this stuff. Now exiting...\n");
      MPI_Abort(comm, 1);
    }
  }

  /* Create hypre vector */
  HYPRE_ParVectorCreate(comm, global_length, partitioning, &Xhyp);
  HYPRE_ParVectorInitialize(Xhyp);

  /* Create hypre ParVector N_Vector wrapper and test */
  X = N_VMake_ParHyp(Xhyp);
  fails += Test_N_VMake(X, local_length, myid);
  if (fails != 0) {
    N_VDestroy(X);
    HYPRE_ParVectorDestroy(Xhyp);
    if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
    MPI_Abort(comm, 1);
  }

  /* Check vector ID */
  fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_PARHYP, myid);

  /* Test clone functions */
  fails += Test_N_VCloneEmpty(X, myid);
  fails += Test_N_VClone(X, local_length, myid);
  fails += Test_N_VCloneEmptyVectorArray(5, X, myid);
  fails += Test_N_VCloneVectorArray(5, X, local_length, myid);

  /* Clone additional vectors for testing */
  Y = N_VClone(X);
  if (Y == NULL) {
    N_VDestroy(X);
    HYPRE_ParVectorDestroy(Xhyp);
    if (myid == 0) printf("FAIL: Unable to create a new vector \n\n");
    MPI_Abort(comm, 1);
  }

  Z = N_VClone(X);
  if (Z == NULL) {
    N_VDestroy(X);
    N_VDestroy(Y);
    HYPRE_ParVectorDestroy(Xhyp);
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
  U = N_VMake_ParHyp(Xhyp);
  retval = N_VEnableFusedOps_ParHyp(U, SUNFALSE);
  if (U == NULL || retval != 0) {
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    HYPRE_ParVectorDestroy(Xhyp);
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
  V = N_VMake_ParHyp(Xhyp);
  retval = N_VEnableFusedOps_ParHyp(V, SUNTRUE);
  if (V == NULL || retval != 0) {
    N_VDestroy(X);
    N_VDestroy(Y);
    N_VDestroy(Z);
    N_VDestroy(U);
    HYPRE_ParVectorDestroy(Xhyp);
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

  /* Free vectors */
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);
  N_VDestroy(U);
  N_VDestroy(V);
  HYPRE_ParVectorDestroy(Xhyp);

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
  int             failure = 0;
  sunindextype    i;
  HYPRE_ParVector Xvec;
  realtype        *Xdata;

  Xvec  = N_VGetVector_ParHyp(X);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xvec));

  /* check vector data */
  for (i = 0; i < local_length; i++) {
    failure += FNEQ(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  /* check if wrapped hypre ParVector is non-null */
  return (N_VGetVector_ParHyp(X) == NULL) ? SUNFALSE : SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  HYPRE_ParVector Xvec;
  realtype        *Xdata;

  /* set i-th element of data array */
  Xvec  = N_VGetVector_ParHyp(X);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xvec));

  Xdata[i] = val;
}

realtype get_element(N_Vector X, sunindextype i)
{
  HYPRE_ParVector Xvec;
  realtype        *Xdata;

  /* get i-th element of data array */
  Xvec  = N_VGetVector_ParHyp(X);
  Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xvec));

  return Xdata[i];
}

double max_time(N_Vector X, double time)
{
  MPI_Comm comm;
  double maxt;

  /* get max time across all MPI ranks */
  comm = hypre_ParVectorComm(N_VGetVector_ParHyp(X));
  (void) MPI_Reduce(&time, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  return(maxt);
}

void sync_device()
{
  /* not running on GPU, just return */
  return;
}
