/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * Test for Trilinos vector kernels used with N_Vector.
 * -----------------------------------------------------------------
 * For unit testing set OMP_PROC_BIND=false and OMP_NUM_THREADS=2
 * -----------------------------------------------------------------
 */

#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <nvector/trilinos/SundialsTpetraVectorInterface.hpp>
#include <nvector/nvector_trilinos.h>
#include "test_nvector.h"

using namespace sundials::trilinos::nvector_tpetra;

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main (int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  /* Define SUNDIALS compatible map and vector types */
  typedef TpetraVectorInterface::vector_type vector_type;
  typedef vector_type::map_type map_type;

  Test_Init(NULL);

  /* Start an MPI session */
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  /* Create Tpetra communicator */
  auto comm = Tpetra::getDefaultComm();

  const int myRank = comm->getRank();

  /* check inputs */
  if (argc < 3) {
    if (myRank == 0)
      printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    Test_Abort(1);
  }

  const sunindextype local_length = (sunindextype) atol(argv[1]);
  if (local_length < 1) {
    if (myRank == 0)
      printf("ERROR: local vector length must be a positive integer \n");
    Test_Abort(1);
  }

  int print_timing = atoi(argv[2]);
  SetTiming(print_timing, myRank);

  /* Make partitioning easy */
  const sunindextype global_length = comm->getSize() * local_length;

  if (myRank == 0) {
    printf("Testing the Trilinos (Tpetra) N_Vector wrapper \n");
    printf("Vector global length %ld \n\n", (long int) global_length);
  }

  /* Choose zero-based (C-style) indexing. */
  const sunindextype index_base = 0;

  /* Construct am MPI Map */
  RCP<const map_type> testMap =
    rcp(new map_type (global_length, index_base, comm,
                      Tpetra::GloballyDistributed));

  /* Construct a Tpetra vector and return refernce counting pointer to it. */
  RCP<vector_type> px = rcp(new vector_type(testMap));

  int fails = 0;       /* counter for test failures */
  int globfails = 0;   /* counter for test failures */

  /* NVector Test */

  /* Create Trilinos (Tpetra) N_Vector wrapper and test */
  N_Vector X = N_VMake_Trilinos(px, sunctx);
  fails += Test_N_VMake(X, local_length, myRank);
  if (fails != 0) {
    N_VDestroy(X);
    px = Teuchos::null;
    if (myRank == 0) printf("FAIL: Unable to create a new vector \n\n");
    Test_Abort(1);
  }

  /* Check vector ID */
  fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_TRILINOS, myRank);

  /* Check vector length */
  fails += Test_N_VGetLength(X, myRank);

  /* Check vector communicator */
#ifdef SUNDIALS_TRILINOS_HAVE_MPI
  auto mpicomm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  fails += Test_N_VGetCommunicatorMPI(X, (MPI_Comm *) mpicomm->getRawMpiComm().get(), myRank);
#else
  fails += Test_N_VGetCommunicator(X, NULL, myRank);
#endif

  /* Test clone functions */
  fails += Test_N_VCloneEmpty(X, myRank);
  fails += Test_N_VClone(X, local_length, myRank);
  fails += Test_N_VCloneEmptyVectorArray(5, X, myRank);
  fails += Test_N_VCloneVectorArray(5, X, local_length, myRank);

  /* Clone additional vectors for testing */
  N_Vector Y = N_VClone(X);
  if (Y == NULL) {
    N_VDestroy(X);
    px = Teuchos::null;
    if (myRank == 0) printf("FAIL: Unable to create a new vector \n\n");
    Test_Abort(1);
  }

  N_Vector Z = N_VClone(X);
  if (Z == NULL) {
    N_VDestroy(X);
    N_VDestroy(Y);
    px = Teuchos::null;
    if (myRank == 0) printf("FAIL: Unable to create a new vector \n\n");
    Test_Abort(1);
  }

  /* Standard vector operation tests */
  if (myRank == 0) printf("\nTesting standard vector operations:\n\n");

  fails += Test_N_VConst(X, local_length, myRank);
  fails += Test_N_VLinearSum(X, Y, Z, local_length, myRank);
  fails += Test_N_VProd(X, Y, Z, local_length, myRank);
  fails += Test_N_VDiv(X, Y, Z, local_length, myRank);
  fails += Test_N_VScale(X, Z, local_length, myRank);
  fails += Test_N_VAbs(X, Z, local_length, myRank);
  fails += Test_N_VInv(X, Z, local_length, myRank);
  fails += Test_N_VAddConst(X, Z, local_length, myRank);
  fails += Test_N_VDotProd(X, Y, local_length, myRank);
  fails += Test_N_VMaxNorm(X, local_length, myRank);
  fails += Test_N_VWrmsNorm(X, Y, local_length, myRank);
  fails += Test_N_VWrmsNormMask(X, Y, Z, local_length, myRank);
  fails += Test_N_VMin(X, local_length, myRank);
  fails += Test_N_VWL2Norm(X, Y, local_length, myRank);
  fails += Test_N_VL1Norm(X, local_length, myRank);
  fails += Test_N_VCompare(X, Z, local_length, myRank);
  fails += Test_N_VInvTest(X, Z, local_length, myRank);
  fails += Test_N_VConstrMask(X, Y, Z, local_length, myRank);
  fails += Test_N_VMinQuotient(X, Y, local_length, myRank);

  /* local reduction operations */
  if (myRank == 0) printf("\nTesting local reduction operations:\n\n");

  fails += Test_N_VDotProdLocal(X, Y, local_length, myRank);
  fails += Test_N_VMaxNormLocal(X, local_length, myRank);
  fails += Test_N_VMinLocal(X, local_length, myRank);
  fails += Test_N_VL1NormLocal(X, local_length, myRank);
  fails += Test_N_VWSqrSumLocal(X, Y, local_length, myRank);
  fails += Test_N_VWSqrSumMaskLocal(X, Y, Z, local_length, myRank);
  fails += Test_N_VInvTestLocal(X, Z, local_length, myRank);
  fails += Test_N_VConstrMaskLocal(X, Y, Z, local_length, myRank);
  fails += Test_N_VMinQuotientLocal(X, Y, local_length, myRank);

  /* Free vectors */
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);

  /* Print result */
  if (fails)
    printf("FAIL: NVector module failed %i tests, Proc %d \n \n", fails, myRank);

  /* Check if any other process failed */
  Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_SUM, fails,
                                Teuchos::outArg (globfails));

  /* Print global result */
  if (myRank == 0) {
    if (globfails)
      printf("FAIL: NVector module failed total of %i tests across all processes \n \n", globfails);
    else
      printf("SUCCESS: NVector module passed all tests on all processes \n \n");
  }

  px = Teuchos::null;

  return(globfails);
}

/* ----------------------------------------------------------------------
 * Check vector
 * --------------------------------------------------------------------*/

/*
 * Checks if all elements of vector X are _ans_
 */
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  Teuchos::RCP<TpetraVectorInterface::vector_type> xv =
    N_VGetVector_Trilinos(X);

  /* Sync the host with the device if needed */
  xv->sync<Kokkos::HostSpace>();
  const auto x_2d = xv->getLocalView<Kokkos::HostSpace>();
  const auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);

  int failure = 0;
  sunindextype i;

  /* check Tpetra vector */
  for (i = 0; i < local_length; ++i){
    failure += SUNRCompare(x_1d(i), ans);
  }

  return (failure > ZERO) ? 1 : 0;
}

/*
 * Checks if there is a Tpetra vector
 */
booleantype has_data(N_Vector X)
{
  if (X->content == NULL)
    return SUNFALSE;
  if (N_VGetVector_Trilinos(X).getRawPtr() == nullptr)
    return SUNFALSE;
  else
    return SUNTRUE;
}

/*
 * Sets ith element of vector X to val
 */
void set_element(N_Vector X, sunindextype i, realtype val)
{
  set_element_range(X, i, i, val);
}

/*
 * Sets elements [is, ie] of vector X to val
 */
void set_element_range(N_Vector X, sunindextype is, sunindextype ie,
                       realtype val)
{
  typedef TpetraVectorInterface::vector_type vector_type;
  typedef vector_type::node_type::memory_space memory_space;

  Teuchos::RCP<vector_type> xv = N_VGetVector_Trilinos(X);

  /* Sync the host with the device if needed */
  xv->sync<Kokkos::HostSpace>();
  const auto x_2d = xv->getLocalView<Kokkos::HostSpace>();
  const auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);

  xv->modify<Kokkos::HostSpace>();

  sunindextype i;
  for(i = is; i <= ie; i++) x_1d(i) = val;

  /* Sync the device with the host */
  xv->sync<memory_space>();
}

/*
 * Returns ith element of vector X
 */
realtype get_element(N_Vector X, sunindextype i)
{
  Teuchos::RCP<TpetraVectorInterface::vector_type> xv =
    N_VGetVector_Trilinos(X);

  /* Sync the host with the device if needed */
  xv->sync<Kokkos::HostSpace>();
  const auto x_2d = xv->getLocalView<Kokkos::HostSpace>();
  const auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);

  return x_1d(i);
}

double max_time(N_Vector X, double time)
{
  double maxtime = 0.0;
  Teuchos::RCP<TpetraVectorInterface::vector_type> xv =
    N_VGetVector_Trilinos(X);
  auto comm = xv->getMap()->getComm();
  Teuchos::reduceAll<int, double> (*comm, Teuchos::REDUCE_SUM, time,
                                   Teuchos::outArg(maxtime));
  return maxtime;
}

void sync_device(N_Vector x)
{
  /* Kokkos should take care of this */
}
