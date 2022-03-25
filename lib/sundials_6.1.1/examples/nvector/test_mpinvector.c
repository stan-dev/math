/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL, Daniel R. Reynolds @ SMU
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
 * These are test functions for an NVECTOR module implementation
 * which have MPI symbols.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <sundials/sundials_nvector.h>

#include "test_nvector.h"

void Test_AbortMPI(void* comm, int code)
{
  Test_Finalize();
  MPI_Abort(*((MPI_Comm*)comm), code);
}

/* ----------------------------------------------------------------------
 * Test_N_VGetCommunicator Test (with MPI dependency).
 * --------------------------------------------------------------------*/
int Test_N_VGetCommunicatorMPI(N_Vector W, void *comm, int myid)
{
  void* wcomm;
  MPI_Comm* Wcomm;
  MPI_Comm* Comm;
  int same;

  /* ask W for its communicator */
  wcomm = NULL;
  wcomm = N_VGetCommunicator(W);

  /* return with success if both are NULL */
  if ((wcomm == NULL) && (comm == NULL))  {
    printf("PASSED test -- N_VGetCommunicator\n");
    return(0);
  }

  /* return with failure if either is NULL */
  if (wcomm == NULL) {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (incorrectly reports NULL comm)\n", myid);
    return(1);
  }
  if (comm == NULL) {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (incorrectly reports non-NULL comm)\n", myid);
    return(1);
  }

  /* call MPI_Comm_compare to check that communicators match or are congruent */
  Wcomm = (MPI_Comm *) wcomm;
  Comm  = (MPI_Comm *) comm;
  if (MPI_Comm_compare(*Comm, *Wcomm, &same) != MPI_SUCCESS) {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (error in MPI_Comm_compare)\n", myid);
    return(1);
  }
  if ((same != MPI_IDENT) && (same != MPI_CONGRUENT)) {
    printf(">>> FAILED test -- N_VGetCommunicator, Proc %d (mismatched comms)\n", myid);
    return(1);
  }
  if (myid == 0)
    printf("PASSED test -- N_VGetCommunicator\n");
  return(0);
}
