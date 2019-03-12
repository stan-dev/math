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
 * This is implementation of SUNDIALS MPI wrapper functions.
 * -----------------------------------------------------------------*/

#include <sundials/sundials_mpi.h>

int SUNMPI_Comm_size(SUNMPI_Comm comm, int *size)
{
#if SUNDIALS_MPI_ENABLED
  return MPI_Comm_size(comm, size);
#else
  *size = 1;
  return 0;
#endif
}

realtype SUNMPI_Allreduce_scalar(realtype d, int op, SUNMPI_Comm comm)
{
  /*
   * This function does a global reduction.  The operation is
   *   sum if op = 1,
   *   max if op = 2,
   *   min if op = 3.
   * The operation is over all processors in the communicator
   */

#if SUNDIALS_MPI_ENABLED

  realtype out;

  switch (op) {
   case 1: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);

#else

  /* If MPI is not enabled don't do reduction */
  return d;

#endif /* ifdef SUNDIALS_MPI_ENABLED */
}


void SUNMPI_Allreduce(realtype *d, int nvec, int op, SUNMPI_Comm comm)
{
  /*
   * This function does a global reduction.  The operation is
   *   sum if op = 1,
   *   max if op = 2,
   *   min if op = 3.
   * The operation is over all processors in the communicator
   */

#if SUNDIALS_MPI_ENABLED

  switch (op) {
   case 1: MPI_Allreduce(MPI_IN_PLACE, d, nvec, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(MPI_IN_PLACE, d, nvec, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(MPI_IN_PLACE, d, nvec, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

#else

  /* If MPI is not enabled don't do reduction */

#endif /* ifdef SUNDIALS_MPI_ENABLED */
}


