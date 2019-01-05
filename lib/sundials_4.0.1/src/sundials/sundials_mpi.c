/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
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


