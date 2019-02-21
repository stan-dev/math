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
 * This header file contains definitions of MPI data types, which
 * are used by MPI parallel vector implementations.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_MPI_H
#define _SUNDIALS_MPI_H

#include <sundials/sundials_types.h>
#include <sundials/sundials_mpi_types.h>


#if SUNDIALS_MPI_ENABLED

#include <mpi.h>
#define SUNMPI_COMM_WORLD MPI_COMM_WORLD

typedef MPI_Comm SUNMPI_Comm;

#else

#define SUNMPI_COMM_WORLD 0

typedef int SUNMPI_Comm;

#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT int SUNMPI_Comm_size(SUNMPI_Comm comm, int *size);
SUNDIALS_EXPORT realtype SUNMPI_Allreduce_scalar(realtype d, int op, SUNMPI_Comm comm);
SUNDIALS_EXPORT void SUNMPI_Allreduce(realtype *d, int nvec, int op, SUNMPI_Comm comm);

#ifdef __cplusplus
}
#endif



#endif /* _SUNDIALS_MPI_H */
