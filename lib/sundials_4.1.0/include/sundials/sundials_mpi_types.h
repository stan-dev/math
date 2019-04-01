/* -----------------------------------------------------------------
 * Programmer(s): Scott Cohen, Alan Hindmarsh, Radu Serban,
 *                Aaron Collier, and Slaven Peles @ LLNL
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

#include <sundials/sundials_types.h>

/* define MPI data types */

#if defined(SUNDIALS_SINGLE_PRECISION)
  #define PVEC_REAL_MPI_TYPE MPI_FLOAT
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  #define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  #define PVEC_REAL_MPI_TYPE MPI_LONG_DOUBLE
#endif

#if defined(SUNDIALS_INT64_T)
  #define PVEC_INTEGER_MPI_TYPE MPI_INT64_T
#elif defined(SUNDIALS_INT32_T)
  #define PVEC_INTEGER_MPI_TYPE MPI_INT32_T
#endif
