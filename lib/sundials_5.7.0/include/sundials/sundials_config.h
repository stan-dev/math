/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos, Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS configuration header file.
 * -----------------------------------------------------------------*/

#include "sundials_export.h"

/* ------------------------------------------------------------------
 * Define SUNDIALS version numbers
 * -----------------------------------------------------------------*/


#define SUNDIALS_VERSION "5.7.0"
#define SUNDIALS_VERSION_MAJOR 5
#define SUNDIALS_VERSION_MINOR 7
#define SUNDIALS_VERSION_PATCH 0
#define SUNDIALS_VERSION_LABEL ""


/* ------------------------------------------------------------------
 * SUNDIALS build information
 * -----------------------------------------------------------------*/


/* Define precision of SUNDIALS data type 'realtype'
 * Depending on the precision level, one of the following
 * three macros will be defined:
 *     #define SUNDIALS_SINGLE_PRECISION 1
 *     #define SUNDIALS_DOUBLE_PRECISION 1
 *     #define SUNDIALS_EXTENDED_PRECISION 1
 */
#define SUNDIALS_DOUBLE_PRECISION 1

/* Define type of vector indices in SUNDIALS 'sunindextype'.
 * Depending on user choice of index type, one of the following
 * two macros will be defined:
 *     #define SUNDIALS_INT64_T 1
 *     #define SUNDIALS_INT32_T 1
 */
#define SUNDIALS_INT64_T 1

/* Define the type of vector indices in SUNDIALS 'sunindextype'.
 * The macro will be defined with a type of the appropriate size.
 */
#define SUNDIALS_INDEX_TYPE int64_t

/* Use generic math functions
 * If it was decided that generic math functions can be used, then
 *     #define SUNDIALS_USE_GENERIC_MATH
 */
#define SUNDIALS_USE_GENERIC_MATH

/* Use POSIX timers if available.
 *     #define SUNDIALS_HAVE_POSIX_TIMERS
 */
/* #undef SUNDIALS_HAVE_POSIX_TIMERS */

/* CVODE should use fused kernels if utilizing
  * the CUDA NVector.
  */
/* #undef SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS */


/* BUILD SUNDIALS with monitoring functionalities
  * the CUDA NVector.
  */
/* #undef SUNDIALS_BUILD_WITH_MONITORING */

/* ------------------------------------------------------------------
 * SUNDIALS TPL macros
 * -----------------------------------------------------------------*/




/* MAGMA backends */
/* #undef SUNDIALS_MAGMA_BACKENDS_HIP */
#define SUNDIALS_MAGMA_BACKENDS_CUDA

/* Set if SUNDIALS is built with MPI support, then
 *     #define SUNDIALS_MPI_ENABLED 1
 * otherwise
 *     #define SUNDIALS_MPI_ENABLED 0
 */
#define SUNDIALS_MPI_ENABLED 0

 /* SUPERLUMT threading type */
/* #undef SUNDIALS_SUPERLUMT_THREAD_TYPE */

 /* Trilinos with MPI is available, then
  *    #define SUNDIALS_TRILINOS_HAVE_MPI
  */
/* #undef SUNDIALS_TRILINOS_HAVE_MPI */

/* RAJA backends */
/* #undef SUNDIALS_RAJA_BACKENDS_HIP */
#define SUNDIALS_RAJA_BACKENDS_CUDA

/* ------------------------------------------------------------------
 * SUNDIALS modules enabled
 * -----------------------------------------------------------------*/


#define SUNDIALS_ARKODE 1
#define SUNDIALS_CVODE 1
#define SUNDIALS_CVODES 1
#define SUNDIALS_IDA 1
#define SUNDIALS_IDAS 1
#define SUNDIALS_KINSOL 1
#define SUNDIALS_NVECTOR_SERIAL 1
#define SUNDIALS_NVECTOR_MANYVECTOR 1
#define SUNDIALS_SUNMATRIX_BAND 1
#define SUNDIALS_SUNMATRIX_DENSE 1
#define SUNDIALS_SUNMATRIX_SPARSE 1
#define SUNDIALS_SUNLINSOL_BAND 1
#define SUNDIALS_SUNLINSOL_DENSE 1
#define SUNDIALS_SUNLINSOL_PCG 1
#define SUNDIALS_SUNLINSOL_SPBCGS 1
#define SUNDIALS_SUNLINSOL_SPFGMR 1
#define SUNDIALS_SUNLINSOL_SPGMR 1
#define SUNDIALS_SUNLINSOL_SPTFQMR 1
#define SUNDIALS_SUNNONLINSOL_NEWTON 1
#define SUNDIALS_SUNNONLINSOL_FIXEDPOINT 1



/* ------------------------------------------------------------------
 * SUNDIALS fortran configuration
 * -----------------------------------------------------------------*/


/* FCMIX: Define Fortran name-mangling macro for C identifiers.
 * Depending on the inferred scheme, one of the following six
 * macros will be defined:
 *     #define SUNDIALS_F77_FUNC(name,NAME) name
 *     #define SUNDIALS_F77_FUNC(name,NAME) name ## _
 *     #define SUNDIALS_F77_FUNC(name,NAME) name ## __
 *     #define SUNDIALS_F77_FUNC(name,NAME) NAME
 *     #define SUNDIALS_F77_FUNC(name,NAME) NAME ## _
 *     #define SUNDIALS_F77_FUNC(name,NAME) NAME ## __
 */


/* FCMIX: Define Fortran name-mangling macro for C identifiers
 *        which contain underscores.
 */


/* FNVECTOR: Allow user to specify different MPI communicator
 * If it was found that the MPI implementation supports MPI_Comm_f2c, then
 *      #define SUNDIALS_MPI_COMM_F2C 1
 * otherwise
 *      #define SUNDIALS_MPI_COMM_F2C 0
 */
#define SUNDIALS_MPI_COMM_F2C 1


/* ------------------------------------------------------------------
 * SUNDIALS inline macros.
 * -----------------------------------------------------------------*/


/* Mark SUNDIALS function as inline.
 */
#ifndef SUNDIALS_CXX_INLINE
#define SUNDIALS_CXX_INLINE inline
#endif

#ifndef SUNDIALS_C_INLINE
#define SUNDIALS_C_INLINE inline
#endif

#ifdef __cplusplus
#define SUNDIALS_INLINE SUNDIALS_CXX_INLINE
#else
#define SUNDIALS_INLINE SUNDIALS_C_INLINE
#endif

/* Mark SUNDIALS function as static inline.
 */
#define SUNDIALS_STATIC_INLINE static SUNDIALS_INLINE
