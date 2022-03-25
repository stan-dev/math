/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_LAPACK_H
#define _SUNDIALS_LAPACK_H

#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ==================================================================
 * Blas and Lapack functions
 * ==================================================================
 */

#if defined(SUNDIALS_F77_FUNC)

#define dgbtrf_f77      SUNDIALS_F77_FUNC(dgbtrf, DGBTRF)
#define dgbtrs_f77      SUNDIALS_F77_FUNC(dgbtrs, DGBTRS)
#define dgetrf_f77      SUNDIALS_F77_FUNC(dgetrf, DGETRF)
#define dgetrs_f77      SUNDIALS_F77_FUNC(dgetrs, DGETRS)

#define sgbtrf_f77      SUNDIALS_F77_FUNC(sgbtrf, SGBTRF)
#define sgbtrs_f77      SUNDIALS_F77_FUNC(sgbtrs, SGBTRS)
#define sgetrf_f77      SUNDIALS_F77_FUNC(sgetrf, SGETRF)
#define sgetrs_f77      SUNDIALS_F77_FUNC(sgetrs, SGETRS)

#else

#define dgbtrf_f77      dgbtrf_
#define dgbtrs_f77      dgbtrs_
#define dgetrf_f77      dgetrf_
#define dgetrs_f77      dgetrs_

#define sgbtrf_f77      sgbtrf_
#define sgbtrs_f77      sgbtrs_
#define sgetrf_f77      sgetrf_
#define sgetrs_f77      sgetrs_

#endif

/* LAPACK */

extern void dgbtrf_f77(const sunindextype *m, const sunindextype *n,
                       const sunindextype *kl, const sunindextype *ku,
                       double *ab, sunindextype *ldab, sunindextype *ipiv,
                       sunindextype *info);

extern void dgbtrs_f77(const char *trans, const sunindextype *n,
                       const sunindextype *kl, const sunindextype *ku,
                       const sunindextype *nrhs, double *ab,
                       const sunindextype *ldab, sunindextype *ipiv,
                       double *b, const sunindextype *ldb, sunindextype *info);


extern void dgetrf_f77(const sunindextype *m, const sunindextype *n, double *a,
                       sunindextype *lda, sunindextype *ipiv,
                       sunindextype *info);

extern void dgetrs_f77(const char *trans, const sunindextype *n,
                       const sunindextype *nrhs, double *a,
                       const sunindextype *lda, sunindextype *ipiv, double *b,
                       const sunindextype *ldb, sunindextype *info);

extern void sgbtrf_f77(const sunindextype *m, const sunindextype *n,
                       const sunindextype *kl, const sunindextype *ku,
                       float *ab, sunindextype *ldab, sunindextype *ipiv,
                       sunindextype *info);

extern void sgbtrs_f77(const char *trans, const sunindextype *n,
                       const sunindextype *kl, const sunindextype *ku,
                       const sunindextype *nrhs, float *ab,
                       const sunindextype *ldab, sunindextype *ipiv,
                       float *b, const sunindextype *ldb, sunindextype *info);

extern void sgetrf_f77(const sunindextype *m, const sunindextype *n, float *a,
                       sunindextype *lda, sunindextype *ipiv,
                       sunindextype *info);

extern void sgetrs_f77(const char *trans, const sunindextype *n,
                       const sunindextype *nrhs, float *a,
                       const sunindextype *lda, sunindextype *ipiv,
                       float *b, const sunindextype *ldb, sunindextype *info);

#ifdef __cplusplus
}
#endif

#endif
