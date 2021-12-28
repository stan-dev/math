/* -----------------------------------------------------------------
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * This is the header file for a generic BAND linear solver
 * package, based on the SUNDlsMat type defined in sundials_direct.h.
 *
 * There are two sets of band solver routines listed in
 * this file: one set uses type SUNDlsMat defined below and the
 * other set uses the type realtype ** for band matrix arguments.
 * Routines that work with the type SUNDlsMat begin with "Band".
 * Routines that work with realtype ** begin with "band".
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_BAND_H
#define _SUNDIALS_BAND_H

#include <sundials/sundials_direct.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_BandGBTRF
 * -----------------------------------------------------------------
 * Usage : ier = SUNDlsMat_BandGBTRF(A, p);
 *         if (ier != 0) ... A is singular
 * -----------------------------------------------------------------
 * SUNDlsMat_BandGBTRF performs the LU factorization of the N by N band
 * matrix A. This is done using standard Gaussian elimination
 * with partial pivoting.
 *
 * A successful LU factorization leaves the "matrix" A and the
 * pivot array p with the following information:
 *
 * (1) p[k] contains the row number of the pivot element chosen
 * at the beginning of elimination step k, k = 0, 1, ..., N-1.
 *
 * (2) If the unique LU factorization of A is given by PA = LU,
 *     where P is a permutation matrix, L is a lower triangular
 *     matrix with all 1's on the diagonal, and U is an upper
 *     triangular matrix, then the upper triangular part of A
 *     (including its diagonal) contains U and the strictly lower
 *     triangular part of A contains the multipliers, I-L.
 *
 * SUNDlsMat_BandGBTRF returns 0 if successful. Otherwise it encountered
 * a zero diagonal element during the factorization. In this case
 * it returns the column index (numbered from one) at which
 * it encountered the zero.
 *
 * Important Note: A must be allocated to accommodate the increase
 * in upper bandwidth that occurs during factorization. If
 * mathematically, A is a band matrix with upper bandwidth mu and
 * lower bandwidth ml, then the upper triangular factor U can
 * have upper bandwidth as big as smu = MIN(n-1,mu+ml). The lower
 * triangular factor L has lower bandwidth ml. Allocate A with
 * call A = BandAllocMat(N,mu,ml,smu), where mu, ml, and smu are
 * as defined above. The user does not have to zero the "extra"
 * storage allocated for the purpose of factorization. This will
 * handled by the SUNDlsMat_BandGBTRF routine.
 *
 * SUNDlsMat_BandGBTRF is only a wrapper around SUNDlsMat_bandGBTRF.
 * All work is done in SUNDlsMat_bandGBTRF, which works directly on the
 * data in the SUNDlsMat A (i.e. in the field A->cols).
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
sunindextype SUNDlsMat_BandGBTRF(SUNDlsMat A, sunindextype* p);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_BandGBTRF instead")
sunindextype BandGBTRF(DlsMat A, sunindextype *p);

SUNDIALS_EXPORT
sunindextype SUNDlsMat_bandGBTRF(realtype **a, sunindextype n,
                                 sunindextype mu, sunindextype ml,
                                 sunindextype smu, sunindextype *p);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_bandGBTRF instead")
sunindextype bandGBTRF(realtype **a, sunindextype n,
                       sunindextype mu, sunindextype ml,
                       sunindextype smu, sunindextype *p);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_BandGBTRS
 * -----------------------------------------------------------------
 * Usage: SUNDlsMat_BandGBTRS(A, p, b);
 * -----------------------------------------------------------------
 * SUNDlsMat_BandGBTRS solves the N-dimensional system A x = b using
 * the LU factorization in A and the pivot information in p computed
 * in SUNDlsMat_BandGBTRF. The solution x is returned in b. This
 * routine cannot fail if the corresponding call to
 * SUNDlsMat_BandGBTRF did not fail.
 *
 * SUNDlsMat_BandGBTRS is only a wrapper around SUNDlsMat_bandGBTRS
 * which does all the work directly on the data in the DlsMat A (i.e.
 * in A->cols).
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_BandGBTRS(SUNDlsMat A, sunindextype *p, realtype *b);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_BandGBTRS instead")
void BandGBTRS(DlsMat A, sunindextype *p, realtype *b);

SUNDIALS_EXPORT
void SUNDlsMat_bandGBTRS(realtype **a, sunindextype n, sunindextype smu,
                         sunindextype ml, sunindextype *p, realtype *b);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_bandGBTRS instead")
void bandGBTRS(realtype **a, sunindextype n, sunindextype smu,
               sunindextype ml, sunindextype *p, realtype *b);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_BandCopy
 * -----------------------------------------------------------------
 * Usage: SUNDlsMat_BandCopy(A, B, copymu, copyml);
 * -----------------------------------------------------------------
 * SUNDlsMat_BandCopy copies the submatrix with upper and lower
 * bandwidths copymu, copyml of the N by N band matrix A into the N by
 * N band matrix B.
 *
 * SUNDlsMat_BandCopy is a wrapper around SUNDlsMat_bandCopy which
 * accesses the data in the DlsMat A and DlsMat B (i.e. the fields
 * cols).
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_BandCopy(SUNDlsMat A, SUNDlsMat B, sunindextype copymu,
                        sunindextype copyml);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_BandCopy instead")
void BandCopy(DlsMat A, DlsMat B, sunindextype copymu,
              sunindextype copyml);

SUNDIALS_EXPORT
void SUNDlsMat_bandCopy(realtype **a, realtype **b, sunindextype n,
                        sunindextype a_smu, sunindextype b_smu,
                        sunindextype copymu, sunindextype copyml);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_bandCopy instead")
void bandCopy(realtype **a, realtype **b, sunindextype n,
              sunindextype a_smu, sunindextype b_smu,
              sunindextype copymu, sunindextype copyml);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_BandScale
 * -----------------------------------------------------------------
 * Usage: SUNDlsMat_BandScale(c, A);
 * -----------------------------------------------------------------
 * A(i,j) <- c*A(i,j), j-(A->mu) < = i < = j+(A->ml).
 *
 * SUNDlsMat_BandScale is a wrapper around SUNDlsMat_bandScale which
 * performs the actual scaling by accessing the data in the
 * SUNDlsMat A (i.e. the field A->cols).
 * -----------------------------------------------------------------
 */

void SUNDlsMat_BandScale(realtype c, SUNDlsMat A);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_BandScale instead")
void BandScale(realtype c, DlsMat A);

SUNDIALS_EXPORT
void SUNDlsMat_bandScale(realtype c, realtype **a, sunindextype n,
                         sunindextype mu, sunindextype ml, sunindextype smu);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_bandScale instead")
void bandScale(realtype c, realtype **a, sunindextype n,
               sunindextype mu, sunindextype ml, sunindextype smu);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_bandAddIdentity
 * -----------------------------------------------------------------
 * SUNDlsMat_bandAddIdentity adds the identity matrix to the n-by-n
 * matrix stored in the realtype** arrays.
 * -----------------------------------------------------------------
 */

void SUNDlsMat_bandAddIdentity(realtype **a, sunindextype n, sunindextype smu);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_bandAddIdentity instead")
void bandAddIdentity(realtype **a, sunindextype n, sunindextype smu);


/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_BandMatvec
 * -----------------------------------------------------------------
 * SUNDlsMat_BandMatvec computes the matrix-vector product y = A*x,
 * where A is an M-by-N band matrix, x is a vector of length N, and y
 * is a vector of length M.  No error checking is performed on the
 * length of the arrays x and y.  Only y is modified in this routine.
 *
 * SUNDlsMat_BandMatvec is a wrapper around SUNDlsMat_bandMatvec which
 * performs the actual product by accessing the data in the SUNDlsMat
 * A.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_BandMatvec(SUNDlsMat A, realtype *x, realtype *y);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_BandMatvec instead")
void BandMatvec(DlsMat A, realtype *x, realtype *y);

SUNDIALS_EXPORT
void SUNDlsMat_bandMatvec(realtype **a, realtype *x, realtype *y,
                          sunindextype n, sunindextype mu,
                          sunindextype ml, sunindextype smu);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_bandMatvec instead")
void bandMatvec(realtype **a, realtype *x, realtype *y,
                sunindextype n, sunindextype mu,
                sunindextype ml, sunindextype smu);

#ifdef __cplusplus
}
#endif

#endif
