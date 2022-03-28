/* -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
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
 * This header file contains definitions and declarations for use by
 * generic direct linear solvers for Ax = b. It defines types for
 * dense and banded matrices and corresponding accessor macros.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_DIRECT_H
#define _SUNDIALS_DIRECT_H

#include <stdio.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *                C O N S T A N T S
 * =================================================================
 */

/*
 *  SUNDIALS_DENSE: dense matrix
 *  SUNDIALS_BAND:  banded matrix
 */

#define SUNDIALS_DENSE 1
#define SUNDIALS_BAND  2

/*
 * ==================================================================
 * Type definitions
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type : SUNDlsMat
 * -----------------------------------------------------------------
 * The type SUNDlsMat is defined to be a pointer to a structure
 * with various sizes, a data field, and an array of pointers to
 * the columns which defines a dense or band matrix for use in
 * direct linear solvers. The M and N fields indicates the number
 * of rows and columns, respectively. The data field is a one
 * dimensional array used for component storage. The cols field
 * stores the pointers in data for the beginning of each column.
 * -----------------------------------------------------------------
 * For DENSE matrices, the relevant fields in SUNDlsMat are:
 *    type  = SUNDIALS_DENSE
 *    M     - number of rows
 *    N     - number of columns
 *    ldim  - leading dimension (ldim >= M)
 *    data  - pointer to a contiguous block of realtype variables
 *    ldata - length of the data array =ldim*N
 *    cols  - array of pointers. cols[j] points to the first element
 *            of the j-th column of the matrix in the array data.
 *
 * The elements of a dense matrix are stored columnwise (i.e. columns
 * are stored one on top of the other in memory).
 * If A is of type SUNDlsMat, then the (i,j)th element of A (with
 * 0 <= i < M and 0 <= j < N) is given by (A->data)[j*n+i].
 *
 * The SUNDLS_DENSE_COL and SUNDLS_DENSE_ELEM macros below allow a
 * user to access efficiently individual matrix elements without
 * writing out explicit data structure references and without knowing
 * too much about the underlying element storage. The only storage
 * assumption needed is that elements are stored columnwise and that a
 * pointer to the jth column of elements can be obtained via the
 * SUNDLS_DENSE_COL macro.
 * -----------------------------------------------------------------
 * For BAND matrices, the relevant fields in SUNDlsMat are:
 *    type  = SUNDIALS_BAND
 *    M     - number of rows
 *    N     - number of columns
 *    mu    - upper bandwidth, 0 <= mu <= min(M,N)
 *    ml    - lower bandwidth, 0 <= ml <= min(M,N)
 *    s_mu  - storage upper bandwidth, mu <= s_mu <= N-1.
 *            The dgbtrf routine writes the LU factors into the storage
 *            for A. The upper triangular factor U, however, may have
 *            an upper bandwidth as big as MIN(N-1,mu+ml) because of
 *            partial pivoting. The s_mu field holds the upper
 *            bandwidth allocated for A.
 *    ldim  - leading dimension (ldim >= s_mu)
 *    data  - pointer to a contiguous block of realtype variables
 *    ldata - length of the data array =ldim*(s_mu+ml+1)
 *    cols  - array of pointers. cols[j] points to the first element
 *            of the j-th column of the matrix in the array data.
 *
 * The SUNDLS_BAND_COL, SUNDLS_BAND_COL_ELEM, and SUNDLS_BAND_ELEM
 * macros below allow a user to access individual matrix elements
 * without writing out explicit data structure references and without
 * knowing too much about the underlying element storage. The only
 * storage assumption needed is that elements are stored columnwise
 * and that a pointer into the jth column of elements can be obtained
 * via the SUNDLS_BAND_COL macro. The SUNDLS_BAND_COL_ELEM macro
 * selects an element from a column which has already been isolated
 * via SUNDLS_BAND_COL. The macro SUNDLS_BAND_COL_ELEM allows the user
 * to avoid the translation from the matrix location (i,j) to the
 * index in the array returned by SUNDLS_BAND_COL at which the (i,j)th
 * element is stored.
 * -----------------------------------------------------------------
 */

typedef struct _DlsMat {
  int type;
  sunindextype M;
  sunindextype N;
  sunindextype ldim;
  sunindextype mu;
  sunindextype ml;
  sunindextype s_mu;
  realtype *data;
  sunindextype ldata;
  realtype **cols;
} *SUNDlsMat; /* DEPRECATED DlsMat: use SUNDlsMat instead */

typedef SUNDlsMat DlsMat;

/*
 * ==================================================================
 * Data accessor macros
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * SUNDLS_DENSE_COL and SUNDLS_DENSE_ELEM
 * -----------------------------------------------------------------
 *
 * SUNDLS_DENSE_COL(A,j) references the jth column of the M-by-N dense
 * matrix A, 0 <= j < N. The type of the expression SUNDLS_DENSE_COL(A,j)
 * is (realtype *). After the assignment col_j = SUNDLS_DENSE_COL(A,j),
 * col_j may be treated as an array indexed from 0 to M-1. The (i,j)-th
 * element of A is thus referenced by * col_j[i].
 *
 * SUNDLS_DENSE_ELEM(A,i,j) references the (i,j)th element of the dense
 * M-by-N matrix A, 0 <= i < M ; 0 <= j < N.
 *
 * -----------------------------------------------------------------
 */

#define SUNDLS_DENSE_COL(A,j) ((A->cols)[j])
#define SUNDLS_DENSE_ELEM(A,i,j) ((A->cols)[j][i])

/* DEPRECATED DENSE_COL: use SUNDLS_DENSE_COL instead */
#define DENSE_COL(A,j) SUNDLS_DENSE_COL(A,j)
/* DEPRECATED DENSE_ELEM: use SUNDLS_DENSE_ELEM instead */
#define DENSE_ELEM(A,i,j) SUNDLS_DENSE_ELEM(A,i,j)

/*
 * -----------------------------------------------------------------
 * SUNDLS_BAND_COL, SUNDLS_BAND_COL_ELEM, and SUNDLS_BAND_ELEM
 * -----------------------------------------------------------------
 *
 * SUNDLS_BAND_COL(A,j) references the diagonal element of the jth
 * column of the N by N band matrix A, 0 <= j <= N-1. The type of the
 * expression SUNDLS_BAND_COL(A,j) is realtype *. The pointer returned
 * by the call SUNDLS_BAND_COL(A,j) can be treated as an array which
 * is indexed from -(A->mu) to (A->ml).
 *
 * SUNDLS_BAND_COL_ELEM references the (i,j)th entry of the band
 * matrix A when used in conjunction with SUNDLS_BAND_COL. The index
 * (i,j) should satisfy j-(A->mu) <= i <= j+(A->ml).
 *
 * SUNDLS_BAND_ELEM(A,i,j) references the (i,j)th element of the
 * M-by-N band matrix A, where 0 <= i,j <= N-1. The location (i,j)
 * should further satisfy j-(A->mu) <= i <= j+(A->ml).
 *
 * -----------------------------------------------------------------
 */

#define SUNDLS_BAND_COL(A,j) (((A->cols)[j])+(A->s_mu))
#define SUNDLS_BAND_COL_ELEM(col_j,i,j) (col_j[(i)-(j)])
#define SUNDLS_BAND_ELEM(A,i,j) ((A->cols)[j][(i)-(j)+(A->s_mu)])

/* DEPRECATED BAND_COL: use SUNDLS_BAND_COL */
#define BAND_COL(A,j) SUNDLS_BAND_COL(A,j)
/* DEPRECATED BAND_COL_ELEM: use SUNDLS_BAND_COL_ELEM */
#define BAND_COL_ELEM(col_j,i,j) SUNDLS_BAND_COL_ELEM(col_j,i,j)
/* DEPRECATED BAND_ELEM: use SUNDLS_BAND_ELEM */
#define BAND_ELEM(A,i,j) SUNDLS_BAND_ELEM(A,i,j)

/*
 * ==================================================================
 * Exported function prototypes (functions working on SUNDlsMat)
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_NewDenseMat
 * -----------------------------------------------------------------
 * SUNDlsMat_NewDenseMat allocates memory for an M-by-N dense matrix
 * and returns the storage allocated (type SUNDlsMat).
 * SUNDlsMat_NewDenseMat returns NULL if the request for matrix
 * storage cannot be satisfied. See the above documentation for the
 * type SUNDlsMat for matrix storage details.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
SUNDlsMat SUNDlsMat_NewDenseMat(sunindextype M, sunindextype N);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_NewDenseMat instead")
DlsMat NewDenseMat(sunindextype M, sunindextype N);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_NewBandMat
 * -----------------------------------------------------------------
 * SUNDlsMat_NewBandMat allocates memory for an M-by-N band matrix
 * with upper bandwidth mu, lower bandwidth ml, and storage upper
 * bandwidth smu. Pass smu as follows depending on whether A will be
 * LU factored:
 *
 * (1) Pass smu = mu if A will not be factored.
 *
 * (2) Pass smu = MIN(N-1,mu+ml) if A will be factored.
 *
 * SUNDlsMat_NewBandMat returns the storage allocated (type SUNDlsMat)
 * or NULL if the request for matrix storage cannot be satisfied. See
 * the documentation for the type SUNDlsMat for matrix storage
 * details.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
SUNDlsMat SUNDlsMat_NewBandMat(sunindextype N, sunindextype mu,
                               sunindextype ml, sunindextype smu);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_NewBandMat instead")
DlsMat NewBandMat(sunindextype N, sunindextype mu,
                  sunindextype ml, sunindextype smu);

/*
 * -----------------------------------------------------------------
 * Functions: SUNDlsMat_DestroyMat
 * -----------------------------------------------------------------
 * SUNDlsMat_DestroyMat frees the memory allocated by
 * SUNDlsMat_NewDenseMat or SUNDlsMat_NewBandMat
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_DestroyMat(DlsMat A);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_DestroyMat instead")
void DestroyMat(DlsMat A);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_NewIntArray
 * -----------------------------------------------------------------
 * SUNDlsMat_NewIntArray allocates memory an array of N int's and
 * returns the pointer to the memory it allocates. If the request for
 * memory storage cannot be satisfied, it returns NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int* SUNDlsMat_NewIntArray(int N);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_NewIntArray instead")
int* NewIntArray(int N);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_NewIndexArray
 * -----------------------------------------------------------------
 * SUNDlsMat_NewIndexArray allocates memory an array of N
 * sunindextype's and returns the pointer to the memory it
 * allocates. If the request for memory storage cannot be satisfied,
 * it returns NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
sunindextype* SUNDlsMat_NewIndexArray(sunindextype N);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_NewIndexArray instead")
sunindextype* NewIndexArray(sunindextype N);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_NewRealArray
 * -----------------------------------------------------------------
 * SUNDlsMat_NewRealArray allocates memory an array of N realtype and
 * returns the pointer to the memory it allocates. If the request for
 * memory storage cannot be satisfied, it returns NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
realtype* SUNDlsMat_NewRealArray(sunindextype N);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_NewRealArray instead")
realtype* NewRealArray(sunindextype N);

/*
 * -----------------------------------------------------------------
 * Function: SUNDlsMat_DestroyArray
 * -----------------------------------------------------------------
 * SUNDlsMat_DestroyArray frees memory allocated by
 * SUNDlsMat_NewIntArray, SUNDlsMat_NewIndexArray, or
 * SUNDlsMat_NewRealArray.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_DestroyArray(void *p);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_DestroyArray instead")
void DestroyArray(void *p);

/*
 * -----------------------------------------------------------------
 * Function : SUNDlsMat_AddIdentity
 * -----------------------------------------------------------------
 * SUNDlsMat_AddIdentity adds 1.0 to the main diagonal (A_ii,
 * i=0,1,...,N-1) of the M-by-N matrix A (M>= N) and stores the result
 * back in A.  SUNDlsMat_AddIdentity is typically used with square
 * matrices.  SUNDlsMat_AddIdentity does not check for M >= N and
 * therefore a segmentation fault will occur if M < N!
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_AddIdentity(SUNDlsMat A);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_AddIdentity instead")
void AddIdentity(DlsMat A);

/*
 * -----------------------------------------------------------------
 * Function : SUNDlsMat_SetToZero
 * -----------------------------------------------------------------
 * SUNDlsMat_SetToZero sets all the elements of the M-by-N matrix A
 * to 0.0.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_SetToZero(SUNDlsMat A);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_SetToZero instead")
void SetToZero(DlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: SUNDlsMat_PrintMat
 * -----------------------------------------------------------------
 * This function prints the M-by-N (dense or band) matrix A to
 * outfile as it would normally appear on paper.
 * It is intended as debugging tools with small values of M and N.
 * The elements are printed using the %g/%lg/%Lg option.
 * A blank line is printed before and after the matrix.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_PrintMat(SUNDlsMat A, FILE *outfile);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_PrintMat")
void PrintMat(DlsMat A, FILE *outfile);

/*
 * ==================================================================
 * Exported function prototypes (functions working on realtype**)
 * ==================================================================
 */

SUNDIALS_EXPORT
realtype** SUNDlsMat_newDenseMat(sunindextype m, sunindextype n);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_newDenseMat instead")
realtype** newDenseMat(sunindextype m, sunindextype n);

SUNDIALS_EXPORT
realtype** SUNDlsMat_newBandMat(sunindextype n, sunindextype smu,
                                sunindextype ml);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_newBandMat instead")
realtype** newBandMat(sunindextype n, sunindextype smu,
                      sunindextype ml);

SUNDIALS_EXPORT
void SUNDlsMat_destroyMat(realtype** a);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_destroyMat instead")
void destroyMat(realtype** a);

SUNDIALS_EXPORT
int* SUNDlsMat_newIntArray(int n);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_newIntArray instead")
int* newIntArray(int n);

SUNDIALS_EXPORT
sunindextype* SUNDlsMat_newIndexArray(sunindextype n);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_newIndexArray instead")
sunindextype* newIndexArray(sunindextype n);

SUNDIALS_EXPORT
realtype* SUNDlsMat_newRealArray(sunindextype m);

SUNDIALS_EXPORT
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_newRealArray instead")
  realtype* newRealArray(sunindextype m);

SUNDIALS_EXPORT
void SUNDlsMat_destroyArray(void* v);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNDlsMat_destroyArray instead")
void destroyArray(void* v);


#ifdef __cplusplus
}
#endif

#endif
