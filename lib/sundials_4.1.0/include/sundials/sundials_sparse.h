/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Carol Woodward, Slaven Peles @ LLNL,
 *             Daniel R. Reynolds @ SMU.
 * -----------------------------------------------------------------
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This header file contains definitions and declarations for use by
 * sparse linear solvers for Ax = b. 
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_SPARSE_H
#define _SUNDIALS_SPARSE_H

#include <stdio.h>

#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ==================================================================
 * Type definitions
 * ==================================================================
 */

#define CSC_MAT 0
#define CSR_MAT 1

/*
 * Type : SlsMat
 */

typedef struct _SlsMat {
  int M;
  int N;
  int NNZ;
  int NP;
  realtype *data;
  int sparsetype;
  int *indexvals;
  int *indexptrs;
  /* CSC indices */
  int **rowvals;
  int **colptrs;
  /* CSR indices */
  int **colvals;
  int **rowptrs;
} *SlsMat;

/*
 * ==================================================================
 * Exported function prototypes (functions working on SlsMat)
 * ==================================================================
 */

SUNDIALS_EXPORT SlsMat SparseNewMat(int M, int N, int NNZ, int sparsetype);

SUNDIALS_EXPORT SlsMat SparseFromDenseMat(const DlsMat A, int sparsetype);

SUNDIALS_EXPORT int SparseDestroyMat(SlsMat A);

SUNDIALS_EXPORT int SparseSetMatToZero(SlsMat A);

SUNDIALS_EXPORT int SparseCopyMat(const SlsMat A, SlsMat B);

SUNDIALS_EXPORT int SparseScaleMat(realtype b, SlsMat A);

SUNDIALS_EXPORT int SparseAddIdentityMat(SlsMat A);

SUNDIALS_EXPORT int SparseAddMat(SlsMat A, const SlsMat B);

SUNDIALS_EXPORT int SparseReallocMat(SlsMat A);

SUNDIALS_EXPORT int SparseMatvec(const SlsMat A, const realtype *x, realtype *y);

SUNDIALS_EXPORT void SparsePrintMat(const SlsMat A, FILE* outfile);


#ifdef __cplusplus
}
#endif

#endif
