/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the main header file for the PETSc vector wrapper
 * for NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type realtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     build configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type booleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum_Petsc(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_PETSC_H
#define _NVECTOR_PETSC_H

#include <mpi.h>
#include <petscvec.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_mpi_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PETSc implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Petsc {
  sunindextype local_length;   /* copy of local vector length  */
  sunindextype global_length;  /* copy of global vector length */
  booleantype own_data;        /* ownership of data            */
  Vec pvec;                    /* the PETSc Vec object         */
  MPI_Comm comm;               /* copy of MPI communicator     */
};

typedef struct _N_VectorContent_Petsc *N_VectorContent_Petsc;

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_petsc
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Petsc(MPI_Comm comm,
                                           sunindextype local_length,
                                           sunindextype global_length,
                                           SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VMake_Petsc(Vec v, SUNContext sunctx);

SUNDIALS_EXPORT realtype *N_VGetArrayPointer_Petsc(N_Vector v);

SUNDIALS_EXPORT Vec N_VGetVector_Petsc(N_Vector v);

SUNDIALS_EXPORT void N_VSetVector_Petsc(N_Vector v, Vec p);

SUNDIALS_EXPORT void N_VPrint_Petsc(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_Petsc(N_Vector v, const char fname[]);

/* nvector API functions */
SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Petsc(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Petsc(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Petsc(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Petsc(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Petsc(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT void N_VSetArrayPointer_Petsc(realtype *v_data, N_Vector v);
SUNDIALS_EXPORT void *N_VGetCommunicator_Petsc(N_Vector v);
SUNDIALS_EXPORT sunindextype N_VGetLength_Petsc(N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Petsc(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Petsc(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Petsc(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Petsc(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Petsc(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Petsc(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Petsc(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Petsc(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Petsc(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Petsc(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Petsc(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Petsc(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Petsc(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Petsc(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Petsc(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Petsc(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Petsc(int nvec, realtype* c,
                                               N_Vector* X, N_Vector z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Petsc(int nvec, realtype* a, N_Vector x,
                                           N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Petsc(int nvec, N_Vector x, N_Vector* Y,
                                          realtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Petsc(int nvec,
                                                  realtype a, N_Vector* X,
                                                  realtype b, N_Vector* Y,
                                                  N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Petsc(int nvec, realtype* c,
                                              N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Petsc(int nvecs, realtype c,
                                              N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Petsc(int nvecs, N_Vector* X,
                                                 N_Vector* W, realtype* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Petsc(int nvec, N_Vector* X,
                                                     N_Vector* W, N_Vector id,
                                                     realtype* nrm);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Petsc(int nvec, int nsum,
                                                      realtype* a,
                                                      N_Vector* X,
                                                      N_Vector** Y,
                                                      N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Petsc(int nvec, int nsum,
                                                          realtype* c,
                                                          N_Vector** X,
                                                          N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT realtype N_VDotProdLocal_Petsc(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNormLocal_Petsc(N_Vector x);
SUNDIALS_EXPORT realtype N_VMinLocal_Petsc(N_Vector x);
SUNDIALS_EXPORT realtype N_VL1NormLocal_Petsc(N_Vector x);
SUNDIALS_EXPORT realtype N_VWSqrSumLocal_Petsc(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWSqrSumMaskLocal_Petsc(N_Vector x, N_Vector w,
                                                   N_Vector id);
SUNDIALS_EXPORT booleantype N_VInvTestLocal_Petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMaskLocal_Petsc(N_Vector c, N_Vector x,
                                                     N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotientLocal_Petsc(N_Vector num,
                                                   N_Vector denom);

/* OPTIONAL single buffer reduction operations */
SUNDIALS_EXPORT int N_VDotProdMultiLocal_Petsc(int nvec, N_Vector x,
                                               N_Vector* Y, realtype* dotprods);
SUNDIALS_EXPORT int N_VDotProdMultiAllReduce_Petsc(int nvec, N_Vector x,
                                                   realtype* sum);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT int N_VBufSize_Petsc(N_Vector x, sunindextype *size);
SUNDIALS_EXPORT int N_VBufPack_Petsc(N_Vector x, void *buf);
SUNDIALS_EXPORT int N_VBufUnpack_Petsc(N_Vector x, void *buf);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Petsc(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Petsc(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Petsc(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableDotProdMulti_Petsc(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Petsc(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Petsc(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Petsc(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Petsc(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Petsc(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Petsc(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Petsc(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableDotProdMultiLocal_Petsc(N_Vector v, booleantype tf);

/*
 * -----------------------------------------------------------------
 * Deprecated functions
 * -----------------------------------------------------------------
 */

/* use N_VCloneVectorArray */
SUNDIALS_DEPRECATED_EXPORT N_Vector *N_VCloneVectorArray_Petsc(int count, N_Vector w);

/* use N_VCloneVectorArrayEmpty */
SUNDIALS_DEPRECATED_EXPORT N_Vector *N_VCloneVectorArrayEmpty_Petsc(int count, N_Vector w);

/* use N_VDestroyVectorArray */
SUNDIALS_DEPRECATED_EXPORT void N_VDestroyVectorArray_Petsc(N_Vector *vs, int count);

#ifdef __cplusplus
}
#endif

#endif /* _NVECTOR_PETSC_H */
