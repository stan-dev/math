/* ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
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
 * This is the main header file for the MPI-enabled implementation
 * of the NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type realtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type booleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum_Parallel(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_PARALLEL_H
#define _NVECTOR_PARALLEL_H

#include <stdio.h>
#include <mpi.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_mpi_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Parallel implementation of N_Vector               
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Parallel {
  sunindextype local_length;   /* local vector length         */
  sunindextype global_length;  /* global vector length        */
  booleantype own_data;        /* ownership of data           */
  realtype *data;              /* local data array            */
  MPI_Comm comm;               /* pointer to MPI communicator */
};

typedef struct _N_VectorContent_Parallel *N_VectorContent_Parallel;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_P, NV_DATA_P, NV_OWN_DATA_P,
 *        NV_LOCLENGTH_P, NV_GLOBLENGTH_P,NV_COMM_P, and NV_Ith_P
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_P(v)    ( (N_VectorContent_Parallel)(v->content) )

#define NV_LOCLENGTH_P(v)  ( NV_CONTENT_P(v)->local_length )

#define NV_GLOBLENGTH_P(v) ( NV_CONTENT_P(v)->global_length )

#define NV_OWN_DATA_P(v)   ( NV_CONTENT_P(v)->own_data )

#define NV_DATA_P(v)       ( NV_CONTENT_P(v)->data )

#define NV_COMM_P(v)       ( NV_CONTENT_P(v)->comm )

#define NV_Ith_P(v,i)      ( NV_DATA_P(v)[i] )

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_parallel
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_Parallel(MPI_Comm comm, 
                                         sunindextype local_length,
                                         sunindextype global_length);

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Parallel(MPI_Comm comm, 
                                              sunindextype local_length,
                                              sunindextype global_length);

SUNDIALS_EXPORT N_Vector N_VMake_Parallel(MPI_Comm comm, 
                                          sunindextype local_length,
                                          sunindextype global_length,
                                          realtype *v_data);

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_Parallel(int count, N_Vector w);

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_Parallel(int count, N_Vector w);

SUNDIALS_EXPORT void N_VDestroyVectorArray_Parallel(N_Vector *vs, int count);

SUNDIALS_EXPORT sunindextype N_VGetLength_Parallel(N_Vector v);

SUNDIALS_EXPORT sunindextype N_VGetLocalLength_Parallel(N_Vector v);

SUNDIALS_EXPORT void N_VPrint_Parallel(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_Parallel(N_Vector v, FILE *outfile);

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Parallel(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Parallel(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Parallel(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Parallel(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Parallel(N_Vector v, sunindextype *lrw,
                                       sunindextype *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_Parallel(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_Parallel(realtype *v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Parallel(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Parallel(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Parallel(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Parallel(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Parallel(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Parallel(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Parallel(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Parallel(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Parallel(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Parallel(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Parallel(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Parallel(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Parallel(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Parallel(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Parallel(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Parallel(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Parallel(int nvec, realtype* c, N_Vector* V,
                                                  N_Vector z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Parallel(int nvec, realtype* a, N_Vector x,
                                              N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Parallel(int nvec, N_Vector x,
                                             N_Vector *Y, realtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Parallel(int nvec, 
                                                     realtype a, N_Vector* X,
                                                     realtype b, N_Vector* Y,
                                                     N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Parallel(int nvec, realtype* c,
                                                 N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Parallel(int nvecs, realtype c,
                                                 N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Parallel(int nvecs, N_Vector* X,
                                                    N_Vector* W, realtype* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Parallel(int nvec, N_Vector* X,
                                                        N_Vector* W, N_Vector id,
                                                        realtype* nrm);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Parallel(int nvec, int nsum,
                                                         realtype* a,
                                                         N_Vector* X,
                                                         N_Vector** Y,
                                                         N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Parallel(int nvec, int nsum,
                                                             realtype* c,
                                                             N_Vector** X,
                                                             N_Vector* Z);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Parallel(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Parallel(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Parallel(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableDotProdMulti_Parallel(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Parallel(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Parallel(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Parallel(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Parallel(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Parallel(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Parallel(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Parallel(N_Vector v, booleantype tf);

#ifdef __cplusplus
}
#endif

#endif
