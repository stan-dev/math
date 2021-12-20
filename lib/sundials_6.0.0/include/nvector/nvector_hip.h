/* -----------------------------------------------------------------
 * Programmer(s): Daniel McGreer and Cody J. Balos @ LLNL
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
 * This is the header file for the hip implementation of the
 * NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be found
 *     in the header file sundials_nvector.h.
 *
 *   - The definitions of the types 'realtype' and 'sunindextype' can
 *     be found in the header file sundials_types.h, and it may be
 *     changed (at the configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype'.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *       N_VLinearSum_Hip(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_HIP_H
#define _NVECTOR_HIP_H

#include <hip/hip_runtime.h>
#include <stdio.h>

#include <sundials/sundials_hip_policies.hpp>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>
#include <sunmemory/sunmemory_hip.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * hip implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Hip
{
  sunindextype       length;
  booleantype        own_exec;
  booleantype        own_helper;
  SUNMemory          host_data;
  SUNMemory          device_data;
  SUNHipExecPolicy*  stream_exec_policy;
  SUNHipExecPolicy*  reduce_exec_policy;
  SUNMemoryHelper    mem_helper;
  void*              priv; /* 'private' data */
};

typedef struct _N_VectorContent_Hip *N_VectorContent_Hip;

/*
 * -----------------------------------------------------------------
 * NVECTOR_HIP implementation specific functions
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_Hip(sunindextype length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewManaged_Hip(sunindextype length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewWithMemHelp_Hip(sunindextype length,
                                               booleantype use_managed_mem,
                                               SUNMemoryHelper helper,
                                               SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewEmpty_Hip(SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Hip(sunindextype length,
                                     realtype *h_vdata,
                                     realtype *d_vdata,
                                     SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMakeManaged_Hip(sunindextype length,
                                            realtype *vdata,
                                            SUNContext sunctx);
SUNDIALS_EXPORT void N_VSetHostArrayPointer_Hip(realtype* h_vdata, N_Vector v);
SUNDIALS_EXPORT booleantype N_VIsManagedMemory_Hip(N_Vector x);
SUNDIALS_EXPORT int N_VSetKernelExecPolicy_Hip(N_Vector x,
                                               SUNHipExecPolicy* stream_exec_policy,
                                               SUNHipExecPolicy* reduce_exec_policy);
SUNDIALS_EXPORT void N_VCopyToDevice_Hip(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Hip(N_Vector v);
SUNDIALS_EXPORT void N_VPrint_Hip(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile_Hip(N_Vector v, FILE *outfile);

SUNDIALS_STATIC_INLINE
sunindextype N_VGetLength_Hip(N_Vector x)
{
  N_VectorContent_Hip content = (N_VectorContent_Hip)x->content;
  return content->length;
}

SUNDIALS_STATIC_INLINE
realtype *N_VGetHostArrayPointer_Hip(N_Vector x)
{
  N_VectorContent_Hip content = (N_VectorContent_Hip)x->content;
  return(content->host_data == NULL ? NULL : (realtype*)content->host_data->ptr);
}

SUNDIALS_STATIC_INLINE
realtype *N_VGetDeviceArrayPointer_Hip(N_Vector x)
{
  N_VectorContent_Hip content = (N_VectorContent_Hip)x->content;
  return(content->device_data == NULL ? NULL : (realtype*)content->device_data->ptr);
}

/*
 * -----------------------------------------------------------------
 * NVECTOR API functions
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Hip(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Hip(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Hip(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Hip(N_Vector v, sunindextype *lrw, sunindextype *liw);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Hip(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Hip(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Hip(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Hip(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Hip(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Hip(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Hip(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Hip(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Hip(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Hip(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Hip(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Hip(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Hip(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Hip(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Hip(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Hip(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Hip(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Hip(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Hip(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Hip(int nvec, realtype* c, N_Vector* X,
                                             N_Vector Z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Hip(int nvec, realtype* c, N_Vector X,
                                         N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Hip(int nvec, N_Vector x, N_Vector* Y,
                                        realtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Hip(int nvec,
                                                realtype a, N_Vector* X,
                                                realtype b, N_Vector* Y,
                                                N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Hip(int nvec, realtype* c, N_Vector* X,
                                            N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Hip(int nvec, realtype c, N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Hip(int nvec, int nsum,
                                                    realtype* a, N_Vector* X,
                                                    N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Hip(int nvec, int nsum,
                                                        realtype* c,
                                                        N_Vector** X,
                                                        N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Hip(int nvec, N_Vector* X,
                                               N_Vector* W, realtype* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Hip(int nvec, N_Vector* X,
                                                   N_Vector* W, N_Vector id,
                                                   realtype* nrm);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT realtype N_VWSqrSumLocal_Hip(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWSqrSumMaskLocal_Hip(N_Vector x, N_Vector w, N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT int N_VBufSize_Hip(N_Vector x, sunindextype *size);
SUNDIALS_EXPORT int N_VBufPack_Hip(N_Vector x, void *buf);
SUNDIALS_EXPORT int N_VBufUnpack_Hip(N_Vector x, void *buf);

/* OPTIONAL operations for debugging */
SUNDIALS_EXPORT void N_VPrint_Hip(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile_Hip(N_Vector v, FILE *outfile);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Hip(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Hip(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Hip(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableDotProdMulti_Hip(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Hip(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Hip(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Hip(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Hip(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Hip(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Hip(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Hip(N_Vector v, booleantype tf);

#ifdef __cplusplus
}
#endif

#endif
