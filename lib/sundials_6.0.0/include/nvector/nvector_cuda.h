/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
 * This is the header file for the CUDA implementation of the
 * NVECTOR module.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_CUDA_H
#define _NVECTOR_CUDA_H

#include <cuda_runtime.h>
#include <stdio.h>

#include <sundials/sundials_cuda_policies.hpp>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>
#include <sunmemory/sunmemory_cuda.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * CUDA implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Cuda
{
  sunindextype       length;
  booleantype        own_exec;
  booleantype        own_helper;
  SUNMemory          host_data;
  SUNMemory          device_data;
  SUNCudaExecPolicy* stream_exec_policy;
  SUNCudaExecPolicy* reduce_exec_policy;
  SUNMemoryHelper    mem_helper;
  void*              priv; /* 'private' data */
};

typedef struct _N_VectorContent_Cuda *N_VectorContent_Cuda;

/*
 * -----------------------------------------------------------------
 * NVECTOR_CUDA implementation specific functions
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Cuda(SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNew_Cuda(sunindextype length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewManaged_Cuda(sunindextype length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewWithMemHelp_Cuda(sunindextype length,
                                                booleantype use_managed_mem,
                                                SUNMemoryHelper helper,
                                                SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Cuda(sunindextype length,
                                      realtype *h_vdata,
                                      realtype *d_vdata,
                                      SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMakeManaged_Cuda(sunindextype length,
                                             realtype *vdata,
                                             SUNContext sunctx);
SUNDIALS_EXPORT void N_VSetHostArrayPointer_Cuda(realtype* h_vdata, N_Vector v);
SUNDIALS_EXPORT void N_VSetDeviceArrayPointer_Cuda(realtype* d_vdata, N_Vector v);
SUNDIALS_EXPORT booleantype N_VIsManagedMemory_Cuda(N_Vector x);
SUNDIALS_EXPORT int N_VSetKernelExecPolicy_Cuda(N_Vector x,
                                                SUNCudaExecPolicy* stream_exec_policy,
                                                SUNCudaExecPolicy* reduce_exec_policy);
SUNDIALS_EXPORT void N_VCopyToDevice_Cuda(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Cuda(N_Vector v);

SUNDIALS_STATIC_INLINE
sunindextype N_VGetLength_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return content->length;
}

SUNDIALS_STATIC_INLINE
realtype *N_VGetHostArrayPointer_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return(content->host_data == NULL ? NULL : (realtype*)content->host_data->ptr);
}

SUNDIALS_STATIC_INLINE
realtype *N_VGetDeviceArrayPointer_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return(content->device_data == NULL ? NULL : (realtype*)content->device_data->ptr);
}

/*
 * -----------------------------------------------------------------
 * NVECTOR API functions
 * -----------------------------------------------------------------
 */

SUNDIALS_STATIC_INLINE
N_Vector_ID N_VGetVectorID_Cuda(N_Vector v)
{
  return SUNDIALS_NVEC_CUDA;
}

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Cuda(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Cuda(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Cuda(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Cuda(N_Vector v, sunindextype *lrw, sunindextype *liw);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Cuda(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Cuda(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Cuda(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Cuda(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Cuda(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Cuda(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Cuda(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Cuda(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Cuda(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Cuda(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Cuda(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Cuda(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Cuda(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Cuda(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Cuda(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Cuda(int nvec, realtype* c, N_Vector* X,
                                              N_Vector Z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Cuda(int nvec, realtype* c, N_Vector X,
                                          N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y,
                                         realtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Cuda(int nvec,
                                                 realtype a, N_Vector* X,
                                                 realtype b, N_Vector* Y,
                                                 N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Cuda(int nvec, realtype* c, N_Vector* X,
                                             N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Cuda(int nvec, realtype c, N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Cuda(int nvec, int nsum,
                                                     realtype* a, N_Vector* X,
                                                     N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Cuda(int nvec, int nsum,
                                                         realtype* c,
                                                         N_Vector** X,
                                                         N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X,
                                                N_Vector* W, realtype* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X,
                                                    N_Vector* W, N_Vector id,
                                                    realtype* nrm);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT realtype N_VWSqrSumLocal_Cuda(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWSqrSumMaskLocal_Cuda(N_Vector x, N_Vector w, N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT int N_VBufSize_Cuda(N_Vector x, sunindextype *size);
SUNDIALS_EXPORT int N_VBufPack_Cuda(N_Vector x, void *buf);
SUNDIALS_EXPORT int N_VBufUnpack_Cuda(N_Vector x, void *buf);

/* OPTIONAL operations for debugging */
SUNDIALS_EXPORT void N_VPrint_Cuda(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile_Cuda(N_Vector v, FILE *outfile);


/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Cuda(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf);

#ifdef __cplusplus
}
#endif

#endif
