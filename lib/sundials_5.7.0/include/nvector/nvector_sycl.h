/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the header file for the SYCL implementation of the
 * NVECTOR module.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_SYCL_H
#define _NVECTOR_SYCL_H

#include <CL/sycl.hpp>
#include <stdio.h>

#include <sundials/sundials_sycl_policies.hpp>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>
#include <sunmemory/sunmemory_sycl.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* -----------------------------------------------------------------
 * SYCL implementation of N_Vector
 * ----------------------------------------------------------------- */


struct _N_VectorContent_Sycl
{
  sunindextype       length;
  booleantype        own_exec;
  booleantype        own_helper;
  SUNMemory          host_data;
  SUNMemory          device_data;
  SUNSyclExecPolicy* stream_exec_policy;
  SUNSyclExecPolicy* reduce_exec_policy;
  SUNMemoryHelper    mem_helper;
  sycl::queue*       queue;
  void*              priv; /* 'private' data */
};

typedef struct _N_VectorContent_Sycl *N_VectorContent_Sycl;


/* -----------------------------------------------------------------
 * NVECTOR_SYCL implementation specific functions
 * ----------------------------------------------------------------- */


SUNDIALS_EXPORT N_Vector N_VNewEmpty_Sycl();
SUNDIALS_EXPORT N_Vector N_VNew_Sycl(sunindextype length,
                                     sycl::queue *Q);
SUNDIALS_EXPORT N_Vector N_VNewManaged_Sycl(sunindextype length,
                                            sycl::queue *Q);
SUNDIALS_EXPORT N_Vector N_VNewWithMemHelp_Sycl(sunindextype length,
                                                booleantype use_managed_mem,
                                                SUNMemoryHelper helper,
                                                sycl::queue *Q);
SUNDIALS_EXPORT N_Vector N_VMake_Sycl(sunindextype length,
                                      realtype *h_vdata,
                                      realtype *d_vdata,
                                      sycl::queue *Q);
SUNDIALS_EXPORT N_Vector N_VMakeManaged_Sycl(sunindextype length,
                                             realtype *vdata,
                                             sycl::queue *Q);

SUNDIALS_EXPORT void N_VSetHostArrayPointer_Sycl(realtype* h_vdata, N_Vector v);
SUNDIALS_EXPORT void N_VSetDeviceArrayPointer_Sycl(realtype* d_vdata,
                                                   N_Vector v);
SUNDIALS_EXPORT booleantype N_VIsManagedMemory_Sycl(N_Vector x);
SUNDIALS_EXPORT int N_VSetKernelExecPolicy_Sycl(N_Vector x,
                                                SUNSyclExecPolicy* stream_exec_policy,
                                                SUNSyclExecPolicy* reduce_exec_policy);
SUNDIALS_EXPORT void N_VCopyToDevice_Sycl(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Sycl(N_Vector v);

SUNDIALS_STATIC_INLINE
sunindextype N_VGetLength_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return content->length;
}

SUNDIALS_STATIC_INLINE
realtype *N_VGetHostArrayPointer_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return(content->host_data == NULL ? NULL : (realtype*)content->host_data->ptr);
}

SUNDIALS_STATIC_INLINE
realtype *N_VGetDeviceArrayPointer_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return(content->device_data == NULL ? NULL : (realtype*)content->device_data->ptr);
}


/* -----------------------------------------------------------------
 * NVECTOR API functions
 * ----------------------------------------------------------------- */


SUNDIALS_STATIC_INLINE
N_Vector_ID N_VGetVectorID_Sycl(N_Vector v)
{
  return SUNDIALS_NVEC_SYCL;
}

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Sycl(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Sycl(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Sycl(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Sycl(N_Vector v, sunindextype *lrw,
                                   sunindextype *liw);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Sycl(realtype a, N_Vector x,
                                       realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Sycl(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Sycl(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Sycl(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Sycl(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Sycl(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Sycl(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Sycl(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Sycl(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Sycl(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Sycl(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Sycl(N_Vector x, N_Vector w,
                                              N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Sycl(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Sycl(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Sycl(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Sycl(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Sycl(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Sycl(N_Vector c, N_Vector x,
                                               N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Sycl(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Sycl(int nvec, realtype* c, N_Vector* X,
                                              N_Vector Z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Sycl(int nvec, realtype* c, N_Vector X,
                                          N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Sycl(int nvec, N_Vector x, N_Vector* Y,
                                         realtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Sycl(int nvec,
                                                 realtype a, N_Vector* X,
                                                 realtype b, N_Vector* Y,
                                                 N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Sycl(int nvec, realtype* c, N_Vector* X,
                                             N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Sycl(int nvec, realtype c, N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Sycl(int nvec, int nsum,
                                                     realtype* a, N_Vector* X,
                                                     N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Sycl(int nvec, int nsum,
                                                         realtype* c,
                                                         N_Vector** X,
                                                         N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Sycl(int nvec, N_Vector* X,
                                                N_Vector* W, realtype* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Sycl(int nvec, N_Vector* X,
                                                    N_Vector* W, N_Vector id,
                                                    realtype* nrm);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT realtype N_VWSqrSumLocal_Sycl(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWSqrSumMaskLocal_Sycl(N_Vector x, N_Vector w,
                                                  N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT int N_VBufSize_Sycl(N_Vector x, sunindextype *size);
SUNDIALS_EXPORT int N_VBufPack_Sycl(N_Vector x, void *buf);
SUNDIALS_EXPORT int N_VBufUnpack_Sycl(N_Vector x, void *buf);

/* OPTIONAL operations for debugging */
SUNDIALS_EXPORT void N_VPrint_Sycl(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile_Sycl(N_Vector v, FILE *outfile);


/* -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * ----------------------------------------------------------------- */


SUNDIALS_EXPORT int N_VEnableFusedOps_Sycl(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Sycl(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Sycl(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableDotProdMulti_Sycl(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Sycl(N_Vector v,
                                                       booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Sycl(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Sycl(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Sycl(N_Vector v,
                                                      booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Sycl(N_Vector v,
                                                          booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Sycl(N_Vector v,
                                                           booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Sycl(N_Vector v,
                                                               booleantype tf);

#ifdef __cplusplus
}
#endif

#endif
