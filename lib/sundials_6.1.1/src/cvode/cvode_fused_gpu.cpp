/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
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
 * This file implements fused CUDA/HIP kernels for CVODE.
 * -----------------------------------------------------------------*/

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include <nvector/nvector_cuda.h>
#include "sundials_cuda_kernels.cuh"
using SUNExecPolicy = SUNCudaExecPolicy;
using NVectorContent = N_VectorContent_Cuda;
constexpr auto gpuDeviceSynchronize = cudaDeviceSynchronize;
constexpr auto gpuGetLastError = cudaGetLastError;
constexpr auto gpuAssert = SUNDIALS_CUDA_Assert;
#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
#define SUNDIALS_DEBUG_GPU_LASTERROR
#endif

#elif USE_HIP
#include <hip/hip_runtime.h>
#include <nvector/nvector_hip.h>
#include "sundials_hip_kernels.hip.hpp"
using SUNExecPolicy = SUNHipExecPolicy;
using NVectorContent = N_VectorContent_Hip;
constexpr auto gpuDeviceSynchronize = hipDeviceSynchronize;
constexpr auto gpuGetLastError = hipGetLastError;
constexpr auto gpuAssert = SUNDIALS_HIP_Assert;
#ifdef SUNDIALS_DEBUG_HIP_LASTERROR
#define SUNDIALS_DEBUG_GPU_LASTERROR
#endif

#else
#error Incompatible GPU option for fused kernels
#endif

/*
 * -----------------------------------------------------------------
 * Compute the ewt vector when the tol type is CV_SS.
 * -----------------------------------------------------------------
 */


__global__
void cvEwtSetSS_kernel(const sunindextype length,
                       const realtype reltol,
                       const realtype Sabstol,
                       const realtype* ycur,
                       realtype* tempv,
                       realtype* weight)
{
  const realtype one = 1.0;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VAbs(ycur, cv_mem->cv_tempv);
    // N_VScale(cv_mem->cv_reltol, cv_mem->cv_tempv, cv_mem->cv_tempv);
    // N_VAddConst(cv_mem->cv_tempv, cv_mem->cv_Sabstol, cv_mem->cv_tempv);
    // N_VInv(cv_mem->cv_tempv, weight);
    realtype tmp = abs(ycur[i]);
    tempv[i] = reltol*tmp + Sabstol;
    weight[i] = one/tempv[i];
  }
}

extern "C"
int cvEwtSetSS_fused(const booleantype atolMin0,
                     const realtype reltol,
                     const realtype Sabstol,
                     const N_Vector ycur,
                     N_Vector tempv,
                     N_Vector weight)
{
  const SUNExecPolicy* exec_policy = ((NVectorContent)weight->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(weight);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvEwtSetSS_kernel<<<grid, block, 0, *(exec_policy->stream())>>>
  (
    N,
    reltol,
    Sabstol,
    N_VGetDeviceArrayPointer(ycur),
    N_VGetDeviceArrayPointer(tempv),
    N_VGetDeviceArrayPointer(weight)
  );

#ifdef SUNDIALS_DEBUG_GPU_LASTERROR
  gpuDeviceSynchronize();
  if (!gpuAssert(gpuGetLastError(), __FILE__, __LINE__)) return -1;
#endif

  return 0;
}


/*
 * -----------------------------------------------------------------
 * Compute the ewt vector when the tol type is CV_SV.
 * -----------------------------------------------------------------
 */


__global__
void cvEwtSetSV_kernel(const sunindextype length,
                       const realtype reltol,
                       const realtype* Vabstol,
                       const realtype* ycur,
                       realtype* tempv,
                       realtype* weight)
{
  const realtype one = 1.0;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VAbs(ycur, cv_mem->cv_tempv);
    // N_VLinearSum(cv_mem->cv_reltol, cv_mem->cv_tempv, ONE,
    //             cv_mem->cv_Vabstol, cv_mem->cv_tempv);
    // N_VInv(cv_mem->cv_tempv, weight);
    realtype tmp = abs(ycur[i]);
    tempv[i] = reltol*tmp + Vabstol[i];
    weight[i] = one/tempv[i];
  }
}

extern "C"
int cvEwtSetSV_fused(const booleantype atolMin0,
                     const realtype reltol,
                     const N_Vector Vabstol,
                     const N_Vector ycur,
                     N_Vector tempv,
                     N_Vector weight)
{
  const SUNExecPolicy* exec_policy = ((NVectorContent)weight->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(weight);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvEwtSetSV_kernel<<<grid, block, 0, *(exec_policy->stream())>>>
  (
    N,
    reltol,
    N_VGetDeviceArrayPointer(Vabstol),
    N_VGetDeviceArrayPointer(ycur),
    N_VGetDeviceArrayPointer(tempv),
    N_VGetDeviceArrayPointer(weight)
  );

#ifdef SUNDIALS_DEBUG_GPU_LASTERROR
  gpuDeviceSynchronize();
  if (!gpuAssert(gpuGetLastError(), __FILE__, __LINE__)) return -1;
#endif

  return 0;
}


/*
 * -----------------------------------------------------------------
 * Determine if the constraints of the problem are satisfied by
 * the proposed step.
 * -----------------------------------------------------------------
 */


__global__
void cvCheckConstraints_kernel(const sunindextype length,
                               const realtype* c,
                               const realtype* ewt,
                               const realtype* y,
                               const realtype* mm,
                               realtype* tempv)
{
  static const realtype zero = 0.0;
  static const realtype pt1 = 0.1;
  static const realtype one = 1.0;
  static const realtype onept5 = 1.5;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VCompare(ONEPT5, cv_mem->cv_constraints, tmp); /* a[i]=1 when |c[i]|=2  */
    // N_VProd(tmp, cv_mem->cv_constraints, tmp);       /* a * c                 */
    // N_VDiv(tmp, cv_mem->cv_ewt, tmp);                /* a * c * wt            */
    // N_VLinearSum(ONE, cv_mem->cv_y, -PT1, tmp, tmp); /* y - 0.1 * a * c * wt  */
    // N_VProd(tmp, mm, tmp);                           /* v = mm*(y-0.1*a*c*wt) */
    realtype tmp = (abs(c[i]) >= onept5) ? one : zero;
    tmp = tmp*c[i];
    tmp = tmp/ewt[i];
    tmp = y[i] - pt1*tmp;
    tempv[i] = tmp*mm[i];
  }
}

extern "C"
int cvCheckConstraints_fused(const N_Vector c,
                             const N_Vector ewt,
                             const N_Vector y,
                             const N_Vector mm,
                             N_Vector tempv)
{
  const SUNExecPolicy* exec_policy = ((NVectorContent)c->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(c);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvCheckConstraints_kernel<<<grid, block, 0, *(exec_policy->stream())>>>
  (
    N,
    N_VGetDeviceArrayPointer(c),
    N_VGetDeviceArrayPointer(ewt),
    N_VGetDeviceArrayPointer(y),
    N_VGetDeviceArrayPointer(mm),
    N_VGetDeviceArrayPointer(tempv)
  );

#ifdef SUNDIALS_DEBUG_GPU_LASTERROR
  gpuDeviceSynchronize();
  if (!gpuAssert(gpuGetLastError(), __FILE__, __LINE__)) return -1;
#endif

  return 0;
}


/*
 * -----------------------------------------------------------------
 * Compute the nonlinear residual.
 * -----------------------------------------------------------------
 */


__global__
void cvNlsResid_kernel(const sunindextype length,
                       const realtype rl1,
                       const realtype ngamma,
                       const realtype* zn1,
                       const realtype* ycor,
                       const realtype* ftemp,
                       realtype* res)
{
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VLinearSum(cv_mem->cv_rl1, cv_mem->cv_zn[1], ONE, ycor, res);
    // N_VLinearSum(-cv_mem->cv_gamma, cv_mem->cv_ftemp, ONE, res, res);
    realtype tmp = rl1*zn1[i] + ycor[i];
    res[i] = ngamma*ftemp[i] + tmp;
  }
}

extern "C"
int cvNlsResid_fused(const realtype rl1,
                     const realtype ngamma,
                     const N_Vector zn1,
                     const N_Vector ycor,
                     const N_Vector ftemp,
                     N_Vector res)
{
  const SUNExecPolicy* exec_policy = ((NVectorContent)res->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(res);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvNlsResid_kernel<<<grid, block, 0, *(exec_policy->stream())>>>
  (
    N,
    rl1,
    ngamma,
    N_VGetDeviceArrayPointer(zn1),
    N_VGetDeviceArrayPointer(ycor),
    N_VGetDeviceArrayPointer(ftemp),
    N_VGetDeviceArrayPointer(res)
  );

#ifdef SUNDIALS_DEBUG_GPU_LASTERROR
  gpuDeviceSynchronize();
  if (!gpuAssert(gpuGetLastError(), __FILE__, __LINE__)) return -1;
#endif

  return 0;
}

/*
 * -----------------------------------------------------------------
 * Form y with perturbation = FRACT*(func. iter. correction)
 * -----------------------------------------------------------------
 */

__global__
void cvDiagSetup_formY_kernel(const sunindextype length,
                              const realtype h,
                              const realtype r,
                              const realtype* fpred,
                              const realtype* zn1,
                              const realtype* ypred,
                              realtype* ftemp,
                              realtype* y)
{
  // N_VLinearSum(h, fpred, -ONE, zn[1], ftemp);
  // N_VLinearSum(r, ftemp, ONE, ypred, y);
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    ftemp[i] = h*fpred[i] - zn1[i];
    y[i] = r*ftemp[i] + ypred[i];
  }
}

extern "C"
int cvDiagSetup_formY(const realtype h,
                      const realtype r,
                      const N_Vector fpred,
                      const N_Vector zn1,
                      const N_Vector ypred,
                      N_Vector ftemp,
                      N_Vector y)
{
  const SUNExecPolicy* exec_policy = ((NVectorContent)y->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(y);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvDiagSetup_formY_kernel<<<grid, block, 0, *(exec_policy->stream())>>>
  (
    N,
    h,
    r,
    N_VGetDeviceArrayPointer(fpred),
    N_VGetDeviceArrayPointer(zn1),
    N_VGetDeviceArrayPointer(ypred),
    N_VGetDeviceArrayPointer(ftemp),
    N_VGetDeviceArrayPointer(y)
  );

#ifdef SUNDIALS_DEBUG_GPU_LASTERROR
  gpuDeviceSynchronize();
  if (!gpuAssert(gpuGetLastError(), __FILE__, __LINE__)) return -1;
#endif

  return 0;
}

/*
 * -----------------------------------------------------------------
 * Construct M = I - gamma*J with J = diag(deltaf_i/deltay_i)
 * protecting against deltay_i being at roundoff level.
 * -----------------------------------------------------------------
 */

__global__
void cvDiagSetup_buildM_kernel(const sunindextype length,
                               const realtype fract,
                               const realtype uround,
                               const realtype h,
                               const realtype* ftemp,
                               const realtype* fpred,
                               const realtype* ewt,
                               realtype* bit,
                               realtype* bitcomp,
                               realtype* y,
                               realtype* M)
{
  static const realtype zero = 0.0;
  static const realtype one = 1.0;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VLinearSum(ONE, M, -ONE, fpred, M);
    // N_VLinearSum(FRACT, ftemp, -h, M, M);
    // N_VProd(ftemp, ewt, y);
    M[i] = fract*ftemp[i] - h*(M[i] - fpred[i]);
    y[i] = ftemp[i] * ewt[i];

    // N_VCompare(uround, y, bit);
    // N_VAddConst(bit, -ONE, bitcomp);
    bool test = (abs(y[i]) > uround);
    bit[i] = test ? one : zero;
    bitcomp[i] = test ? zero : -one;

    // N_VProd(ftemp, bit, y);
    // N_VLinearSum(FRACT, y, -ONE, bitcomp, y);
    // N_VDiv(M, y, M);
    // N_VProd(M, bit, M);
    // N_VLinearSum(ONE, M, -ONE, bitcomp, M);
    y[i] = fract*ftemp[i]*bit[i] - bitcomp[i];
    M[i] = M[i]/y[i] * bit[i] - bitcomp[i];
  }
}

extern "C"
int cvDiagSetup_buildM(const realtype fract,
                       const realtype uround,
                       const realtype h,
                       const N_Vector ftemp,
                       const N_Vector fpred,
                       const N_Vector ewt,
                       N_Vector bit,
                       N_Vector bitcomp,
                       N_Vector y,
                       N_Vector M)
{
  const SUNExecPolicy* exec_policy = ((NVectorContent)M->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(M);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvDiagSetup_buildM_kernel<<<grid, block, 0, *(exec_policy->stream())>>>
  (
    N,
    fract,
    uround,
    h,
    N_VGetDeviceArrayPointer(ftemp),
    N_VGetDeviceArrayPointer(fpred),
    N_VGetDeviceArrayPointer(ewt),
    N_VGetDeviceArrayPointer(bit),
    N_VGetDeviceArrayPointer(bitcomp),
    N_VGetDeviceArrayPointer(y),
    N_VGetDeviceArrayPointer(M)
  );

#ifdef SUNDIALS_DEBUG_GPU_LASTERROR
  gpuDeviceSynchronize();
  if (!gpuAssert(gpuGetLastError(), __FILE__, __LINE__)) return -1;
#endif

  return 0;
}


/*
 * -----------------------------------------------------------------
 *  Update M with changed gamma so that M = I - gamma*J.
 * -----------------------------------------------------------------
 */


 __global__
void cvDiagSolve_updateM_kernel(const sunindextype length, const realtype r, realtype* M)
{
  static const realtype one = 1.0;
  GRID_STRIDE_XLOOP(sunindextype, i, length)
  {
    // N_VInv(M, M);
    // N_VAddConst(M, -ONE, M);
    // N_VScale(r, M, M);
    // N_VAddConst(M, ONE, M);
    realtype a = one/M[i] - one;
    M[i] = r*a + one;
  }
}


extern "C"
int cvDiagSolve_updateM(const realtype r, N_Vector M)
{
  const SUNExecPolicy* exec_policy = ((NVectorContent)M->content)->stream_exec_policy;
  const sunindextype N = N_VGetLength(M);
  size_t block = exec_policy->blockSize(N);
  size_t grid  = exec_policy->gridSize(N);

  cvDiagSolve_updateM_kernel<<<grid, block, 0, *(exec_policy->stream())>>>
  (
    N,
    r,
    N_VGetDeviceArrayPointer(M)
  );

#ifdef SUNDIALS_DEBUG_GPU_LASTERROR
  gpuDeviceSynchronize();
  if (!gpuAssert(gpuGetLastError(), __FILE__, __LINE__)) return -1;
#endif

  return 0;
}
