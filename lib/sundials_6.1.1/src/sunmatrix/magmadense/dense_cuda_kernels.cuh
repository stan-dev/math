/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the implementation file for the dense matrix CUDA kernels
 * for the SUNMATRIX package based on MAGMA.
 * -----------------------------------------------------------------*/

#ifndef _SUNGPUDENSE_MATRIX_KERNELS_CUH_
#define _SUNGPUDENSE_MATRIX_KERNELS_CUH_

#include <cuda_runtime.h>
#include <sunmatrix/sunmatrix_magmadense.h>

namespace sundials
{
namespace sunmatrix_gpudense
{
namespace cuda
{

  template <typename T, typename I, typename Lambda>
  __device__ __forceinline__ void
  block_col_row(I nblocks, I m, I n, Lambda&& fn)
  {
    for (I block = blockIdx.x*blockDim.x + threadIdx.x;
         block < nblocks;
         block += blockDim.x*gridDim.x)
    {
      for (I col = blockIdx.y*blockDim.y + threadIdx.y;
           col < n;
           col += blockDim.y*gridDim.y)
      {
        for (I row = blockIdx.z*blockDim.z + threadIdx.z;
             row < m;
             row += blockDim.z*gridDim.z)
        {
          fn(block*m*n+(col*m + row), row, col);
        }
      }
    }
  }

  template <typename T, typename I>
  __global__ void
  getBlockPointers(I m, I n, I nblocks, T* A, T** Ablocks,
                   T* x, T** xblocks, T* y, T** yblocks)
  {
    for (I block = blockIdx.x*blockDim.x + threadIdx.x;
         block < nblocks;
         block += blockDim.x*gridDim.x)
    {
      Ablocks[block] = &A[block*m*n];
      xblocks[block] = &x[block*n];
      yblocks[block] = &y[block*m];
    };
  }

  template <typename T, typename I>
  __global__ void
  zeroKernel(I m, I n, I nblocks, T* A)
  {
    block_col_row<T,I>(nblocks, m, n,
      [=] __device__ (I kij, I row, I col) {
        A[kij] = 0.0;
      });
  }

  template <typename T, typename I>
  __global__ void
  copyKernel(I m, I n, I nblocks, const T* A, T* B)
  {
    block_col_row<T,I>(nblocks, m, n,
      [=] __device__ (I kij, I row, I col) {
        B[kij] = A[kij];
      });
  }

  template <typename T, typename I>
  __global__ void
  scaleAddIKernel(I m, I n, I nblocks, T c, T* A)
  {
    block_col_row<T,I>(nblocks, m, n,
      [=] __device__ (I kij, I row, I col) {
        if (row == col) A[kij] = c*A[kij] + 1.0;
        else            A[kij] = c*A[kij];
      });
  }

  template <typename T, typename I>
  __global__ void
  scaleAddKernel(I m, I n, I nblocks, T c, T* A, const T* B)
  {
    block_col_row<T,I>(nblocks, m, n,
      [=] __device__ (I kij, I row, I col) {
        A[kij] = c*A[kij] + B[kij];
      });
  }

} // namespace cuda
} // namespace sunmatrix_gpudense
} // namespace sundials

#endif
