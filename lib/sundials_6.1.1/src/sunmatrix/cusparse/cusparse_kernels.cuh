/*
 * -----------------------------------------------------------------
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
 * This is the header file is for the cuSPARSE implementation of the
 * SUNMATRIX module.
 * -----------------------------------------------------------------
 */


#ifndef _SUNCUSPARSE_MATRIX_KERNELS_CUH_
#define _SUNCUSPARSE_MATRIX_KERNELS_CUH_

#include <cuda_runtime.h>
#include <sunmatrix/sunmatrix_cusparse.h>

namespace sundials
{
namespace sunmatrix_cusparse
{

template <typename T, typename I>
__global__ void
scaleAddIKernelCSR(I m, T c, T* A, const I* rowptr, const I* colind)
{
    // REQUIRES THE DIAGONAL TO BE PRESENT!

    // Each thread loops over one row of the matrix so memory accesses by a thread are stride-1.
    // If there aren't enough threads to cover all rows, then some threads will be reused for
    // more than one row.
    for (I row = blockIdx.x*blockDim.x + threadIdx.x;
         row < m;
         row += blockDim.x * gridDim.x)
    {
        I tmp = rowptr[row];
        I rownnz = rowptr[row+1] - tmp;
        I idx = tmp;
        for (I j = 0; j < rownnz; j++)
        {
            if (colind[idx+j] == row) A[idx+j] = c*A[idx+j] + 1.0;
            else                      A[idx+j] = c*A[idx+j];
        }
    }
}

template <typename T, typename I>
__global__ void
scaleAddIKernelBCSR(I m, I nblocks, I blocknnz, T c, T* A, const I* rowptr, const I* colind)
{
    // REQUIRES THE DIAGONAL TO BE PRESENT!

    // Ideally each thread block will be in charge of one block of the matrix.
    for (I block = blockIdx.x;
         block < nblocks;
         block += gridDim.x)
    {
        // Each thread loops over one row of the matrix so memory accesses by a thread are stride-1.
        // If there aren't enough threads to cover all rows, then some threads will be reused for
        // more than one row.
        for (I row = threadIdx.x;
             row < m;
             row += blockDim.x)
        {
            I tmp = rowptr[row];
            I rownnz = rowptr[row+1] - tmp;
            I idxl = tmp;
            I idxg = block*blocknnz + tmp;
            for (I j = 0; j < rownnz; j++)
            {
                if (colind[idxl+j] == row) A[idxg+j] = c*A[idxg+j] + 1.0;
                else                       A[idxg+j] = c*A[idxg+j];
            }
        }
    }
}

template <typename T, typename I>
__global__ void
scaleAddKernelCSR(I nnz, T c, T* A, const T* B)
{
    // REQUIRES A AND B TO HAVE THE SAME SPARSITY PATTERN
    for (I i = blockIdx.x * blockDim.x + threadIdx.x;
         i < nnz;
         i += blockDim.x * gridDim.x)
    {
        A[i] = c*A[i] + B[i];
    }
}

template <typename T, typename I>
__global__ void
matvecBCSR(I m, I nblocks, I blocknnz, const T* A, const I* rowptr, const I* colind, const T* x,  T* y)
{
    // Zero out result vector
    for (I i = blockIdx.x * blockDim.x + threadIdx.x;
         i < nblocks*blocknnz;
         i += blockDim.x * gridDim.x)
    {
        y[i] = 0.0;
    }

    __syncthreads();

    // Ideally each thread block will be in charge of one block of the matrix.
    for (I block = blockIdx.x;
         block < nblocks;
         block += gridDim.x)
    {
        // Each thread loops over one row of the matrix so memory accesses by a thread are stride-1.
        // If there aren't enough threads to cover all rows, then some threads will be reused for
        // more than one row.
        for (I row = threadIdx.x;
             row < m;
             row += blockDim.x)
        {
            I tmp = rowptr[row];
            I rownnz = rowptr[row+1] - tmp; // number of nnz in this row
            I idxl = tmp;                   // local (to this block) starting nonzero index
            I idxg = block*blocknnz + tmp;  // global (overall matrix) starting nonzero index
            I rowg = block*m+row;           // global (overall matrix) row
            I colg = block*m;               // global (overall matrix) starting column
            for (I j = 0; j < rownnz; j++)
            {
                y[rowg] += A[idxg+j] * x[ colg+colind[idxl+j] ];
            }
        }
    }
}

// kernels for debugging
#ifdef SUNDIALS_DEBUG

template <typename T, typename I>
__global__ void
print_kernel(I m, I nnz, I blocknnz, T* A, const I* rowptr, const I* colind)
{
    for (I i = blockIdx.x * blockDim.x + threadIdx.x;
         i < nnz;
         i += blockDim.x * gridDim.x)
    {
        printf("A[%d] = %f\n", i, A[i]);
    }
    for (I i = blockIdx.x * blockDim.x + threadIdx.x;
        i < m+1;
        i += blockDim.x * gridDim.x)
    {
        printf("rowptr[%d] = %d\n", i, rowptr[i]);
    }
    for (I i = blockIdx.x * blockDim.x + threadIdx.x;
        i < blocknnz;
        i += blockDim.x * gridDim.x)
    {
        printf("colind[%d] = %d\n", i, colind[i]);
    }
}

#endif

} // namespace sunmatrix_cusparse
} // namespace sundials

#endif