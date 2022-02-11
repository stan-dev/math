/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Cody J. Balos @ LLNL
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
 */

#ifndef _NVECTOR_HIP_KERNELS_HIP_HPP_
#define _NVECTOR_HIP_KERNELS_HIP_HPP_

#include <limits>
#include <hip/hip_runtime.h>

#include "sundials_hip_kernels.hip.hpp"

using namespace sundials::hip;

namespace sundials
{
namespace nvector_hip
{

/* -----------------------------------------------------------------
 * The namespace for HIP kernels
 *
 * Reduction HIP kernels in nvector are based in part on "reduction"
 * example in NVIDIA Corporation CUDA Samples, and parallel reduction
 * examples in textbook by J. Cheng at al. "CUDA C Programming".
 * -----------------------------------------------------------------
 */

/*
 * Sets all elements of the vector X to constant value a.
 *
 */

template <typename T, typename I>
__global__ void
setConstKernel(T a, T *X, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    X[i] = a;
  }
}


/*
 * Computes linear sum (combination) of two vectors.
 *
 */

template <typename T, typename I>
__global__ void
linearSumKernel(T a, const T *X, T b, const T *Y, T *Z, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    Z[i] = a*X[i] + b*Y[i];
  }
}


/*
 * Elementwise product of two vectors.
 *
 */

template <typename T, typename I>
__global__ void
prodKernel(const T *X, const T *Y, T *Z, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    Z[i] = X[i]*Y[i];
  }
}


/*
 * Elementwise division of two vectors.
 *
 */

template <typename T, typename I>
__global__ void
divKernel(const T *X, const T *Y, T *Z, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    Z[i] = X[i]/Y[i];
  }
}


/*
 * Scale vector with scalar value 'a'.
 *
 */

template <typename T, typename I>
__global__ void
scaleKernel(T a, const T *X, T *Z, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    Z[i] = a*X[i];
  }
}


/*
 * Stores absolute values of vector X elements into vector Z.
 *
 */

template <typename T, typename I>
__global__ void
absKernel(const T *X, T *Z, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    Z[i] = abs(X[i]);
  }
}


/*
 * Elementwise inversion.
 *
 */

template <typename T, typename I>
__global__ void
invKernel(const T *X, T *Z, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    Z[i] = 1.0/(X[i]);
  }
}


/*
 * Add constant 'c' to each vector element.
 *
 */

template <typename T, typename I>
__global__ void
addConstKernel(T a, const T *X, T *Z, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    Z[i] = a + X[i];
  }
}


/*
 * Compare absolute values of vector 'X' with constant 'c'.
 *
 */

template <typename T, typename I>
__global__ void
compareKernel(T c, const T *X, T *Z, I n)
{
  GRID_STRIDE_XLOOP(I, i, n)
  {
    Z[i] = (abs(X[i]) >= c) ? 1.0 : 0.0;
  }
}


/*
 * Dot product of two vectors.
 *
 */
template <typename T, typename I>
__global__ void
dotProdKernel(const T *x, const T *y, T *out, I n)
{
  T sum = 0.0;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    sum += x[i] * y[i];
  }
  sum = blockReduce<T, RSUM>(sum, 0.0);

  // Copy reduction result for each block to global memory
  if (threadIdx.x == 0) atomicAdd(out, sum);
}


/*
 * Finds max norm the vector.
 *
 */
template <typename T, typename I>
__global__ void
maxNormKernel(const T *x, T *out, I n)
{
  T maximum = 0.0;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    maximum = max(abs(x[i]), maximum);
  }
  maximum = blockReduce<T, RMAX>(maximum, 0.0);

  // Maximum of reduction result for each block
  if (threadIdx.x == 0) AtomicMax(out, maximum);
}


/*
 * Weighted L2 norm squared.
 *
 */
template <typename T, typename I>
__global__ void
wL2NormSquareKernel(const T *x, const T *w, T *out, I n)
{
  T sum = 0.0;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    sum += x[i] * w[i] * x[i] * w[i];
  }
  sum = blockReduce<T, RSUM>(sum, 0.0);

  // Copy reduction result for each block to global memory
  if (threadIdx.x == 0) atomicAdd(out, sum);
}

/*
 * Weighted L2 norm squared with mask. Vector id specifies the mask.
 *
 */
template <typename T, typename I>
__global__ void
wL2NormSquareMaskKernel(const T *x, const T *w, const T *id, T *out, I n)
{
  T sum = 0.0;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    if(id[i] > 0.0) sum += x[i] * w[i] * x[i] * w[i];
  }
  sum = blockReduce<T, RSUM>(sum, 0.0);

  // Copy reduction result for each block to global memory
  if (threadIdx.x == 0) atomicAdd(out, sum);
}


/*
 * Finds min value in the vector.
 *
 */
template <typename T, typename I>
__global__ void
findMinKernel(T MAX_VAL, const T *x, T *out, I n)
{
  T minimum = MAX_VAL;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    minimum = min(x[i], minimum);
  }
  minimum = blockReduce<T, RMIN>(minimum, MAX_VAL);

  // minimum of reduction result for each block
  if (threadIdx.x == 0) AtomicMin(out, minimum);
}


/*
 * Computes L1 norm of vector
 *
 */
template <typename T, typename I>
__global__ void
L1NormKernel(const T *x, T *out, I n)
{
  T sum = 0.0;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    sum += abs(x[i]);
  }
  sum = blockReduce<T, RSUM>(sum, 0.0);

  // Copy reduction result for each block to global memory
  if (threadIdx.x == 0) atomicAdd(out, sum);
}

/*
 * Vector inverse  z[i] = 1/x[i] with check for zeros. Reduction is performed
 * to flag the result if any x[i] = 0.
 *
 */
template <typename T, typename I>
__global__ void
invTestKernel(const T *x, T *z, T *out, I n)
{
  T flag = 0.0;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    if (x[i] == static_cast<T>(0.0))
      flag += 1.0;
    else
      z[i] = 1.0/x[i];
  }
  flag = blockReduce<T, RSUM>(flag, 0.0);

  // Copy reduction result for each block to global memory
  if (threadIdx.x == 0) atomicAdd(out, flag);
}


/*
 * Checks if inequality constraints are satisfied. Constraint check
 * results are stored in vector 'm'. A sum reduction over all elements
 * of 'm' is performed to find if any of the constraints is violated.
 * If all constraints are satisfied sum == 0.
 *
 */
template <typename T, typename I>
__global__ void
constrMaskKernel(const T *c, const T *x, T *m, T *out, I n)
{
  T sum = 0.0;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    // test = true if constraints violated
    bool test = (std::abs(c[i]) > 1.5 && c[i]*x[i] <= 0.0) ||
                (std::abs(c[i]) > 0.5 && c[i]*x[i] <  0.0);
    m[i] = test ? 1.0 : 0.0;
    sum = m[i];
  }
  sum = blockReduce<T, RSUM>(sum, 0.0);

  // Copy reduction result for each block to global memory
  if (threadIdx.x == 0) atomicAdd(out, sum);
}


/*
 * Finds minimum component-wise quotient.
 *
 */
template <typename T, typename I>
__global__ void
minQuotientKernel(const T MAX_VAL, const T *num, const T *den, T *min_quotient, I n)
{
  T minimum = MAX_VAL;
  T quotient = 0.0;
  GRID_STRIDE_XLOOP(I, i, n)
  {
    quotient = (den[i] == static_cast<T>(0.0)) ? MAX_VAL : num[i]/den[i];
    minimum = min(quotient, minimum);
  }
  minimum = blockReduce<T, RMIN>(minimum, MAX_VAL);

  // minimum of reduction result for each block
  if (threadIdx.x == 0) AtomicMin(min_quotient, minimum);
}

} // namespace nvector_hip
} // namespace sundials

#endif // _NVECTOR_HIP_KERNELS_HIP_HPP_
