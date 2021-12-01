/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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

#ifndef _SUNDIALS_HIP_KERNELS_HIP_HPP
#define _SUNDIALS_HIP_KERNELS_HIP_HPP

#include "sundials_hip.h"

#define GRID_STRIDE_XLOOP(type, iter, max)  \
  for (type iter = blockDim.x * blockIdx.x + threadIdx.x; \
       iter < max; \
       iter += blockDim.x * gridDim.x)

/* Uses correct __shfl_down depending on architecture being used */
#if defined(__CUDA_ARCH__)
  #define _SHFL_DOWN(val,offset) (__shfl_down_sync(0xFFFFFFFF, val, offset))
#else
  #define _SHFL_DOWN(val,offset) (__shfl_down(val, offset))
#endif

namespace sundials
{
namespace hip
{

/* The atomic functions below are implemented using the atomic compare and swap
   function atomicCAS which performs an atomic version of
   (*address == assumed) ? (assumed + val) : *address. Since *address could change
   between when the value is loaded and the atomicCAS call the operation is repeated
   until *address does not change between the read and the compare and swap operation. */

typedef enum { RSUM, RMAX, RMIN } BinaryReductionOp;

#if defined(__CUDA_ARCH__) and __CUDA_ARCH__ < 600
__forceinline__ __device__
double atomicAdd(double* address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
                      __double_as_longlong(val +
                              __longlong_as_double(assumed)));
  // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);

  return __longlong_as_double(old);
}
#endif

/*
 * Compute the maximum of 2 double-precision floating point values using an atomic operation
 * "address" is the address of the reference value which might get updated with the maximum
 * "value" is the value that is compared to the reference in order to determine the maximum
 */
__forceinline__ __device__
void AtomicMax(double* const address, const double value)
{
  if (*address >= value)
  {
    return;
  }

  unsigned long long * const address_as_i = (unsigned long long *)address;
  unsigned long long old = * address_as_i, assumed;

  do
  {
    assumed = old;
    if (__longlong_as_double(assumed) >= value)
    {
      break;
    }
    old = atomicCAS(address_as_i, assumed, __double_as_longlong(value));
  } while (assumed != old);
}

/*
 * Compute the maximum of 2 single-precision floating point values using an atomic operation
 * "address" is the address of the reference value which might get updated with the maximum
 * "value" is the value that is compared to the reference in order to determine the maximum
 */
 __forceinline__ __device__
void AtomicMax(float* const address, const float value)
{
  if (*address >= value)
  {
    return;
  }

  unsigned int* const address_as_i = (unsigned int *)address;
  unsigned int old = *address_as_i, assumed;

  do
  {
    assumed = old;
    if (__int_as_float(assumed) >= value)
    {
      break;
    }
    old = atomicCAS(address_as_i, assumed, __float_as_int(value));
  } while (assumed != old);
}

/*
 * Compute the minimum of 2 double-precision floating point values using an atomic operation
 * "address" is the address of the reference value which might get updated with the minimum
 * "value" is the value that is compared to the reference in order to determine the minimum
 */
__forceinline__ __device__
void AtomicMin(double* const address, const double value)
{
  if (*address <= value)
  {
    return;
  }

  unsigned long long* const address_as_i = (unsigned long long *)address;
  unsigned long long old = *address_as_i, assumed;

  do
  {
    assumed = old;
    if (__longlong_as_double(assumed) <= value)
    {
      break;
    }
    old = atomicCAS(address_as_i, assumed, __double_as_longlong(value));
  } while (assumed != old);
}

/*
 * Compute the minimum of 2 single-precision floating point values using an atomic operation
 * "address" is the address of the reference value which might get updated with the minimum
 * "value" is the value that is compared to the reference in order to determine the minimum
 */
__forceinline__ __device__
void AtomicMin(float* const address, const float value)
{
  if (*address <= value)
  {
    return;
  }

  unsigned int* const address_as_i = (unsigned int *)address;
  unsigned int old = *address_as_i, assumed;

  do
  {
    assumed = old;
    if (__int_as_float(assumed) <= value)
    {
      break;
    }
    old = atomicCAS(address_as_i, assumed, __float_as_int(value));
  } while (assumed != old);
}

/*
 * Perform a reduce on the warp to get the sum.
 */
template <typename T>
__inline__ __device__
T warpReduceSum(T val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    val += _SHFL_DOWN(val, offset);
  return val;
}

/*
 * Perform a reduce on the warp to get the maximum value.
 */
template<typename T>
__inline__ __device__
T warpReduceMax(T val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    val = max(_SHFL_DOWN(val, offset), val);
  return val;
}

/*
 * Perform a reduce on the warp to get the minimum value.
 */
template<typename T>
__inline__ __device__
T warpReduceMin(T val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    val = min(_SHFL_DOWN(val, offset), val);
  return val;
}

/*
 * Reduce value across the thread block.
 */
template <typename T, BinaryReductionOp op>
__inline__ __device__
T blockReduce(T val, T default_value)
{
  // Shared memory for the partial sums
  static __shared__ T shared[warpSize];

  int lane = threadIdx.x % warpSize; // thread lane within warp
  int wid = threadIdx.x / warpSize;  // warp ID

  // Each warp performs partial reduction
  switch(op)
  {
    case RSUM:
      val = warpReduceSum<T>(val);
      break;
    case RMAX:
      val = warpReduceMax<T>(val);
      break;
    case RMIN:
      val = warpReduceMin<T>(val);
      break;
    default:
      asm("trap;"); // illegal instruction
      break;
  }

  // Write reduced value from each warp to shared memory
  if (lane == 0) shared[wid] = val;

  // Wait for all partial reductions to complete
  __syncthreads();

  // Read from shared memory only if that warp existed
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : default_value;

  // Final reduce within first warp
  if (wid == 0)
  {
    switch(op)
    {
      case RSUM:
        val = warpReduceSum<T>(val);
        break;
      case RMAX:
        val = warpReduceMax<T>(val);
        break;
      case RMIN:
        val = warpReduceMin<T>(val);
        break;
      default:
        asm("trap;"); // illegal instruction
        break;
    }
  }

  return val;
}

} // namespace hip
} // namespace sundials

#endif // _SUNDIALS_HIP_KERNELS_HIP_HPP
