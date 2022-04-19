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
 */

#ifndef _SUNDIALS_HIP_KERNELS_CUH
#define _SUNDIALS_HIP_KERNELS_CUH

#define SUNDIALS_HOST_DEVICE __host__ __device__
#define SUNDIALS_DEVICE_INLINE __forceinline__
#include "sundials_reductions.hpp"

#define GRID_STRIDE_XLOOP(type, iter, max)  \
  for (type iter = blockDim.x * blockIdx.x + threadIdx.x; \
       iter < max; \
       iter += blockDim.x * gridDim.x)

#include "sundials_hip.h"

namespace sundials
{
namespace hip
{
namespace impl
{

template <typename T>
__forceinline__ __device__ T shfl_xor_sync(T var, int laneMask);

template <typename T>
__forceinline__ __device__ T shfl_sync(T var, int srcLane);

template <typename T>
__forceinline__ __device__ T shfl_down_sync(T var, int srcLane);

template <>
__forceinline__ __device__ int shfl_xor_sync<int>(int var, int laneMask)
{
  return ::__shfl_xor(var, laneMask);
}

template <>
__forceinline__ __device__ float shfl_xor_sync<float>(float var, int laneMask)
{
  return ::__shfl_xor(var, laneMask);
}

template <>
__forceinline__ __device__ double shfl_xor_sync<double>(double var, int laneMask)
{
  return ::__shfl_xor(var, laneMask);
}

template <>
__forceinline__ __device__ int shfl_sync<int>(int var, int srcLane)
{
  return ::__shfl(var, srcLane);
}

template <>
__forceinline__ __device__ float shfl_sync<float>(float var, int srcLane)
{
  return ::__shfl(var, srcLane);
}

template <>
__forceinline__ __device__ double shfl_sync<double>(double var, int srcLane)
{
  return ::__shfl(var, srcLane);
}

template<>
__forceinline__ __device__ float shfl_down_sync(float val, int srcLane)
{
  return ::__shfl_down(val, srcLane);
}

template<>
__forceinline__ __device__ double shfl_down_sync(double val, int srcLane)
{
  return ::__shfl_down(val, srcLane);
}


/* The atomic functions below are implemented using the atomic compare and swap
   function atomicCAS which performs an atomic version of
   (*address == assumed) ? (assumed + val) : *address. Since *address could change
   between when the value is loaded and the atomicCAS call the operation is repeated
   until *address does not change between the read and the compare and swap operation. */

__forceinline__ __device__
double atomicAdd(double* address, double val)
{
#ifndef __HIP_ARCH_HAS_DOUBLE_ATOMIC_ADD__
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
#else
  return ::atomicAdd(address, val);
#endif
}

__forceinline__ __device__
float atomicAdd(float* address, float val)
{
#ifndef __HIP_ARCH_HAS_FLOAT_ATOMIC_ADD__
  unsigned int* address_as_ull = (unsigned int*)address;
  unsigned int old = *address_as_ull, assumed;

  do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
                      __float_as_int(val +
                              __int_as_float(assumed)));
  // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);

  return __int_as_float(old);
#else
  return ::atomicAdd(address, val);
#endif
}

/*
 * Compute the maximum of 2 double-precision floating point values using an atomic operation
 * "address" is the address of the reference value which might get updated with the maximum
 * "value" is the value that is compared to the reference in order to determine the maximum
 */
__forceinline__ __device__
void atomicMax(double* const address, const double value)
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
void atomicMax(float* const address, const float value)
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
void atomicMin(double* const address, const double value)
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
void atomicMin(float* const address, const float value)
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

//
// Atomic specializations of sundials::reductions operators
//

template<typename BinaryReductionOp>
struct atomic;

template<typename T>
struct atomic<sundials::reductions::impl::plus<T>> {
  __device__ __forceinline__ void operator()(T* out, const T val)
  {
    atomicAdd(out, val);
  }
};

template<typename T>
struct atomic<sundials::reductions::impl::maximum<T>> {
  __device__ __forceinline__ void operator()(T* out, const T val)
  {
    atomicMax(out, val);
  }
};

template<typename T>
struct atomic<sundials::reductions::impl::minimum<T>> {
  __device__ __forceinline__ void operator()(T* out, const T val)
  {
    atomicMin(out, val);
  }
};


/*
 * Perform a reduce on the warp to get the operation result.
 */
template <typename T, typename BinaryReductionOp>
__inline__ __device__
T warpReduceShflDown(T val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
  {
    T rhs = shfl_down_sync<T>(val, offset);
    val = BinaryReductionOp{}(val, rhs);
  }
  return val;
}

/*
 * Reduce value across the thread block.
 */
template <typename T, typename op>
__inline__ __device__
T blockReduceShflDown(T val, T identity)
{
  // Shared memory for the partial sums
  static __shared__ T shared[MAX_WARPS];

  int numThreads = blockDim.x * blockDim.y * blockDim.z;

  int threadId = threadIdx.x + blockDim.x * threadIdx.y +
                 (blockDim.x * blockDim.y) * threadIdx.z;

  int warpId = threadId / WARP_SIZE;
  int warpLane = threadId % WARP_SIZE;

  // Each warp performs partial reduction
  val = warpReduceShflDown<T, op>(val);

  // Write reduced value from each warp to shared memory
  if (warpLane == 0) shared[warpId] = val;

  // Wait for all partial reductions to complete
  __syncthreads();

  // Read per warp values from shared memory only if that warp existed
  val = (threadId < numThreads / warpSize) ? shared[warpLane] : identity;

  // Final reduce within first warp
  if (warpId == 0)
    val = warpReduceShflDown<T, op>(val);

  return val;
}

/*
 * Warp reduce + block reduce using shfl instead of shfl_down.
 */
template <typename T, typename BinaryReductionOp>
__inline__ __device__
T blockReduceShfl(T val, T identity)
{
  int numThreads = blockDim.x * blockDim.y * blockDim.z;

  int threadId = threadIdx.x + blockDim.x * threadIdx.y +
                 (blockDim.x * blockDim.y) * threadIdx.z;

  int warpId = threadId / WARP_SIZE;
  int warpLane = threadId % WARP_SIZE;

  T temp = val;

  // Reduce each warp
  if (numThreads % WARP_SIZE == 0)
  {
    for (int i = 1; i < WARP_SIZE; i *= 2)
    {
      T rhs = shfl_xor_sync<T>(temp, i);
      temp = BinaryReductionOp{}(temp, rhs);
    }
  }
  else
  {
    for (int i = 1; i < WARP_SIZE; i *= 2)
    {
      int srcLane = threadId ^ i;
      T rhs = shfl_sync<T>(temp, srcLane);
      // Only add from threads that exist to avoid double counting
      if (srcLane < numThreads)
        temp = BinaryReductionOp{}(temp, rhs);
    }
  }

  // Reduce per warp values
  if (numThreads > WARP_SIZE)
  {
    static_assert(MAX_WARPS <= WARP_SIZE, "max warps must be <= warp size for this algorithm to work");

    __shared__ T shared[MAX_WARPS];

    // Write per warp values to shared memory
    if (warpLane == 0)
      shared[warpId] = temp;

    __syncthreads();

    if (warpId == 0)
    {
      // Read per warp values only if the warp existed
      temp = (warpLane * WARP_SIZE < numThreads) ? shared[warpLane] : identity;

      // Final reduction
      for (int i = 1; i < MAX_WARPS; i *= 2)
      {
        T rhs = shfl_xor_sync<T>(temp, i);
        temp = BinaryReductionOp{}(temp, rhs);
      }
    }

    __syncthreads();
  }

  return temp;
}

/*
 * Reduce values into thread 0 of the last running thread block.
 * Output value is device_mem[0].
 */
template <typename T, typename BinaryReductionOp>
__device__ __forceinline__ void gridReduce(T val,
                                           T identity,
                                           T* device_mem,
                                           unsigned int* device_count)
{
  int numBlocks = gridDim.x * gridDim.y * gridDim.z;
  int numThreads = blockDim.x * blockDim.y * blockDim.z;
  unsigned int wrap_around = numBlocks - 1;

  int blockId = blockIdx.x + gridDim.x * blockIdx.y +
                (gridDim.x * gridDim.y) * blockIdx.z;

  int threadId = threadIdx.x + blockDim.x * threadIdx.y +
                 (blockDim.x * blockDim.y) * threadIdx.z;

  // Each block reduces a subset of the input
  T temp = blockReduceShfl<T, BinaryReductionOp>(val, identity);

  __shared__ bool isLastBlockDone;
  if (threadId == 0)
  {
    // One thread per block stores the partial reductions to global memory
    device_mem[blockId] = temp;

    // Ensure write visible to all threads
    __threadfence();

    // Increment counter, (wraps back to zero if old count == wrap_around)
    unsigned int old_count = atomicInc(device_count, wrap_around);
    isLastBlockDone = (old_count == wrap_around) ? 1 : 0;
  }

  // Synchronize to ensure that each thread reads the
  // correct value of isLastBlockDone.
  __syncthreads();

  // The last block reduces values in device_mem
  if (isLastBlockDone)
  {
    // Reduce thread_i in each block into temp
    temp = identity;
    for (int i = threadId; i < numBlocks; i += numThreads)
      temp = BinaryReductionOp{}(temp, device_mem[i]);

    // Compute the final block partial reductions
    temp = blockReduceShfl<T, BinaryReductionOp>(temp, identity);

    // One thread returns the final value
    if (threadId == 0)
      device_mem[0] = temp;
  }
}

template<typename T, typename BinaryReductionOp>
__device__ __forceinline__ void gridReduceAtomic(T val,
                                                 T identity,
                                                 T* device_mem)
{
  int threadId = threadIdx.x + blockDim.x * threadIdx.y +
                 (blockDim.x * blockDim.y) * threadIdx.z;
  val = blockReduceShflDown<T, BinaryReductionOp>(val, identity);
  // Final reduction of all block values into the output device_mem
  if (threadId == 0)
    atomic<BinaryReductionOp>{}(device_mem, val);
}

template<typename T, typename BinaryReductionOp>
struct GridReducerLDS
{
  __device__ __forceinline__ void operator()(T val,
                                             T identity,
                                             T* device_mem,
                                             unsigned int* device_count)
  {
    return impl::gridReduce<T, BinaryReductionOp>(val, identity, device_mem, device_count);
  }
};


template<typename T, typename BinaryReductionOp>
struct GridReducerAtomic
{
  __device__ __forceinline__ void operator()(T val,
                                             T identity,
                                             T* device_mem,
                                             unsigned int* device_count)
  {
    return impl::gridReduceAtomic<T, BinaryReductionOp>(val, identity, device_mem);
  }
};

} // namespace impl
} // namespace hip
} // namespace sundials

#endif // _SUNDIALS_HIP_KERNELS_CUH
