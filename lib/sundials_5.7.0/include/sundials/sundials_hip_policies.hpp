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
 * This header files defines the HipExecPolicy classes which
 * are utilized to determine HIP kernel launch paramaters.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_HIPEXECPOLICIES_HPP
#define _SUNDIALS_HIPEXECPOLICIES_HPP

#include <cstdio>
#include <stdexcept>
#include <hip/hip_runtime.h>

namespace sundials
{

class HipExecPolicy
{
public:
  virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const = 0;
  virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const = 0;
  virtual const hipStream_t* stream() const = 0;
  virtual HipExecPolicy* clone() const = 0;
  virtual ~HipExecPolicy() {}
};


/*
 * A kernel execution policy that maps each thread to a work unit.
 * The number of threads per block (blockSize) can be set to anything.
 * The grid size will be chosen so that there are enough threads for one
 * thread per element. If a stream is provided, it will be used to
 * execute the kernel.
 */
class HipThreadDirectExecPolicy : public HipExecPolicy
{
public:
  HipThreadDirectExecPolicy(const size_t blockDim, const hipStream_t stream = 0)
    : blockDim_(blockDim), stream_(stream)
  {}

  HipThreadDirectExecPolicy(const HipThreadDirectExecPolicy& ex)
    : blockDim_(ex.blockDim_), stream_(ex.stream_)
  {}

  virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const
  {
    /* ceil(n/m) = floor((n + m - 1) / m) */
    return (numWorkUnits + blockSize() - 1) / blockSize();
  }

  virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const
  {
    return blockDim_;
  }

  virtual const hipStream_t* stream() const
  {
    return &stream_;
  }

  virtual HipExecPolicy* clone() const
  {
    return static_cast<HipExecPolicy*>(new HipThreadDirectExecPolicy(*this));
  }

private:
  const hipStream_t stream_;
  const size_t blockDim_;
};

/*
 * A kernel execution policy for kernels that use grid stride loops.
 * The number of threads per block (blockSize) can be set to anything.
 * The number of blocks (gridSize) can be set to anything. If a stream
 * is provided, it will be used to execute the kernel.
 */
class HipGridStrideExecPolicy : public HipExecPolicy
{
public:
  HipGridStrideExecPolicy(const size_t blockDim, const size_t gridDim, const hipStream_t stream = 0)
    : blockDim_(blockDim), gridDim_(gridDim), stream_(stream)
  {}

  HipGridStrideExecPolicy(const HipGridStrideExecPolicy& ex)
    : blockDim_(ex.blockDim_), gridDim_(ex.gridDim_), stream_(ex.stream_)
  {}

  virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const
  {
    return gridDim_;
  }

  virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const
  {
    return blockDim_;
  }

  virtual const hipStream_t* stream() const
  {
    return &stream_;
  }

  virtual HipExecPolicy* clone() const
  {
    return static_cast<HipExecPolicy*>(new HipGridStrideExecPolicy(*this));
  }

private:
  const hipStream_t stream_;
  const size_t blockDim_;
  const size_t gridDim_;
};


/*
 * A kernel execution policy for performing a reduction across indvidual thread
 * blocks. The number of threads per block (blockSize) can be set to any valid
 * multiple of the HIP warp size. The number of blocks (gridSize) can be set to
 * any value greater or equal to 0. If it is set to 0, then the grid size will
 * be chosen so that there are at most two work units per thread. If a stream is
 * provided, it will be used to execute the kernel.
 */
class HipBlockReduceExecPolicy : public HipExecPolicy
{
public:
  HipBlockReduceExecPolicy(const size_t blockDim, const size_t gridDim = 0, const hipStream_t stream = 0)
    : blockDim_(blockDim), gridDim_(gridDim), stream_(stream)
  {
    if (blockDim < 1 || blockDim % warpSize)
    {
      throw std::invalid_argument("the block size must be a multiple of the HIP warp size");
    }
  }

  HipBlockReduceExecPolicy(const HipBlockReduceExecPolicy& ex)
    : blockDim_(ex.blockDim_), gridDim_(ex.gridDim_), stream_(ex.stream_)
  {}

  virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const
  {
    if (gridDim_ == 0)
    {
      return (numWorkUnits + (blockSize() * 2 - 1)) / (blockSize() * 2);
    }
    return gridDim_;
  }

  virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const
  {
    return blockDim_;
  }

  virtual const hipStream_t* stream() const
  {
    return &stream_;
  }

  virtual HipExecPolicy* clone() const
  {
    return static_cast<HipExecPolicy*>(new HipBlockReduceExecPolicy(*this));
  }

private:
  const hipStream_t stream_;
  const size_t blockDim_;
  const size_t gridDim_;
};

} // namespace sundials

typedef sundials::HipExecPolicy SUNHipExecPolicy;
typedef sundials::HipThreadDirectExecPolicy SUNHipThreadDirectExecPolicy;
typedef sundials::HipGridStrideExecPolicy SUNHipGridStrideExecPolicy;
typedef sundials::HipBlockReduceExecPolicy SUNHipBlockReduceExecPolicy;

#endif
