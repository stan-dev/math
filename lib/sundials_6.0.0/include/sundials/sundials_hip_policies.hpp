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
 * This header files defines the ExecPolicy classes which
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
namespace hip
{

class ExecPolicy
{
public:
  virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const = 0;
  virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const = 0;
  virtual const hipStream_t* stream() const = 0;
  virtual ExecPolicy* clone() const = 0;
  virtual ~ExecPolicy() {}
};


/*
 * A kernel execution policy that maps each thread to a work unit.
 * The number of threads per block (blockSize) can be set to anything.
 * The grid size will be chosen so that there are enough threads for one
 * thread per element. If a stream is provided, it will be used to
 * execute the kernel.
 */
class ThreadDirectExecPolicy : public ExecPolicy
{
public:
  ThreadDirectExecPolicy(const size_t blockDim, const hipStream_t stream = 0)
    : blockDim_(blockDim), stream_(stream)
  {}

  ThreadDirectExecPolicy(const ThreadDirectExecPolicy& ex)
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

  virtual ExecPolicy* clone() const
  {
    return static_cast<ExecPolicy*>(new ThreadDirectExecPolicy(*this));
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
class GridStrideExecPolicy : public ExecPolicy
{
public:
  GridStrideExecPolicy(const size_t blockDim, const size_t gridDim, const hipStream_t stream = 0)
    : blockDim_(blockDim), gridDim_(gridDim), stream_(stream)
  {}

  GridStrideExecPolicy(const GridStrideExecPolicy& ex)
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

  virtual ExecPolicy* clone() const
  {
    return static_cast<ExecPolicy*>(new GridStrideExecPolicy(*this));
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
class BlockReduceExecPolicy : public ExecPolicy
{
public:
  BlockReduceExecPolicy(const size_t blockDim, const size_t gridDim = 0, const hipStream_t stream = 0)
    : blockDim_(blockDim), gridDim_(gridDim), stream_(stream)
  {
    if (blockDim < 1 || blockDim % warpSize)
    {
      throw std::invalid_argument("the block size must be a multiple of the HIP warp size");
    }
  }

  BlockReduceExecPolicy(const BlockReduceExecPolicy& ex)
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

  virtual ExecPolicy* clone() const
  {
    return static_cast<ExecPolicy*>(new BlockReduceExecPolicy(*this));
  }

private:
  const hipStream_t stream_;
  const size_t blockDim_;
  const size_t gridDim_;
};

} // namespace hip
} // namespace sundials

typedef sundials::hip::ExecPolicy SUNHipExecPolicy;
typedef sundials::hip::ThreadDirectExecPolicy SUNHipThreadDirectExecPolicy;
typedef sundials::hip::GridStrideExecPolicy SUNHipGridStrideExecPolicy;
typedef sundials::hip::BlockReduceExecPolicy SUNHipBlockReduceExecPolicy;

#endif
