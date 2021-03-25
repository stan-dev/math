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
 * This header files defines the SyclExecPolicy classes which
 * are utilized to determine SYCL kernel launch paramaters.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_SYCLEXECPOLICIES_HPP
#define _SUNDIALS_SYCLEXECPOLICIES_HPP

#include <cstdio>
#include <stdexcept>
#include <CL/sycl.hpp>

namespace sundials
{

class SyclExecPolicy
{
public:
  virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const = 0;
  virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const = 0;
  virtual SyclExecPolicy* clone() const = 0;
  virtual ~SyclExecPolicy() {}
};


/*
 * A kernel execution policy that maps each thread to a work unit.
 * The number of threads per block (blockSize) can be set to anything.
 * The grid size will be chosen so that there are enough threads for one
 * thread per element.
 */
class SyclThreadDirectExecPolicy : public SyclExecPolicy
{
public:
  SyclThreadDirectExecPolicy(const size_t blockDim)
    : blockDim_(blockDim)
  {}

  SyclThreadDirectExecPolicy(const SyclThreadDirectExecPolicy& ex)
    : blockDim_(ex.blockDim_)
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

  virtual SyclExecPolicy* clone() const
  {
    return static_cast<SyclExecPolicy*>(new SyclThreadDirectExecPolicy(*this));
  }

private:
  const size_t blockDim_;
};

/*
 * A kernel execution policy for kernels that use grid stride loops.
 * The number of threads per block (blockSize) can be set to anything.
 * The number of blocks (gridSize) can be set to anything.
 */
class SyclGridStrideExecPolicy : public SyclExecPolicy
{
public:
  SyclGridStrideExecPolicy(const size_t blockDim, const size_t gridDim)
    : blockDim_(blockDim), gridDim_(gridDim)
  {}

  SyclGridStrideExecPolicy(const SyclGridStrideExecPolicy& ex)
    : blockDim_(ex.blockDim_), gridDim_(ex.gridDim_)
  {}

  virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const
  {
    return gridDim_;
  }

  virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const
  {
    return blockDim_;
  }

  virtual SyclExecPolicy* clone() const
  {
    return static_cast<SyclExecPolicy*>(new SyclGridStrideExecPolicy(*this));
  }

private:
  const size_t blockDim_;
  const size_t gridDim_;
};


/*
 * A kernel execution policy for performing a reduction across indvidual thread
 * blocks. The number of threads per block (blockSize) can be set to anything.
 * The number of blocks (gridSize) can be set to any value greater than or equal
 * to 0. If it is set to 0, then the grid size will be chosen so that there are
 * at most two work units per thread.
 */
class SyclBlockReduceExecPolicy : public SyclExecPolicy
{
public:
  SyclBlockReduceExecPolicy(const size_t blockDim, const size_t gridDim = 0)
    : blockDim_(blockDim), gridDim_(gridDim)
  {}

  SyclBlockReduceExecPolicy(const SyclBlockReduceExecPolicy& ex)
    : blockDim_(ex.blockDim_), gridDim_(ex.gridDim_)
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

  virtual SyclExecPolicy* clone() const
  {
    return static_cast<SyclExecPolicy*>(new SyclBlockReduceExecPolicy(*this));
  }

private:
  const size_t blockDim_;
  const size_t gridDim_;
};

} // namespace sundials

typedef sundials::SyclExecPolicy SUNSyclExecPolicy;
typedef sundials::SyclThreadDirectExecPolicy SUNSyclThreadDirectExecPolicy;
typedef sundials::SyclGridStrideExecPolicy SUNSyclGridStrideExecPolicy;
typedef sundials::SyclBlockReduceExecPolicy SUNSyclBlockReduceExecPolicy;

#endif
