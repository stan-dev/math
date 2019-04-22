/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 */



#ifndef _THREAD_PARTITIONING_HPP_
#define _THREAD_PARTITIONING_HPP_

#include <iostream>
#include <cuda_runtime.h>

namespace suncudavec
{

template<class T, class I>
class ThreadPartitioning
{
public:
  ThreadPartitioning()
  : block_(1),
    grid_(1),
    shMemSize_(0),
    stream_(0),
    bufferSize_(0),
    d_buffer_(nullptr),
    h_buffer_(nullptr)
  {}

  ThreadPartitioning(unsigned block)
  : block_(block),
    grid_(1),
    shMemSize_(0),
    stream_(0),
    bufferSize_(0),
    d_buffer_(nullptr),
    h_buffer_(nullptr)
  {}

  explicit ThreadPartitioning(ThreadPartitioning<T, I>& p)
  : block_(p.block_),
    grid_(p.grid_),
    shMemSize_(p.shMemSize_),
    stream_(p.stream_)
  {}

  virtual ~ThreadPartitioning(){}

  unsigned grid() const
  {
    return grid_;
  }

  unsigned block() const
  {
    return block_;
  }

  unsigned shmem() const
  {
    return shMemSize_;
  }

  cudaStream_t stream() const
  {
    return stream_;
  }

  unsigned int buffSize()
  {
    return bufferSize_;
  }

  T* devBuffer()
  {
    return d_buffer_;
  }

  const T* devBuffer() const
  {
    return d_buffer_;
  }

  T* hostBuffer()
  {
    return h_buffer_;
  }

  const T* hostBuffer() const
  {
    return h_buffer_;
  }

  void setStream(const cudaStream_t& stream)
  {
    stream_ = stream;
  }

  virtual void copyFromDevBuffer(unsigned int n) const
  {
    std::cerr << "Trying to copy buffer from base class!\n";
  }

  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, unsigned& shMemSize,
                              cudaStream_t& stream)
  {
    block = 1;
    grid  = 1;
    shMemSize = 0;
    stream = 0;
    std::cerr << "Trying to set partitioning from base class!\n";

    return 0;
  }

  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, unsigned& shMemSize)
  {
    block = 1;
    grid  = 1;
    shMemSize = 0;
    std::cerr << "Trying to set partitioning from base class!\n";

    return 0;
  }
  
  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, cudaStream_t& stream)
  {
    block = 1;
    grid  = 1;
    stream = 0;
    std::cerr << "Trying to set partitioning from base class!\n";

    return 0;
  }
  
  virtual int setPartitioning(I N, unsigned& grid, unsigned& block)
  {
    block = 1;
    grid  = 1;
    std::cerr << "Trying to set partitioning from base class!\n";

    return 0;
  }


protected:
  unsigned block_;
  unsigned grid_;
  unsigned shMemSize_;
  unsigned bufferSize_;
  cudaStream_t stream_;
  T* d_buffer_;
  T* h_buffer_;

}; // class ThreadPartitioning



template<class T, class I>
class StreamPartitioning : public ThreadPartitioning<T, I>
{
  using ThreadPartitioning<T, I>::block_;
  using ThreadPartitioning<T, I>::grid_;
  using ThreadPartitioning<T, I>::stream_;

public:
  StreamPartitioning(I N, unsigned block, cudaStream_t stream)
  : ThreadPartitioning<T, I>(block)
  {
    grid_ = (N + block_ - 1) / block_;
    stream_ = stream;
  }
  
  StreamPartitioning(I N, unsigned block)
  : ThreadPartitioning<T, I>(block)
  {
    grid_ = (N + block_ - 1) / block_;
  }

  explicit StreamPartitioning(StreamPartitioning<T, I>& p)
  : ThreadPartitioning<T, I>(p)
  {
  }

  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, unsigned& shMemSize,
                              cudaStream_t& stream)
  {
    block = block_;
    grid  = (N + block_ - 1) / block_;
    shMemSize = 0;
    stream = stream_;

    return 0;
  }
  
  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, unsigned& shMemSize)
  {
    block = block_;
    grid  = (N + block_ - 1) / block_;
    shMemSize = 0;

    return 0;
  }

  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, cudaStream_t& stream)
  {
    block = block_;
    grid  = (N + block_ - 1) / block_;
    stream = stream_;

    return 0;
  }
  
  virtual int setPartitioning(I N, unsigned& grid, unsigned& block)
  {
    block = block_;
    grid  = (N + block_ - 1) / block_;

    return 0;
  }

}; // class StreamPartitioning


template<class T, class I=int>
class ReducePartitioning : public ThreadPartitioning<T, I>
{
  using ThreadPartitioning<T, I>::block_;
  using ThreadPartitioning<T, I>::grid_;
  using ThreadPartitioning<T, I>::shMemSize_;
  using ThreadPartitioning<T, I>::stream_;
  using ThreadPartitioning<T, I>::bufferSize_;
  using ThreadPartitioning<T, I>::d_buffer_;
  using ThreadPartitioning<T, I>::h_buffer_;

public:
  ReducePartitioning(I N, unsigned block, cudaStream_t stream)
  : ThreadPartitioning<T, I>(block)
  {
    grid_ = (N + (block_ * 2 - 1)) / (block_ * 2);
    shMemSize_ = block_*sizeof(T);
    stream_ = stream;
    allocateBuffer();
  }
  
  ReducePartitioning(I N, unsigned block)
  : ThreadPartitioning<T, I>(block)
  {
    grid_ = (N + (block_ * 2 - 1)) / (block_ * 2);
    shMemSize_ = block_*sizeof(T);
    allocateBuffer();
  }

  explicit ReducePartitioning(ReducePartitioning<T, I>& p)
  : ThreadPartitioning<T, I>(p)
  {
    shMemSize_ = p.shMemSize_;
    allocateBuffer();
  }

  ~ReducePartitioning()
  {
    cudaError_t err;
    if (bufferSize_ > 0)
      free(h_buffer_);
    if (bufferSize_ > 0)
    {
      err = cudaFree(d_buffer_);
      if(err != cudaSuccess)
        std::cerr << "Failed to free device vector (error code " << err << ")!\n";
    }
  }

  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, unsigned& shMemSize,
                              cudaStream_t& stream)
  {
    block = block_;
    grid  = (N + (block_ * 2 - 1)) / (block_ * 2);
    shMemSize = block_ * sizeof(T);
    stream = stream_;

    return 0;
  }
  
  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, unsigned& shMemSize)
  {
    block = block_;
    grid  = (N + (block_ * 2 - 1)) / (block_ * 2);
    shMemSize = block_ * sizeof(T);

    return 0;
  }

  virtual int setPartitioning(I N, unsigned& grid, unsigned& block, cudaStream_t& stream)
  {
    block = block_;
    grid  = (N + (block_ * 2 - 1)) / (block_ * 2);
    stream = stream_;

    return 0;
  }
  
  virtual int setPartitioning(I N, unsigned& grid, unsigned& block)
  {
    block = block_;
    grid  = (N + (block_ * 2 - 1)) / (block_ * 2);

    return 0;
  }

  virtual void copyFromDevBuffer(unsigned int n) const
  {
    cudaError_t err = cudaMemcpy(h_buffer_, d_buffer_, n*sizeof(T), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess)
      std::cerr << "Failed to copy vector from device to host (error code " << err << ")!\n";
  }

private:
  int allocateBuffer()
  {
    bufferSize_ = grid_ * sizeof(T);
    h_buffer_ = static_cast<T*>(malloc(bufferSize_));
    if(h_buffer_ == NULL)
      std::cerr << "Failed to allocate host vector!\n";

    cudaError_t err;
    err = cudaMalloc((void**) &d_buffer_, bufferSize_);
    if(err != cudaSuccess)
      std::cerr << "Failed to allocate device vector (error code " << err << ")!\n";

    return 0;
  }

}; // class ReducePartitioning


} // namespace suncudavec

#endif // _THREAD_PARTITIONING_HPP_
