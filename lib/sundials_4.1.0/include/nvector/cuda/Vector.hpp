/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, and Cody J. Balos @ LLNL
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


/*
 * Vector class
 *
 * Manages vector data layout for CUDA implementation of N_Vector.
 *
 */

#ifndef _NVECTOR_CUDA_HPP_
#define _NVECTOR_CUDA_HPP_

#include <cstdlib>
#include <iostream>

#include <cuda_runtime.h>
#include "ThreadPartitioning.hpp"

#include <sundials/sundials_config.h>
#include <sundials/sundials_mpi.h>

#if SUNDIALS_MPI_ENABLED
#include <nvector/nvector_mpicuda.h>
#else
#include <nvector/nvector_cuda.h>
#endif

namespace suncudavec
{

template <typename T, typename I>
class Vector : public _N_VectorContent_Cuda
{
public:
  Vector(I N, bool use_managed_memory = false, bool allocate_data = true, T* const h_vec = nullptr, T* const d_vec = nullptr)
  : size_(N),
    mem_size_(N*sizeof(T)),
    global_size_(N),
    ownPartitioning_(true),
    ownData_(allocate_data),
    managed_mem_(use_managed_memory),
    h_vec_(h_vec),
    d_vec_(d_vec),
    comm_(0)
  {
    // Set partitioning
    partStream_ = new StreamPartitioning<T, I>(N, 256);
    partReduce_ = new ReducePartitioning<T, I>(N, 256);

    // Allocate data arrays
    if (allocate_data)
      allocate();
  }

  Vector(I N, cudaStream_t stream,
         bool use_managed_memory = false, bool allocate_data = true, T* const h_vec = nullptr, T* const d_vec = nullptr)
  : size_(N),
    mem_size_(N*sizeof(T)),
    global_size_(N),
    ownPartitioning_(true),
    ownData_(allocate_data),
    managed_mem_(use_managed_memory),
    h_vec_(h_vec),
    d_vec_(d_vec),
    comm_(0)
  {
    // Set partitioning
    partStream_ = new StreamPartitioning<T, I>(N, 256, stream);
    partReduce_ = new ReducePartitioning<T, I>(N, 256, stream);

    // Allocate data arrays
    if (allocate_data)
      allocate();
  }
  
  Vector(SUNMPI_Comm comm, I N, I Nglobal,
         bool use_managed_memory = false, bool allocate_data = true, T* const h_vec = nullptr, T* const d_vec = nullptr)
  : size_(N),
    mem_size_(N*sizeof(T)),
    global_size_(Nglobal),
    ownPartitioning_(true),
    ownData_(allocate_data),
    managed_mem_(use_managed_memory),
    h_vec_(h_vec),
    d_vec_(d_vec),
    comm_(comm)
  {
    // Set partitioning
    partStream_ = new StreamPartitioning<T, I>(N, 256);
    partReduce_ = new ReducePartitioning<T, I>(N, 256);

    // Allocate data arrays
    if (allocate_data)
      allocate();
  }

  Vector(SUNMPI_Comm comm, I N, I Nglobal, cudaStream_t stream,
         bool use_managed_memory = false, bool allocate_data = true, T* const h_vec = nullptr, T* const d_vec = nullptr) 
  : size_(N),
    mem_size_(N*sizeof(T)),
    global_size_(Nglobal),
    ownPartitioning_(true),
    ownData_(allocate_data),
    managed_mem_(use_managed_memory),
    h_vec_(h_vec),
    d_vec_(d_vec),
    comm_(comm)
  {
    // Set partitioning
    partStream_ = new StreamPartitioning<T, I>(N, 256, stream);
    partReduce_ = new ReducePartitioning<T, I>(N, 256, stream);

    // Allocate data arrays
    if (allocate_data)
      allocate();
  }
  
  // Copy constructor does not copy data array values
  explicit Vector(const Vector& v)
  : size_(v.size()),
    mem_size_(size_*sizeof(T)),
    global_size_(v.global_size_),
    partStream_(v.partStream_),
    partReduce_(v.partReduce_),
    ownPartitioning_(false),
    ownData_(true),
    managed_mem_(v.managed_mem_),
    h_vec_(nullptr),
    d_vec_(nullptr),
    comm_(v.comm_)
  {
    allocate();
  }

  ~Vector()
  {
    cudaError_t err;
    
    if (ownPartitioning_) {
      delete partReduce_;
      delete partStream_;
    }
    
    if (ownData_) {
      if (!managed_mem_)
        free(h_vec_);
      
      err = cudaFree(d_vec_);
      if(err != cudaSuccess)
        std::cerr << "Failed to free device vector (error code " << err << ")!\n";
    
      d_vec_ = nullptr;
      h_vec_ = nullptr;
    }
  }

  void allocate()
  {
    if (managed_mem_) {
      allocateManaged();
    } else {
      allocateUnmanaged();
    }
  }

  void allocateManaged()
  {
    cudaError_t err;
    err = cudaMallocManaged((void**) &d_vec_, mem_size_);
    if (err != cudaSuccess)
      std::cerr << "Failed to allocate managed vector (error code " << err << ")!\n";
    h_vec_ = d_vec_;
  }

  void allocateUnmanaged()
  {
    cudaError_t err;
    
    h_vec_ = static_cast<T*>(malloc(mem_size_));
    if(h_vec_ == nullptr)
      std::cerr << "Failed to allocate host vector!\n";
    
    err = cudaMalloc((void**) &d_vec_, mem_size_);
    if(err != cudaSuccess)
      std::cerr << "Failed to allocate device vector (error code " << err << ")!\n";
  }
  
  int size() const
  {
    return size_;
  }

  int sizeGlobal() const
  {
    return global_size_;
  }

  SUNMPI_Comm comm() const
  {
    return comm_;
  }

  T* host()
  {
    return h_vec_;
  }

  const T* host() const
  {
    return h_vec_;
  }

  T* device()
  {
    return d_vec_;
  }

  const T* device() const
  {
    return d_vec_;
  }

  bool isManaged() const
  {
    return managed_mem_;
  }

  void copyToDev()
  {
    cudaError_t err = cudaMemcpy(d_vec_, h_vec_, mem_size_, cudaMemcpyHostToDevice);
    if(err != cudaSuccess)
      std::cerr << "Failed to copy vector from host to device (error code " << err << ")!\n";
  }

  void copyFromDev()
  {
    cudaError_t err = cudaMemcpy(h_vec_, d_vec_, mem_size_, cudaMemcpyDeviceToHost);
    if(err != cudaSuccess)
      std::cerr << "Failed to copy vector from device to host (error code " << err << ")!\n";
  }

  ThreadPartitioning<T, I>& partStream() const
  {
    return *partStream_;
  }

  ThreadPartitioning<T, I>& partReduce() const
  {
    return *partReduce_;
  }


private:
  I size_;
  I mem_size_;
  I global_size_;
  T* h_vec_;
  T* d_vec_;
  ThreadPartitioning<T, I>* partStream_;
  ThreadPartitioning<T, I>* partReduce_;
  bool ownPartitioning_;
  bool ownData_;
  bool managed_mem_;
  SUNMPI_Comm comm_;
  
};


} // namespace suncudavec




#endif // _NVECTOR_CUDA_HPP_
