/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 */


/*
 * Vector class
 *
 * Manages vector data layout for RAJA implementation of N_Vector.
 *
 */

#ifndef _NVECTOR_RAJA_HPP_
#define _NVECTOR_RAJA_HPP_

#include <cstdlib>
#include <iostream>

#include <sundials/sundials_config.h>
#include <sundials/sundials_mpi.h>

#if SUNDIALS_MPI_ENABLED
#include <nvector/nvector_mpiraja.h>
#else
#include <nvector/nvector_raja.h>
#endif

namespace sunrajavec
{

template <typename T, typename I>
class Vector : public _N_VectorContent_Raja
{
public:
  Vector(I N)
  : size_(N),
    mem_size_(N*sizeof(T)),
    global_size_(N),
    comm_(0)
  {
    allocate();
  }

  Vector(SUNMPI_Comm comm, I N, I Nglobal)
  : size_(N),
    mem_size_(N*sizeof(T)),
    global_size_(Nglobal),
    comm_(comm)
  {
    allocate();
  }

  // Copy constructor does not copy values
  explicit Vector(const Vector& v)
  : size_(v.size()),
    mem_size_(size_*sizeof(T)),
    global_size_(v.global_size_),
    comm_(v.comm_)
  {
    allocate();
  }

  ~Vector()
  {
    cudaError_t err;
    free(h_vec_);
    err = cudaFree(d_vec_);
    if(err != cudaSuccess)
      std::cout << "Failed to free device vector (error code " << err << ")!\n";
  }


  void allocate()
  {
    cudaError_t err;
    h_vec_ = static_cast<T*>(malloc(mem_size_));
    if(h_vec_ == NULL)
      std::cout << "Failed to allocate host vector!\n";
    err = cudaMalloc((void**) &d_vec_, mem_size_);
    if(err != cudaSuccess)
      std::cout << "Failed to allocate device vector (error code " << err << ")!\n";
  }

  int size() const
  {
    return size_;
  }

  int sizeGlobal() const
  {
    return global_size_;
  }

  SUNMPI_Comm comm()
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

private:
  I size_;
  I mem_size_;
  I global_size_;
  T* h_vec_;
  T* d_vec_;
  SUNMPI_Comm comm_;
};


} // namespace sunrajavec



#endif // _NVECTOR_RAJA_HPP_
