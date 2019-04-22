/*
 * -----------------------------------------------------------------
 * Programmer(s): David Gardner @ LLNL
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


#ifndef _VECTOR_ARRAY_KERNELS_CUH_
#define _VECTOR_ARRAY_KERNELS_CUH_

#include <limits>
#include <cuda_runtime.h>


namespace suncudavec
{


/* -----------------------------------------------------------------
 * The namespace for CUDA kernels
 *
 * Reduction CUDA kernels in nvector are based in part on "reduction"
 * example in NVIDIA Corporation CUDA Samples, and parallel reduction
 * examples in textbook by J. Cheng at al. "CUDA C Programming".
 * -----------------------------------------------------------------
 */
namespace math_kernels
{


/*
 * -----------------------------------------------------------------------------
 * fused vector operation kernels
 * -----------------------------------------------------------------------------
 */

/*
 * Computes the linear combination of nv vectors
 */
template <typename T, typename I>
__global__ void
linearCombinationKernel(int nv, T* c, T** xd, T* zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n) {
    zd[i] = c[0]*xd[0][i];
    for (int j=1; j<nv; j++)
      zd[i] += c[j]*xd[j][i];
  }
}

/*
 * Computes the scaled sum of one vector with nv other vectors
 */
template <typename T, typename I>
__global__ void
scaleAddMultiKernel(int nv, T* c, T* xd, T** yd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int j=0; j<nv; j++)
      zd[j][i] = c[j] * xd[i] + yd[j][i];
}


/*
 * Dot product of one vector with nv other vectors.
 *
 */
template <typename T, typename I>
__global__ void
dotProdMultiKernel(int nv, T* xd, T** yd, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // Initialize shared memory to zero
  for (int k=0; k<nv; k++)
    shmem[tid + k*blockDim.x] = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] = xd[i] * yd[k][i];
  if (i + blockDim.x < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] += (xd[i + blockDim.x] * yd[k][i + blockDim.x]);

  __syncthreads();

  // Perform blockwise reduction in shared memory
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j)
      for (int k=0; k<nv; k++)
        shmem[tid + k*blockDim.x] += shmem[tid + j + k*blockDim.x];

    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    for (int k=0; k<nv; k++)
      out[blockIdx.x + k*gridDim.x] = shmem[k*blockDim.x];
}


/*
 * Sums all elements of the vector.
 *
 */
template <typename T, typename I>
__global__ void
sumReduceVectorKernel(int nv, T* x, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // First reduction step before storing data in shared memory.
  if (i < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] = x[i];
  if (i + blockDim.x < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] += x[i+blockDim.x];

  __syncthreads();

  // Perform reduction block-wise in shared memory.
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j)
      for (int k=0; k<nv; k++)
        shmem[tid + k*blockDim.x] += shmem[tid + j + k*blockDim.x];

    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    for (int k=0; k<nv; k++)
      out[blockIdx.x + k*gridDim.x] = shmem[k*blockDim.x];
}



/*
 * -----------------------------------------------------------------------------
 * vector array operation kernels
 * -----------------------------------------------------------------------------
 */

/*
 * Computes the linear sum of multiple vectors
 */
template <typename T, typename I>
__global__ void
linearSumVectorArrayKernel(int nv, T a, T** xd, T b, T** yd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int j=0; j<nv; j++)
      zd[j][i] = a * xd[j][i] + b * yd[j][i];
}


/*
 * Scales multiple vectors
 */
template <typename T, typename I>
__global__ void
scaleVectorArrayKernel(int nv, T* c, T** xd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int j=0; j<nv; j++)
      zd[j][i] = c[j] * xd[j][i];
}


/*
 * Sets multiple vectors equal to a constant
 */
template <typename T, typename I>
__global__ void
constVectorArrayKernel(int nv, T c, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int j=0; j<nv; j++)
      zd[j][i] = c;
}


/*
 * WRMS norm of nv vectors.
 *
 */
template <typename T, typename I>
__global__ void
wL2NormSquareVectorArrayKernel(int nv, T** xd, T** wd, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // Initialize shared memory to zero
  for (int k=0; k<nv; k++)
    shmem[tid + k*blockDim.x] = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] = xd[k][i] * wd[k][i] * xd[k][i] * wd[k][i];
  if (i + blockDim.x < n)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] += (xd[k][i + blockDim.x] * wd[k][i + blockDim.x]
                                    * xd[k][i + blockDim.x] * wd[k][i + blockDim.x]);

  __syncthreads();

  // Perform blockwise reduction in shared memory
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j)
      for (int k=0; k<nv; k++)
        shmem[tid + k*blockDim.x] += shmem[tid + j + k*blockDim.x];

    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    for (int k=0; k<nv; k++)
      out[blockIdx.x + k*gridDim.x] = shmem[k*blockDim.x];
}


/*
 * Masked WRMS norm of nv vectors.
 *
 */
template <typename T, typename I>
__global__ void
wL2NormSquareMaskVectorArrayKernel(int nv, T** xd, T** wd, T* id, T* out, I n)
{
  extern __shared__ T shmem[];

  I tid = threadIdx.x;
  I i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  // Initialize shared memory to zero
  for (int k=0; k<nv; k++)
    shmem[tid + k*blockDim.x] = 0.0;

  // First reduction step before storing data in shared memory.
  if (i < n && id[i] > 0.0)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] = xd[k][i] * wd[k][i] * xd[k][i] * wd[k][i];
  if (i + blockDim.x < n && id[i + blockDim.x] > 0.0)
    for (int k=0; k<nv; k++)
      shmem[tid + k*blockDim.x] += (xd[k][i + blockDim.x] * wd[k][i + blockDim.x]
                                    * xd[k][i + blockDim.x] * wd[k][i + blockDim.x]);

  __syncthreads();

  // Perform blockwise reduction in shared memory
  for (I j = blockDim.x/2; j > 0; j >>= 1) {
    if (tid < j)
      for (int k=0; k<nv; k++)
        shmem[tid + k*blockDim.x] += shmem[tid + j + k*blockDim.x];

    __syncthreads();
  }

  // Copy reduction result for each block to global memory
  if (tid == 0)
    for (int k=0; k<nv; k++)
      out[blockIdx.x + k*gridDim.x] = shmem[k*blockDim.x];
}


/*
 * Computes the scaled sum of a vector array with multiple other vector arrays
 */
template <typename T, typename I>
__global__ void
scaleAddMultiVectorArrayKernel(int nv, int ns, T* c, T** xd, T** yd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n)
    for (int k=0; k<nv; k++)
      for (int j=0; j<ns; j++)
        zd[k*ns+j][i] = c[j] * xd[k][i] + yd[k*ns+j][i];
}


/*
 * Computes the scaled sum of a vector array with multiple other vector arrays
 */
template <typename T, typename I>
__global__ void
linearCombinationVectorArrayKernel(int nv, int ns, T* c, T** xd, T** zd, I n)
{
  I i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n) {
    for (int k=0; k<nv; k++) {
      zd[k][i] = c[0]*xd[k*ns][i];
      for (int j=1; j<ns; j++) {
        zd[k][i] += c[j]*xd[k*ns+j][i];
      }
    }
  }
}

} // namespace math_kernels






/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */

template <typename T, typename I>
inline cudaError_t linearCombination(int nvec, T* c, Vector<T,I>** X, Vector<T,I>* Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = X[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();
  const cudaStream_t stream   = p.stream();

  math_kernels::linearCombinationKernel<<<grid, block, 0, stream>>>(
      nvec,
      d_c,
      d_Xd,
      Z->device(),
      Z->size()
  );

  // Free host array
  delete[] h_Xd;

  // Free device arrays
  err = cudaFree(d_c);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t scaleAddMulti(int nvec, T* c, Vector<T,I>* X,
                                 Vector<T,I>** Y, Vector<T,I>** Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Yd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = Y[i]->device();

  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();
  const cudaStream_t stream   = p.stream();

  math_kernels::scaleAddMultiKernel<<<grid, block, 0, stream>>>(
      nvec,
      d_c,
      X->device(),
      d_Yd,
      d_Zd,
      X->size()
  );

  // Free host array
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_c);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t dotProdMulti(int nvec, Vector<T,I>* x, Vector<T,I>** Y,
                                T* dots)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Yd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = Y[i]->device();

  // Copy array of device pointers to device from host
  T** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = x->partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = nvec*block*sizeof(T);
  const cudaStream_t stream   = p.stream();

  // Allocate reduction buffer on device
  T* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();

  math_kernels::dotProdMultiKernel<T,I><<<grid, block, shMemSize, stream>>>(
      nvec,
      x->device(),
      d_Yd,
      d_buff,
      x->size()
  );

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax) {

    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    math_kernels::sumReduceVectorKernel<T,I><<<grid, block, shMemSize, stream>>>(
        nvec,
        d_buff,
        d_buff,
        n
    );

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  T* h_buff = new T[nvec*n*sizeof(T)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(T), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) return cudaGetLastError();

  for (int k=0; k<nvec; k++) {
    dots[k] = h_buff[k*n];
    for (int i=1; i<n; i++){
      dots[k] += h_buff[i + k*n];
    }
  }

  // Free host array
  delete[] h_Yd;
  delete[] h_buff;

  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_buff);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

template <typename T, typename I>
inline cudaError_t linearSumVectorArray(int nvec, T a, Vector<T,I>** X, T b,
                                        Vector<T,I>** Y, Vector<T,I>** Z)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Yd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = Y[i]->device();

  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();
  const cudaStream_t stream   = p.stream();

  math_kernels::linearSumVectorArrayKernel<<<grid, block, 0, stream>>>(
      nvec,
      a,
      d_Xd,
      b,
      d_Yd,
      d_Zd,
      Z[0]->size()
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t scaleVectorArray(int nvec, T* c, Vector<T,I>** X,
                                    Vector<T,I>** Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();
  const cudaStream_t stream   = p.stream();

  math_kernels::scaleVectorArrayKernel<<<grid, block, 0, stream>>>(
      nvec,
      d_c,
      d_Xd,
      d_Zd,
      Z[0]->size()
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t constVectorArray(int nvec, T c, Vector<T,I>** Z)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();
  const cudaStream_t stream   = p.stream();

  math_kernels::constVectorArrayKernel<<<grid, block, 0, stream>>>(
      nvec,
      c,
      d_Zd,
      Z[0]->size()
  );

  // Free host array
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t wL2NormSquareVectorArray(int nvec, Vector<T,I>** X,
                                            Vector<T,I>** W, T* nrm)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Wd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = W[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Wd;
  err = cudaMalloc((void**) &d_Wd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Wd, h_Wd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = X[0]->partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = nvec*block*sizeof(T);
  const cudaStream_t stream   = p.stream();

  // Allocate reduction buffer on device
  T* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();

  math_kernels::wL2NormSquareVectorArrayKernel<<<grid, block, shMemSize, stream>>>(
      nvec,
      d_Xd,
      d_Wd,
      d_buff,
      X[0]->size()
  );

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax) {

    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    math_kernels::sumReduceVectorKernel<T,I><<<grid, block, shMemSize, stream>>>(
        nvec,
        d_buff,
        d_buff,
        n
    );

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  T* h_buff = new T[nvec*n*sizeof(T)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(T), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) return cudaGetLastError();

  for (int k=0; k<nvec; k++) {
    nrm[k] = h_buff[k*n];
    for (int i=1; i<n; i++){
      nrm[k] += h_buff[i + k*n];
    }
  }

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;
  delete[] h_buff;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Wd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_buff);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t wL2NormSquareMaskVectorArray(int nvec, Vector<T,I>** X,
                                           Vector<T,I>** W, Vector<T,I>* ID,
                                           T* nrm)
{
  cudaError_t err;

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Wd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = W[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Wd;
  err = cudaMalloc((void**) &d_Wd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Wd, h_Wd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = X[0]->partReduce();
  unsigned grid               = p.grid();
  unsigned block              = p.block();
  unsigned shMemSize          = nvec*block*sizeof(T);
  const cudaStream_t stream   = p.stream();

  // Allocate reduction buffer on device
  T* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();

  math_kernels::wL2NormSquareMaskVectorArrayKernel<<<grid, block, shMemSize, stream>>>(
      nvec,
      d_Xd,
      d_Wd,
      ID->device(),
      d_buff,
      X[0]->size()
  );

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax) {

    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    math_kernels::sumReduceVectorKernel<T,I><<<grid, block, shMemSize, stream>>>(
        nvec,
        d_buff,
        d_buff,
        n
    );

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  T* h_buff = new T[nvec*n*sizeof(T)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(T), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) return cudaGetLastError();

  for (int k=0; k<nvec; k++) {
    nrm[k] = h_buff[k*n];
    for (int i=1; i<n; i++){
      nrm[k] += h_buff[i + k*n];
    }
  }

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;
  delete[] h_buff;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Wd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_buff);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t scaleAddMultiVectorArray(int nvec, int nsum, T* c,
                                            Vector<T,I>** X, Vector<T,I>** Y,
                                            Vector<T,I>** Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nsum*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Xd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Yd = new T*[nsum*nvec];
  for (int i=0; i<nsum*nvec; i++)
    h_Yd[i] = Y[i]->device();

  T** h_Zd = new T*[nsum*nvec];
  for (int i=0; i<nsum*nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nsum*nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nsum*nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nsum*nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nsum*nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();
  const cudaStream_t stream   = p.stream();

  math_kernels::scaleAddMultiVectorArrayKernel<<<grid, block, 0, stream>>>(
      nvec,
      nsum,
      d_c,
      d_Xd,
      d_Yd,
      d_Zd,
      Z[0]->size()
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


template <typename T, typename I>
inline cudaError_t linearCombinationVectorArray(int nvec, int nsum, T* c,
                                                Vector<T,I>** X, Vector<T,I>** Z)
{
  cudaError_t err;

  // Copy c array to device
  T* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(T));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nsum*sizeof(T), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  T** h_Xd = new T*[nsum*nvec];
  for (int i=0; i<nsum*nvec; i++)
    h_Xd[i] = X[i]->device();

  T** h_Zd = new T*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = Z[i]->device();

  // Copy array of device pointers to device from host
  T** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nsum*nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nsum*nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  T** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(T*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(T*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Set partitioning
  ThreadPartitioning<T, I>& p = Z[0]->partStream();
  const I grid                = p.grid();
  const unsigned block        = p.block();
  const cudaStream_t stream   = p.stream();

  math_kernels::linearCombinationVectorArrayKernel<<<grid, block, 0, stream>>>(
      nvec,
      nsum,
      d_c,
      d_Xd,
      d_Zd,
      Z[0]->size()
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return cudaGetLastError();
}


} // namespace nvec



#endif // _VECTOR_ARRAY_KERNELS_CUH_
