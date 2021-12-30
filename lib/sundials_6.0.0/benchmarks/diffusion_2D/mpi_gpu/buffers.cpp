/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * CUDA/HIP functions for exchange buffers
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

#if defined(USE_HIP)
#include <hip/hip_runtime.h>
#define BLOCK_SIZE   256
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#elif defined(USE_CUDA)
#include <cuda_runtime.h>
#define BLOCK_SIZE   256
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#else
#error Define USE_CUDA or USE_HIP
#endif

// Pack exchange buffers kernel
__global__ void pack_buffers_kernel(const sunindextype nx_loc,
                                    const sunindextype ny_loc,
                                    const realtype *u,
                                    realtype *wbuf, realtype *ebuf,
                                    realtype *sbuf, realtype *nbuf)
{
  // Thread ID
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  // West and East faces
  if (i < ny_loc)
  {
    if (wbuf) wbuf[i] = u[i * nx_loc];
    if (ebuf) ebuf[i] = u[(i + 1) * nx_loc - 1];
  }

  // South and North faces
  if (i < nx_loc)
  {
    if (sbuf) sbuf[i] = u[i];
    if (nbuf) nbuf[i] = u[i + (ny_loc - 1) * nx_loc];
  }
}


// Pack exchange buffers
int UserData::pack_buffers(const N_Vector u)
{
  // Access data array
  const realtype *uarray = N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(u));
  if (check_flag((void *) uarray, "N_VGetDeviceArrayPointer", 0)) return -1;

  sunindextype maxdim = max(nx_loc, ny_loc);
  dim3 block(BLOCK_SIZE);
  dim3 grid(ICEIL(maxdim, BLOCK_SIZE));

  pack_buffers_kernel<<<grid,block>>>(nx_loc, ny_loc, uarray,
                                      Wsend, Esend, Ssend, Nsend);

  return 0;
}


// Allocate exchange buffers
int UserData::allocate_buffers()
{
  if (HaveNbrW)
  {
#if defined(USE_HIP)
    hipMalloc(&Wrecv, ny_loc * sizeof(realtype));
    hipMalloc(&Wsend, ny_loc * sizeof(realtype));
#elif defined(USE_CUDA)
    cudaMalloc(&Wrecv, ny_loc * sizeof(realtype));
    cudaMalloc(&Wsend, ny_loc * sizeof(realtype));
#else
#error Define USE_CUDA or USE_HIP
#endif
  }

  if (HaveNbrE)
  {
#if defined(USE_HIP)
    hipMalloc(&Erecv, ny_loc * sizeof(realtype));
    hipMalloc(&Esend, ny_loc * sizeof(realtype));
#elif defined(USE_CUDA)
    cudaMalloc(&Erecv, ny_loc * sizeof(realtype));
    cudaMalloc(&Esend, ny_loc * sizeof(realtype));
#else
#error Define USE_CUDA or USE_HIP
#endif
  }

  if (HaveNbrS)
  {
#if defined(USE_HIP)
    hipMalloc(&Srecv, nx_loc * sizeof(realtype));
    hipMalloc(&Ssend, nx_loc * sizeof(realtype));
#elif defined(USE_CUDA)
    cudaMalloc(&Srecv, nx_loc * sizeof(realtype));
    cudaMalloc(&Ssend, nx_loc * sizeof(realtype));
#else
#error Define USE_CUDA or USE_HIP
#endif
  }

  if (HaveNbrN)
  {
#if defined(USE_HIP)
    hipMalloc(&Nrecv, nx_loc * sizeof(realtype));
    hipMalloc(&Nsend, nx_loc * sizeof(realtype));
#elif defined(USE_CUDA)
    cudaMalloc(&Nrecv, nx_loc * sizeof(realtype));
    cudaMalloc(&Nsend, nx_loc * sizeof(realtype));
#else
#error Define USE_CUDA or USE_HIP
#endif
  }

  return 0;
}


// Free exchange buffers
int UserData::free_buffers()
{
#if defined(USE_HIP)
  if (Wrecv != NULL)  hipFree(Wrecv);
  if (Wsend != NULL)  hipFree(Wsend);
  if (Erecv != NULL)  hipFree(Erecv);
  if (Esend != NULL)  hipFree(Esend);
  if (Srecv != NULL)  hipFree(Srecv);
  if (Ssend != NULL)  hipFree(Ssend);
  if (Nrecv != NULL)  hipFree(Nrecv);
  if (Nsend != NULL)  hipFree(Nsend);
#elif defined(USE_CUDA)
  if (Wrecv != NULL)  cudaFree(Wrecv);
  if (Wsend != NULL)  cudaFree(Wsend);
  if (Erecv != NULL)  cudaFree(Erecv);
  if (Esend != NULL)  cudaFree(Esend);
  if (Srecv != NULL)  cudaFree(Srecv);
  if (Ssend != NULL)  cudaFree(Ssend);
  if (Nrecv != NULL)  cudaFree(Nrecv);
  if (Nsend != NULL)  cudaFree(Nsend);
#else
#error Define USE_CUDA or USE_HIP
#endif

  return 0;
}
