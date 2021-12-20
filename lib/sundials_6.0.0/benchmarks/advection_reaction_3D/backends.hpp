/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * ---------------------------------------------------------------------------*/

#ifndef ADVECTION_REACTION_3D_BACKENDS_HPP
#define ADVECTION_REACTION_3D_BACKENDS_HPP

#include <RAJA/RAJA.hpp>

#include "nvector/nvector_mpiplusx.h"               /* MPI+X N_Vector            */
#if defined(USE_RAJA_NVEC)
#include "nvector/nvector_raja.h"                   /* RAJA N_Vector             */
#elif defined(USE_OMPDEV_NVEC)
#include <omp.h>
#include "nvector/nvector_openmpdev.h"              /* OpenMPDEV N_Vector        */
#elif defined(USE_HIP_NVEC)
#include "nvector/nvector_hip.h"                    /* HIP N_Vector              */
#elif defined(USE_CUDA_NVEC) || defined(USE_CUDAUVM_NVEC)
#include "nvector/nvector_cuda.h"                   /* CUDA N_Vector             */
#else
#include "nvector/nvector_serial.h"                 /* serial N_Vector           */
#include "sunmemory/sunmemory_system.h"
#endif

#include "check_retval.h"

#if defined(USE_CUDA_NVEC) || defined(USE_CUDAUVM_NVEC) || defined(USE_RAJA_NVEC)
#define USE_CUDA
#define USE_GPU
#elif defined(USE_OMPDEV_NVEC)
#define USE_GPU
#elif defined(USE_HIP_NVEC)
#define USE_HIP
#define USE_GPU
#endif

#if defined(USE_CUDA) || defined(USE_HIP)
#define USE_CUDA_OR_HIP
#endif

#if defined(USE_CUDA)
#define GPU_PREFIX(a) cuda##a
#elif defined(USE_HIP)
#define GPU_PREFIX(a) hip##a
#endif

using XYZ_KERNEL_SERIAL_POLICY =
    RAJA::KernelPolicy< RAJA::statement::For<2, RAJA::seq_exec,
                          RAJA::statement::For<1, RAJA::seq_exec,
                            RAJA::statement::For<0, RAJA::seq_exec,
                              RAJA::statement::Lambda<0>
                            >
                          >
                        >
                      >;

#ifdef USE_RAJA_NVEC
#define NVECTOR_ID_STRING "RAJA-unmanaged"
using EXEC_POLICY = RAJA::cuda_exec< 256, false >;
constexpr auto LocalNvector = N_VNew_Raja;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_Raja;

#elif USE_OMPDEV_NVEC
#define NVECTOR_ID_STRING "OpenMPDEV"
using EXEC_POLICY = RAJA::omp_target_parallel_for_exec< 256 >;
constexpr auto LocalNvector = N_VNew_OpenMPDEV;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_OpenMPDEV;

#elif USE_CUDAUVM_NVEC
#define NVECTOR_ID_STRING "CUDA-managed"
using EXEC_POLICY = RAJA::cuda_exec< 256, false >;
using XYZ_KERNEL_POL =
    RAJA::KernelPolicy<
                        RAJA::CudaKernel<
                          RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
                            RAJA::statement::For<1, RAJA::cuda_thread_y_loop,
                              RAJA::statement::For<0, RAJA::cuda_thread_z_loop,
                                RAJA::statement::Lambda<0>
                              >
                            >
                          >
                        >
                      >;

constexpr auto LocalNvector = N_VNewManaged_Cuda;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_Cuda;

#elif USE_CUDA_NVEC
#define NVECTOR_ID_STRING "CUDA-unmanaged"
using EXEC_POLICY = RAJA::cuda_exec< 256, false >;
using XYZ_KERNEL_POL =
    RAJA::KernelPolicy<
                        RAJA::statement::CudaKernel<
                          RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
                            RAJA::statement::For<1, RAJA::cuda_thread_y_loop,
                              RAJA::statement::For<0, RAJA::cuda_thread_z_loop,
                                RAJA::statement::Lambda<0>
                              >
                            >
                          >
                        >
                      >;


constexpr auto LocalNvector = N_VNew_Cuda;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_Cuda;

#elif USE_HIP_NVEC
#define NVECTOR_ID_STRING "HIP-unmanaged"
using EXEC_POLICY = RAJA::hip_exec< 512, false >;
constexpr auto LocalNvector = N_VNew_Hip;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_Hip;

using XYZ_KERNEL_POL =
    RAJA::KernelPolicy<
                        RAJA::statement::HipKernel<
                          RAJA::statement::For<2, RAJA::hip_thread_x_loop,
                            RAJA::statement::For<1, RAJA::hip_thread_y_loop,
                              RAJA::statement::For<0, RAJA::hip_thread_z_loop,
                                RAJA::statement::Lambda<0>
                              >
                            >
                          >
                        >
                      >;

#else
#define NVECTOR_ID_STRING "Serial"
using EXEC_POLICY = RAJA::seq_exec;
using XYZ_KERNEL_POL =
    RAJA::KernelPolicy< RAJA::statement::For<2, RAJA::loop_exec,
                          RAJA::statement::For<1, RAJA::loop_exec,
                            RAJA::statement::For<0, RAJA::loop_exec,
                              RAJA::statement::Lambda<0>
                            >
                          >
                        >
                      >;
constexpr auto LocalNvector = N_VNew_Serial;
#define CopyVecFromDevice(v)

#endif // USE_RAJA_NVEC

#ifdef USE_CUDA_OR_HIP
#define DEVICE_FUNC __device__
#define GPU_SAFE_CALL(val) { gpuAssert((val), __FILE__, __LINE__, 1); }
#else
#define DEVICE_FUNC
#define GPU_SAFE_CALL
#endif


/* Get the vector data array pointer for the device
   if using the GPU, or host if not. */
static realtype* GetVecData(N_Vector y)
{
#ifdef USE_GPU
  if (N_VGetVectorID(y) == SUNDIALS_NVEC_MPIPLUSX)
    return N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(y));
  else
    return N_VGetDeviceArrayPointer(y);
#else
  if (N_VGetVectorID(y) == SUNDIALS_NVEC_MPIPLUSX)
    return N_VGetArrayPointer(N_VGetLocalVector_MPIPlusX(y));
  else
    return N_VGetArrayPointer(y);
#endif
}

/* Turn on fused vector ops for y */
static int EnableFusedVectorOps(N_Vector y)
{
  int myid;
  int retval = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  retval = N_VEnableFusedOps_MPIPlusX(y, 1);
  if (check_retval(&retval, "N_VEnableFusedOps_MPIPlusX", 1, myid)) return(-1);

#if defined(USE_CUDA_NVEC) || defined(USE_CUDAUVM_NVEC)
  retval = N_VEnableFusedOps_Cuda(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_Cuda", 1, myid)) return(-1);
#elif defined(USE_HIP_NVEC)
  retval = N_VEnableFusedOps_Hip(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_Hip", 1, myid)) return(-1);
#elif defined(USE_RAJA_NVEC)
  retval = N_VEnableFusedOps_Raja(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_Raja", 1, myid)) return(-1);
#elif defined(USE_OMPDEV_NVEC)
  retval = N_VEnableFusedOps_OpenMPDEV(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_OpenMPDEV", 1, myid)) return(-1);
#else
  retval = N_VEnableFusedOps_Serial(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_Serial", 1, myid)) return(-1);
#endif

  return(0);
}

#ifdef USE_CUDA_OR_HIP
static void gpuAssert(GPU_PREFIX(Error_t) code, const char *file, int line, int abort)
{
   if (code != GPU_PREFIX(Success))
   {
      fprintf(stderr, "GPU ERROR: %s %s %d\n", GPU_PREFIX(GetErrorString)(code), file, line);
      if (abort) MPI_Abort(MPI_COMM_WORLD, -1);
   }
}
#endif


#endif
