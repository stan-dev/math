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
 * CUDA/HIP utility functions
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

int DeviceSynchronize()
{
#if defined(USE_HIP)
  hipDeviceSynchronize();
#elif defined(USE_CUDA)
  cudaDeviceSynchronize();
#endif

  return 0;
}

int CopyDataFromDevice(N_Vector y)
{
#if defined(USE_HIP)
  N_VCopyFromDevice_Hip(N_VGetLocalVector_MPIPlusX(y));
#elif defined(USE_CUDA)
  N_VCopyFromDevice_Cuda(N_VGetLocalVector_MPIPlusX(y));
#endif

  return 0;
}
