# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# SUNMemory level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
endif()

add_subdirectory(system)

if(ENABLE_CUDA)
  add_subdirectory(cuda)
endif()

if(ENABLE_HIP)
  add_subdirectory(hip)
endif()

if(ENABLE_SYCL)
  add_subdirectory(sycl)
endif()
