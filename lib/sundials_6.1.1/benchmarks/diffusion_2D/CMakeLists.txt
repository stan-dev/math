# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

if(BUILD_ARKODE OR BUILD_CVODE OR BUILD_IDA)

  # Shared sources
  set(shared_sources
    diffusion_2D.hpp
    diffusion_2D.cpp
    preconditioner_jacobi.cpp)

  # Benchmark prefix
  set(benchmark_prefix ${SUNDIALS_SOURCE_DIR}/benchmarks/diffusion_2D/)

  add_prefix(${benchmark_prefix} shared_sources)

  # MPI + Serial
  add_subdirectory(mpi_serial)

  # MPI + CUDA/HIP
  if(BUILD_NVECTOR_CUDA OR BUILD_NVECTOR_HIP)
    add_subdirectory(mpi_gpu)
  endif()

  install(FILES README.md
    DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/diffusion_2D")

endif()
