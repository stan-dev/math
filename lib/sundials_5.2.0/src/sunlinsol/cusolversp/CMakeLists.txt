# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the cuda cuSolverSp SUNLinearSolver
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_CUSOLVERSP\n\")")

# Add variable cusolversp_SOURCES with the sources for the cusolversp lib
set(cusolversp_SOURCES sunlinsol_cusolversp_batchqr.cu)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECCUDA library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  )

# Add variable cusolversp_HEADERS with the exported cusolversp header files
set(cusolversp_HEADERS sunlinsol_cusolversp_batchqr.h)
add_prefix(${sundials_SOURCE_DIR}/include/sunlinsol/ cusolversp_HEADERS)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

if(BUILD_STATIC_LIBS)
  add_library(sundials_sunlinsolcusolversp_static STATIC ${cusolversp_SOURCES} ${shared_SOURCES})
  target_link_libraries(sundials_sunlinsolcusolversp_static
                        PUBLIC
                        sundials_nveccuda_static
                        sundials_sunmatrixcusparse_static
                        ${CUDA_cusolver_LIBRARY}
                        PRIVATE
                        ${CUDA_cusparse_LIBRARY})
  target_include_directories(sundials_sunlinsolcusolversp_static
                             PRIVATE . ${sundials_SOURCE_DIR}/src/sundials)
  set_target_properties(sundials_sunlinsolcusolversp_static PROPERTIES
                        OUTPUT_NAME sundials_sunlinsolcusolversp
                        CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_sunlinsolcusolversp_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

if(BUILD_SHARED_LIBS)
  add_library(sundials_sunlinsolcusolversp_shared SHARED ${cusolversp_SOURCES} ${shared_SOURCES})
  target_link_libraries(sundials_sunlinsolcusolversp_shared
                        PUBLIC
                        sundials_nveccuda_shared
                        sundials_sunmatrixcusparse_shared
                        ${CUDA_cusolver_LIBRARY}
                        PRIVATE
                        ${CUDA_cusparse_LIBRARY})
  target_include_directories(sundials_sunlinsolcusolversp_shared
                             PRIVATE . ${sundials_SOURCE_DIR}/src/sundials)
  set_target_properties(sundials_sunlinsolcusolversp_shared PROPERTIES
                        OUTPUT_NAME sundials_sunlinsolcusolversp
                        CLEAN_DIRECT_OUTPUT 1
                        VERSION ${sunlinsollib_VERSION}
                        SOVERSION ${sunlinsollib_SOVERSION})
  install(TARGETS sundials_sunlinsolcusolversp_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the CUDA NVector header files
install(FILES ${cusolversp_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_CUSOLVERSP module")
