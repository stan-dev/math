# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Radu Serban, and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the CVODE library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall CVODE\n\")")

# Add variable cvode_SOURCES with the sources for the CVODE library
set(cvode_SOURCES
  cvode.c
  cvode_bandpre.c
  cvode_bbdpre.c
  cvode_diag.c
  cvode_direct.c
  cvode_io.c
  cvode_ls.c
  cvode_nls.c
  cvode_proj.c
  cvode_spils.c
  )

# Add variable cvode_HEADERS with the exported CVODE header files
set(cvode_HEADERS
  cvode.h
  cvode_bandpre.h
  cvode_bbdpre.h
  cvode_diag.h
  cvode_direct.h
  cvode_ls.h
  cvode_proj.h
  cvode_spils.h
  )

# Add prefix with complete path to the CVODE header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/cvode/ cvode_HEADERS)

# Build fused kernel libraries
if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)

  if(BUILD_NVECTOR_CUDA)

    set_source_files_properties(cvode_fused_gpu.cpp PROPERTIES
      LANGUAGE CUDA)

    sundials_add_library(sundials_cvode_fused_cuda
      SOURCES
        cvode_fused_gpu.cpp
      COMPILE_DEFINITIONS
        PRIVATE USE_CUDA
      LINK_LIBRARIES
        PRIVATE sundials_nveccuda
      OUTPUT_NAME
        sundials_cvode_fused_cuda
      VERSION
        ${cvodelib_VERSION}
      SOVERSION
        ${cvodelib_SOVERSION}
      )
  endif()

  if(BUILD_NVECTOR_HIP)
    sundials_add_library(sundials_cvode_fused_hip
      SOURCES
        cvode_fused_gpu.cpp
      COMPILE_DEFINITIONS
        PRIVATE USE_HIP
      LINK_LIBRARIES
        PRIVATE sundials_nvechip
      OUTPUT_NAME
        sundials_cvode_fused_hip
      VERSION
        ${cvodelib_VERSION}
      SOVERSION
        ${cvodelib_SOVERSION}
      )
  endif()

  sundials_add_library(sundials_cvode_fused_stubs
    SOURCES
      cvode_fused_stubs.c
    OUTPUT_NAME
      sundials_cvode_fused_stubs
    VERSION
      ${cvodelib_VERSION}
    SOVERSION
      ${cvodelib_SOVERSION}
  )

  set(_fused_link_lib sundials_cvode_fused_stubs)
endif()

if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  set(_include_caliper_if_needed PUBLIC ${caliper_INCLUDE_DIR})
  set(_link_caliper_if_needed PUBLIC caliper)
endif()

# Create the library
sundials_add_library(sundials_cvode
  SOURCES
    ${cvode_SOURCES}
  HEADERS
    ${cvode_HEADERS}
  INCLUDE_SUBDIR
    cvode
  INCLUDE_DIRECTORIES
    ${_include_caliper_if_needed}
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_sunmemsys_obj
    sundials_nvecserial_obj
    sundials_sunmatrixband_obj
    sundials_sunmatrixdense_obj
    sundials_sunmatrixsparse_obj
    sundials_sunlinsolband_obj
    sundials_sunlinsoldense_obj
    sundials_sunlinsolspbcgs_obj
    sundials_sunlinsolspfgmr_obj
    sundials_sunlinsolspgmr_obj
    sundials_sunlinsolsptfqmr_obj
    sundials_sunlinsolpcg_obj
    sundials_sunnonlinsolnewton_obj
    sundials_sunnonlinsolfixedpoint_obj
  LINK_LIBRARIES
    # Link to stubs so examples work.
    PRIVATE ${_fused_link_lib}
    # Link to caliper if conditions are met.
    ${_link_caliper_if_needed}
  OUTPUT_NAME
    sundials_cvode
  VERSION
    ${cvodelib_VERSION}
  SOVERSION
    ${cvodelib_SOVERSION}
)

# Finished CVODE
message(STATUS "Added CVODE module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()
