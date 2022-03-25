# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup MAGMA correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# TPL_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_MAGMA_INCLUDED)
  set(SUNDIALS_MAGMA_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

if(SUNDIALS_PRECISION MATCHES "extended")
  print_error("SUNDIALS MAGMA interface is not compatible with extended precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(MAGMA REQUIRED)

message(STATUS "MAGMA_VERSION:     ${MAGMA_VERSION}")
message(STATUS "MAGMA_LIBRARIES:   ${MAGMA_LIBRARIES}")
message(STATUS "MAGMA_INCLUDE_DIR: ${MAGMA_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA" AND NOT ENABLE_CUDA)
  print_error("SUNDIALS_MAGMA_BACKENDS includes CUDA but CUDA is not enabled. Set ENABLE_CUDA=ON or change the backend.")
endif()
if(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP" AND NOT ENABLE_HIP)
  print_error("SUNDIALS_MAGMA_BACKENDS includes HIP but HIP is not enabled. Set ENABLE_HIP=ON or change the backend.")
endif()

if(MAGMA_FOUND AND (NOT MAGMA_WORKS))
  # Create the MAGMA_TEST directory
  set(MAGMA_TEST_DIR ${PROJECT_BINARY_DIR}/CMakeFiles/MAGMA_TEST)
  file(MAKE_DIRECTORY ${MAGMA_TEST_DIR})

  if(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP")
    set(lang CXX)
    set(ext cxx)
    set(define_have "\#define HAVE_HIP")
    set(lib hip::host)
  elseif(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA")
    set(lang CUDA)
    set(ext cu)
    set(define_have "\#define HAVE_CUBLAS")
    set(lib )
  endif()

  file(WRITE ${MAGMA_TEST_DIR}/ltest.${ext}
  "${define_have}\n"
  "\#include \"magma_v2.h\"\n"
  "int main(){\n"
  "magma_int_t a=0;\n"
  "return(a);\n"
  "}\n")

  try_compile(COMPILE_OK ${MAGMA_TEST_DIR} ${MAGMA_TEST_DIR}/ltest.${ext}
    CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES=${MAGMA_INCLUDE_DIR}"
    LINK_LIBRARIES ${MAGMA_LIBRARIES} ${lib}
    OUTPUT_VARIABLE COMPILE_OUTPUT
    ${lang}_STANDARD ${CMAKE_${lang}_STANDARD}
  )

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if MAGMA works... OK")
    set(MAGMA_WORKS TRUE CACHE BOOL "MAGMA works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if MAGMA works... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to MAGMA is not functional.")
  endif()
elseif(MAGMA_FOUND AND MAGMA_WORKS)
  message(STATUS "Skipped MAGMA tests, assuming MAGMA works with SUNDIALS.")
endif()
