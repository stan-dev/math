# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Configures the SUNDIALS config header files:
#  sundials_config.h and sundials_fconfig.h
# ---------------------------------------------------------------

# ============================================================================
# Generate macros and substitution variables related to the build options.
# ============================================================================

# prepare substitution variable PRECISION_LEVEL for sundials_config.h
string(TOUPPER ${SUNDIALS_PRECISION} SUNDIALS_PRECISION)
set(PRECISION_LEVEL "#define SUNDIALS_${SUNDIALS_PRECISION}_PRECISION 1")

# prepare substitution variable INDEX_TYPE for sundials_config.h
set(INDEX_TYPE "#define SUNDIALS_INT${SUNDIALS_INDEX_SIZE}_T 1")

# Prepare substitution variable SUNDIALS_EXPORT for sundials_config.h
# When building shared SUNDIALS libraries under Windows, use
#      #define SUNDIALS_EXPORT __declspec(dllexport)
# When linking to shared SUNDIALS libraries under Windows, use
#      #define SUNDIALS_EXPORT __declspec(dllimport)
# In all other cases (other platforms or static libraries
# under Windows), the SUNDIALS_EXPORT macro is empty.
# See https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/BuildingWinDLL
if(BUILD_SHARED_LIBS)
  set(SUNDIALS_EXPORT_MACRO
"#if defined(_WIN32)
  #if defined(SUNDIALS_EXPORT)
    #define SUNDIALS_EXPORT __declspec(dllexport)
  #else
    #define SUNDIALS_EXPORT __declspec(dllimport)
  #endif
#else
  #define SUNDIALS_EXPORT
#endif")
endif()

# prepare substitution variable SUNDIALS_USE_GENERIC_MATH for sundials_config.h
if(USE_GENERIC_MATH)
  set(SUNDIALS_USE_GENERIC_MATH TRUE)
endif()

# ============================================================================
# Generate macros and substitution variables related to TPLs
# that SUNDIALS is being built with.
# ============================================================================

# prepare substitution variables for modules that have been built
set(SUNDIALS_CONFIGH_BUILDS "")
foreach(_item ${SUNDIALS_BUILD_LIST})
  if(${${_item}})
    string(REPLACE "BUILD_" "" _module ${_item})
    string(APPEND SUNDIALS_CONFIGH_BUILDS "#define SUNDIALS_${_module} 1\n")
  endif()
endforeach()

# prepare substitution variable SUNDIALS_MPI_ENABLED for sundials_config.h
if(ENABLE_MPI)
  set(SUNDIALS_MPI_ENABLED TRUE)
endif()

# prepare substitution variable SUNDIALS_TRILINOS_HAVE_MPI for sundials_config.h
if(Trilinos_MPI)
  set(SUNDIALS_TRILINOS_HAVE_MPI TRUE)
endif()

# prepare substitution variable(s) SUNDIALS_RAJA_BACKENDS_*
foreach(backend ${SUNDIALS_RAJA_BACKENDS})
  set(SUNDIALS_RAJA_BACKENDS_${backend} TRUE)
endforeach()

# prepare substitution variable(s) SUNDIALS_MAGMA_BACKENDS_*
foreach(backend ${SUNDIALS_MAGMA_BACKENDS})
  set(SUNDIALS_MAGMA_BACKENDS_${backend} TRUE)
endforeach()

# prepare substitution variable SUNDIALS_HAVE_POSIX_TIMERS for sundials_config.h
if(POSIX_TIMERS_TEST_OK) # set in SundialsPOSIXTimers.cmake
  set(SUNDIALS_HAVE_POSIX_TIMERS TRUE)
endif()

# ============================================================================
# Generate macros and substitution variables for the FCMIX interface.
# ============================================================================

# prepare substitution variable FPRECISION_LEVEL for sundials_fconfig.h
if(SUNDIALS_PRECISION MATCHES "SINGLE")
  set(FPRECISION_LEVEL "4")
endif(SUNDIALS_PRECISION MATCHES "SINGLE")
if(SUNDIALS_PRECISION MATCHES "DOUBLE")
  set(FPRECISION_LEVEL "8")
endif(SUNDIALS_PRECISION MATCHES "DOUBLE")
if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  set(FPRECISION_LEVEL "16")
endif(SUNDIALS_PRECISION MATCHES "EXTENDED")

# always define FMPI_COMM_F2C substitution variable in sundials_fconfig.h
set(F77_MPI_COMM_F2C "#define SUNDIALS_MPI_COMM_F2C 1")
set(FMPI_COMM_F2C ".true.")

# =============================================================================
# All required substitution variables should be available at this point.
# Generate the header file and place it in the binary dir.
# =============================================================================

configure_file(
  ${PROJECT_SOURCE_DIR}/include/sundials/sundials_config.in
  ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  )
configure_file(
  ${PROJECT_SOURCE_DIR}/include/sundials/sundials_fconfig.in
  ${PROJECT_BINARY_DIR}/include/sundials/sundials_fconfig.h
  )
