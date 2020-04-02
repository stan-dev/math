# ---------------------------------------------------------------
# Programmer(s): Radu Serban, David J. Gardner, Cody J. Balos,
#                and Slaven Peles @ LLNL
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
# Top level CMakeLists.txt for SUNDIALS (for cmake build system)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Initial commands
# ---------------------------------------------------------------

# Require CMake 3.5 by default, but certain options will force us
# to require CMake 3.10
cmake_minimum_required(VERSION 3.5)

# Project SUNDIALS (initially only C supported)
# sets PROJECT_SOURCE_DIR and PROJECT_BINARY_DIR variables
project(sundials C)

# Set some variables with info on the SUNDIALS project
set(PACKAGE_BUGREPORT "woodward6@llnl.gov")
set(PACKAGE_NAME "SUNDIALS")
set(PACKAGE_STRING "SUNDIALS 5.2.0")
set(PACKAGE_TARNAME "sundials")

# set SUNDIALS version numbers
# (use "" for the version label if none is needed)
set(PACKAGE_VERSION_MAJOR "5")
set(PACKAGE_VERSION_MINOR "2")
set(PACKAGE_VERSION_PATCH "0")
set(PACKAGE_VERSION_LABEL "")

if(PACKAGE_VERSION_LABEL)
  set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-${PACKAGE_VERSION_LABEL}")
else()
  set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Prohibit in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source build prohibited.")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Hide some cache variables
mark_as_advanced(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)

# Always show the C compiler and flags
mark_as_advanced(CLEAR
  CMAKE_C_COMPILER
  CMAKE_C_FLAGS)

# Specify the VERSION and SOVERSION for shared libraries

set(arkodelib_VERSION "4.2.0")
set(arkodelib_SOVERSION "4")

set(cvodelib_VERSION "5.2.0")
set(cvodelib_SOVERSION "5")

set(cvodeslib_VERSION "5.2.0")
set(cvodeslib_SOVERSION "5")

set(idalib_VERSION "5.2.0")
set(idalib_SOVERSION "5")

set(idaslib_VERSION "4.2.0")
set(idaslib_SOVERSION "4")

set(kinsollib_VERSION "5.2.0")
set(kinsollib_SOVERSION "5")

set(cpodeslib_VERSION "0.0.0")
set(cpodeslib_SOVERSION "0")

set(nveclib_VERSION "5.2.0")
set(nveclib_SOVERSION "5")

set(sunmatrixlib_VERSION "3.2.0")
set(sunmatrixlib_SOVERSION "3")

set(sunlinsollib_VERSION "3.2.0")
set(sunlinsollib_SOVERSION "3")

set(sunnonlinsollib_VERSION "2.2.0")
set(sunnonlinsollib_SOVERSION "2")

# Specify the location of additional CMAKE modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# Get correct build paths automatically, but expose CMAKE_INSTALL_LIBDIR
# as a regular cache variable so that a user can more easily see what
# the library dir was set to be by GNUInstallDirs.
include(GNUInstallDirs)
mark_as_advanced(CLEAR CMAKE_INSTALL_LIBDIR)

# ---------------------------------------------------------------
# Which modules to build?
# ---------------------------------------------------------------

# For each SUNDIALS solver available (i.e. for which we have the
# sources), give the user the option of enabling/disabling it.

if(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/arkode")
  option(BUILD_ARKODE "Build the ARKODE library" ON)
else()
  set(BUILD_ARKODE OFF)
endif()

if(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvode")
  option(BUILD_CVODE "Build the CVODE library" ON)
else()
  set(BUILD_CVODE OFF)
endif()

if(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvodes")
  option(BUILD_CVODES "Build the CVODES library" ON)
else()
  set(BUILD_CVODES OFF)
endif()

if(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/ida")
  option(BUILD_IDA "Build the IDA library" ON)
else()
  set(BUILD_IDA OFF)
endif()

if(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/idas")
  option(BUILD_IDAS "Build the IDAS library" ON)
else()
  set(BUILD_IDAS OFF)
endif()

if(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/kinsol")
  option(BUILD_KINSOL "Build the KINSOL library" ON)
else()
  set(BUILD_KINSOL OFF)
endif()

# CPODES is always OFF for now.  (commented out for Release); ToDo: better way to do this?
#IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cpodes")
#  OPTION(BUILD_CPODES  "Build the CPODES library"  OFF)
#ELSE()
#  SET(BUILD_CPODES OFF)
#ENDIF()

# ---------------------------------------------------------------
# MACRO definitions
# ---------------------------------------------------------------
include(SundialsCMakeMacros)

# ---------------------------------------------------------------
# Check for deprecated SUNDIALS CMake options/variables
# ---------------------------------------------------------------
include(SundialsDeprecated)

# ---------------------------------------------------------------
# xSDK specific options
# ---------------------------------------------------------------
include(SundialsXSDK)

# ---------------------------------------------------------------
# Build specific C flags
# ---------------------------------------------------------------

# Hide all build type specific flags
mark_as_advanced(FORCE
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_RELWITHDEBINFO)

# Only show flags for the current build type if it is set
# NOTE: Build specific flags are appended those in CMAKE_C_FLAGS
if(CMAKE_BUILD_TYPE)
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    message("Appending C debug flags")
    mark_as_advanced(CLEAR CMAKE_C_FLAGS_DEBUG)
  elseif(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
    message("Appending C min size release flags")
    mark_as_advanced(CLEAR CMAKE_C_FLAGS_MINSIZEREL)
  elseif(CMAKE_BUILD_TYPE MATCHES "Release")
    message("Appending C release flags")
    mark_as_advanced(CLEAR CMAKE_C_FLAGS_RELEASE)
  elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    message("Appending C release with debug info flags")
    mark_as_advanced(CLEAR CMAKE_C_FLAGS_RELWITHDEBINFO)
  endif()
endif()

# ---------------------------------------------------------------
# Option to specify precision (realtype)
# ---------------------------------------------------------------

set(DOCSTR "single, double, or extended")
show_variable(SUNDIALS_PRECISION STRING "${DOCSTR}" "double")

# prepare substitution variable PRECISION_LEVEL for sundials_config.h
string(TOUPPER ${SUNDIALS_PRECISION} SUNDIALS_PRECISION)
set(PRECISION_LEVEL "#define SUNDIALS_${SUNDIALS_PRECISION}_PRECISION 1")

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

# ---------------------------------------------------------------
# Option to specify index type
# ---------------------------------------------------------------

set(DOCSTR "Signed 64-bit (64) or signed 32-bit (32) integer")
show_variable(SUNDIALS_INDEX_SIZE STRING "${DOCSTR}" "64")
set(DOCSTR "Integer type to use for indices in SUNDIALS")
show_variable(SUNDIALS_INDEX_TYPE STRING "${DOCSTR}" "")
mark_as_advanced(SUNDIALS_INDEX_TYPE)
include(SundialsIndexSize)

# ---------------------------------------------------------------
# Enable Fortran interface?
# ---------------------------------------------------------------

# Fortran interface is disabled by default
set(DOCSTR "Enable Fortran 77 interfaces")
option(F77_INTERFACE_ENABLE "${DOCSTR}" OFF)

# Check that at least one solver with a Fortran 77 interface is built
if(NOT BUILD_ARKODE AND NOT BUILD_CVODE AND NOT BUILD_IDA AND NOT BUILD_KINSOL)
  if(F77_INTERFACE_ENABLE)
    print_warning("Enabled packages do not support Fortran 77 interface" "Disabling F77 interface")
    force_variable(F77_INTERFACE_ENABLE BOOL "${DOCSTR}" OFF)
  endif()
  hide_variable(F77_INTERFACE_ENABLE)
endif()

# Fortran 2003 interface is disabled by default
set(DOCSTR "Enable Fortran 2003 interfaces")
option(F2003_INTERFACE_ENABLE "${DOCSTR}" OFF)

if(F2003_INTERFACE_ENABLE)
  # F2003 interface only supports double precision
  if(NOT (SUNDIALS_PRECISION MATCHES "DOUBLE"))
    print_warning("F2003 interface is not compatible with ${SUNDIALS_PRECISION} precision"
                  "Disabling F2003 interface")
    force_variable(F2003_INTERFACE_ENABLE BOOL "${DOCSTR}" OFF)
  endif()

  # F2003 interface only supports 64-bit indices
  if(NOT (SUNDIALS_INDEX_SIZE MATCHES "64"))
    print_warning("F2003 interface is not compatible with ${SUNDIALS_INDEX_SIZE}-bit indicies"
                  "Disabling F2003 interface")
    force_variable(F2003_INTERFACE_ENABLE BOOL "${DOCSTR}" OFF)
  endif()

  # Put all F2003 modules into one build directory
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran")

  # Allow a user to set where the Fortran modules will be installed
  set(DOCSTR "Directory where Fortran module files are installed")
  show_variable(Fortran_INSTALL_MODDIR STRING "${DOCSTR}" "fortran")
endif()

# ---------------------------------------------------------------
# Options to build static and/or shared libraries
# ---------------------------------------------------------------

option(BUILD_STATIC_LIBS "Build static libraries" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Prepare substitution variable SUNDIALS_EXPORT for sundials_config.h
# When building shared SUNDIALS libraries under Windows, use
#      #define SUNDIALS_EXPORT __declspec(dllexport)
# When linking to shared SUNDIALS libraries under Windows, use
#      #define SUNDIALS_EXPORT __declspec(dllimport)
# In all other cases (other platforms or static libraries
# under Windows), the SUNDIALS_EXPORT macro is empty

if(BUILD_SHARED_LIBS AND WIN32)
  set(SUNDIALS_EXPORT
    "#ifdef BUILD_SUNDIALS_LIBRARY
#define SUNDIALS_EXPORT __declspec(dllexport)
#else
#define SUNDIALS_EXPORT __declspec(dllimport)
#endif")
else(BUILD_SHARED_LIBS AND WIN32)
  set(SUNDIALS_EXPORT "#define SUNDIALS_EXPORT")
endif(BUILD_SHARED_LIBS AND WIN32)

# Make sure we build at least one type of libraries
if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
  print_warning("Both static and shared library generation were disabled"
                "Building static libraries was re-enabled")
  force_variable(BUILD_STATIC_LIBS BOOL "Build static libraries" ON)
endif(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)

# ---------------------------------------------------------------
# Option to use the generic math libraries (UNIX only)
# ---------------------------------------------------------------

if(UNIX)
  option(USE_GENERIC_MATH "Use generic (std-c) math libraries" ON)
  # executables will be linked against -lm
  set(EXTRA_LINK_LIBS -lm)
  if(USE_GENERIC_MATH)
    # prepare substitution variable for sundials_config.h
    set(SUNDIALS_USE_GENERIC_MATH TRUE)
  endif(USE_GENERIC_MATH)
endif(UNIX)

# ---------------------------------------------------------------
# Check for POSIX timers
# ---------------------------------------------------------------
include(SundialsPOSIXTimers)

# ===============================================================
# Options for languages and parallelism
# ===============================================================

# ---------------------------------------------------------------
# Enable MPI support?
# ---------------------------------------------------------------
option(MPI_ENABLE "Enable MPI support" OFF)

# ---------------------------------------------------------------
# Enable OpenMP support?
# ---------------------------------------------------------------
option(OPENMP_ENABLE "Enable OpenMP support" OFF)

# provide OPENMP_DEVICE_ENABLE option
option(OPENMP_DEVICE_ENABLE "Enable OpenMP device offloading support" OFF)

# Advanced option to skip OpenMP device offloading support check.
# This is needed for a specific compiler that doesn't correctly
# report its OpenMP spec date (with CMake >= 3.9).
option(SKIP_OPENMP_DEVICE_CHECK "Skip the OpenMP device offloading support check" OFF)
mark_as_advanced(FORCE SKIP_OPENMP_DEVICE_CHECK)

if(OPENMP_DEVICE_ENABLE)
  if(CMAKE_VERSION VERSION_LESS 3.10)
    print_error("CMake >= 3.10 is required for OPENMP_DEVICE_ENABLE.")
  endif()
  cmake_minimum_required(VERSION 3.10)
endif()

# ---------------------------------------------------------------
# Enable Pthread support?
# ---------------------------------------------------------------
option(PTHREAD_ENABLE "Enable Pthreads support" OFF)

# -------------------------------------------------------------
# Enable CUDA support?
# -------------------------------------------------------------
option(CUDA_ENABLE "Enable CUDA support" OFF)

if(CUDA_ENABLE)
  if(CMAKE_VERSION VERSION_LESS 3.10)
    print_error("CMake >= 3.10 is required for CUDA_ENABLE.")
  endif()
  cmake_minimum_required(VERSION 3.10)
endif()

# ===============================================================
# Options for external packages
# ===============================================================

# ---------------------------------------------------------------
# Enable LAPACK support?
# ---------------------------------------------------------------
option(LAPACK_ENABLE "Enable Lapack support" OFF)

# LAPACK does not support extended precision
if(LAPACK_ENABLE AND SUNDIALS_PRECISION MATCHES "EXTENDED")
  print_error("LAPACK is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# ---------------------------------------------------------------
# Enable SuperLU_DIST support?
# ---------------------------------------------------------------
sundials_option(SUPERLUDIST_ENABLE BOOL "Enable SuperLU_DIST support" OFF
                XSDK_NAME TPL_ENABLE_SUPERLUDIST)

sundials_option(SUPERLUDIST_INCLUDE_DIR PATH "SuperLU_DIST include directory" "${SUPERLUDIST_INCLUDE_DIR}"
                DEPENDS_ON SUPERLUDIST_ENABLE
                XSDK_NAME TPL_SUPERLUDIST_INCLUDE_DIRS)

sundials_option(SUPERLUDIST_LIBRARY_DIR PATH "SuperLU_DIST library directory" "${SUPERLUDIST_LIBRARY_DIR}"
                DEPENDS_ON SUPERLUDIST_ENABLE
                XSDK_HIDE)

sundials_option(SUPERLUDIST_LIBRARIES STRING "Semi-colon separated list of additional libraries needed for SuperLU_DIST." "${SUPERLUDIST_LIBRARIES}"
                DEPENDS_ON SUPERLUDIST_ENABLE
                XSDK_NAME TPL_SUPERLUDIST_LIBRARIES)

sundials_option(SUPERLUDIST_OpenMP BOOL "Enable SUNDIALS support for SuperLU_DIST OpenMP on-node parallelism" OFF
                DEPENDS_ON SUPERLUDIST_ENABLE
                XSDK_NAME TPL_SUPERLUDIST_OpenMP)

# ---------------------------------------------------------------
# Enable SuperLU_MT support?
# ---------------------------------------------------------------
sundials_option(SUPERLUMT_ENABLE BOOL "Enable SuperLU_MT support" OFF
                XSDK_NAME TPL_ENABLE_SUPERLUMT)

sundials_option(SUPERLUMT_INCLUDE_DIR PATH "SuperLU_MT include directory" "${SUPERLUMT_INCLUDE_DIR}"
                DEPENDS_ON SUPERLUMT_ENABLE
                XSDK_NAME TPL_SUPERLUMT_INCLUDE_DIRS)

sundials_option(SUPERLUMT_LIBRARY_DIR PATH "SuperLU_MT library directory" "${SUPERLUMT_LIBRARY_DIR}"
                DEPENDS_ON SUPERLUMT_ENABLE
                XSDK_HIDE)

sundials_option(SUPERLUMT_LIBRARIES STRING "Semi-colon separated list of additional libraries needed for SuperLU_MT." "${SUPERLUMT_LIBRARIES}"
                DEPENDS_ON SUPERLUMT_ENABLE
                XSDK_NAME TPL_SUPERLUMT_LIBRARIES)

sundials_option(SUPERLUMT_THREAD_TYPE BOOL "SuperLU_MT threading type: OPENMP or PTHREAD" "PTHREAD"
                DEPENDS_ON SUPERLUMT_ENABLE
                XSDK_NAME TPL_SUPERLUMT_THREAD_TYPE)

# ---------------------------------------------------------------
# Enable KLU support?
# ---------------------------------------------------------------
option(KLU_ENABLE "Enable KLU support" OFF)

# KLU does not support single or extended precision
if(KLU_ENABLE AND
  (SUNDIALS_PRECISION MATCHES "SINGLE" OR SUNDIALS_PRECISION MATCHES "EXTENDED"))
  print_warning("KLU is not compatible with ${SUNDIALS_PRECISION} precision"
                "Disabling KLU")
  force_variable(KLU_ENABLE BOOL "KLU is disabled" OFF)
endif()

# ---------------------------------------------------------------
# Enable hypre Vector support?
# ---------------------------------------------------------------
option(HYPRE_ENABLE "Enable hypre support" OFF)

# Using hypre requres building with MPI enabled
if(HYPRE_ENABLE AND NOT MPI_ENABLE)
  print_warning("MPI not enabled - Disabling hypre"
                "Set MPI_ENABLE to ON to use parhyp")
  force_variable(HYPRE_ENABLE BOOL "Enable hypre support" OFF)
endif()

# ---------------------------------------------------------------
# Enable PETSc support?
# ---------------------------------------------------------------

sundials_option(PETSC_ENABLE BOOL "Enable PETSc support" OFF
                XSDK_NAME TPL_ENABLE_PETSC)

sundials_option(PETSC_DIR PATH "Path to the root of a PETSc installation" "${PETSC_DIR}"
                DEPENDS_ON PETSC_ENABLE
                XSDK_NAME TPL_PETSC_DIR)

sundials_option(PETSC_ARCH STRING "PETSc architecture (optional)" "${PETSC_ARCH}"
                DEPENDS_ON PETSC_ENABLE
                XSDK_NAME TPL_PETSC_ARCH)

sundials_option(PETSC_LIBRARIES STRING "Semi-colon separated list of PETSc link libraries" "${PETSC_LIBRARIES}"
                DEPENDS_ON PETSC_ENABLE
                ADVANCED)

sundials_option(PETSC_INCLUDES STRING "Semi-colon separated list of PETSc include directories" "${PETSC_INCLUDES}"
                DEPENDS_ON PETSC_ENABLE
                ADVANCED)

sundials_option(PETSC_WORKS BOOL "Set to ON to force CMake to accept a given PETSc configuration" OFF
                DEPENDS_ON PETSC_ENABLE
                ADVANCED)

# Using PETSc requres building with MPI enabled
if(PETSC_ENABLE AND NOT MPI_ENABLE)
  print_error("MPI is required for PETSc support. Set MPI_ENABLE to ON.")
endif()


# ---------------------------------------------------------------
# Enable Trilinos support?
# ---------------------------------------------------------------
option(Trilinos_ENABLE "Enable Trilinos support" OFF)

# -------------------------------------------------------------
# Enable RAJA support?
# -------------------------------------------------------------
sundials_option(RAJA_ENABLE BOOL "Enable RAJA support" OFF
                XSDK_NAME TPL_ENABLE_RAJA)

if(RAJA_ENABLE AND (NOT CUDA_ENABLE))
  print_error("CUDA is required for RAJA support. Please enable CUDA and RAJA.")
endif()

# ===============================================================
# Options for examples
# ===============================================================

# ---------------------------------------------------------------
# Enable examples?
# ---------------------------------------------------------------

# Enable C examples (on by default)
option(EXAMPLES_ENABLE_C "Build SUNDIALS C examples" ON)

# C++ examples (off by default, unless Trilinos is enabled)
set(DOCSTR "Build C++ examples")
option(EXAMPLES_ENABLE_CXX "${DOCSTR}" ${Trilinos_ENABLE})

# F77/F90 examples (on by default) are an option only if
# the Fortran 77 interface is enabled.
if(F77_INTERFACE_ENABLE)

  set(DOCSTR "Build SUNDIALS FORTRAN 77 examples")
  show_variable(EXAMPLES_ENABLE_F77 BOOL "${DOCSTR}" ON)
  # Fortran 77 examples do not support single or extended precision
  if(EXAMPLES_ENABLE_F77 AND (SUNDIALS_PRECISION MATCHES "EXTENDED" OR SUNDIALS_PRECISION MATCHES "SINGLE"))
    print_warning("F77 examples are not compatible with ${SUNDIALS_PRECISION} precision"
                  "EXAMPLES_ENABLE_F77")
    force_variable(EXAMPLES_ENABLE_F77 BOOL "${DOCSTR}" OFF)
  endif()

  set(DOCSTR "Build SUNDIALS FORTRAN 90 examples")
  show_variable(EXAMPLES_ENABLE_F90 BOOL "${DOCSTR}" ON)
  # Fortran 90 examples do not support extended precision
  if(EXAMPLES_ENABLE_F90 AND (SUNDIALS_PRECISION MATCHES "EXTENDED"))
    print_warning("F90 examples are not compatible with ${SUNDIALS_PRECISION} precision"
                  "Disabling EXAMPLES_ENABLE_F90")
    force_variable(EXAMPLES_ENABLE_F90 BOOL "${DOCSTR}" OFF)
  endif()

else()

  # set back to OFF (in case was ON)
  if(EXAMPLES_ENABLE_F77)
    print_warning("EXAMPLES_ENABLE_F77 is ON but F77_INTERFACE_ENABLE is OFF"
                  "Disabling EXAMPLES_ENABLE_F77")
    force_variable(EXAMPLES_ENABLE_F77 BOOL "${DOCSTR}" OFF)
  endif()
  hide_variable(EXAMPLES_ENABLE_F77)

  # set back to OFF (in case was ON)
  if(EXAMPLES_ENABLE_F90)
    print_warning("EXAMPLES_ENABLE_F90 is ON but F77_INTERFACE_ENABLE is  OFF"
                  "Disabling EXAMPLES_ENABLE_F90")
    force_variable(EXAMPLES_ENABLE_F90 BOOL "${DOCSTR}" OFF)
  endif()
  hide_variable(EXAMPLES_ENABLE_F90)

endif()

# F2003 examples (on by default) are an option only if the
# Fortran 2003 interface is enabled.
if(F2003_INTERFACE_ENABLE)

  set(DOCSTR "Build SUNDIALS Fortran 2003 examples")
  show_variable(EXAMPLES_ENABLE_F2003 BOOL "${DOCSTR}" ON)

  # Fortran 2003 examples only support double precision
  if(EXAMPLES_ENABLE_F2003 AND (NOT (SUNDIALS_PRECISION MATCHES "DOUBLE")))
    print_warning("F2003 examples are not compatible with ${SUNDIALS_PRECISION} precision"
                  "Disabling EXAMPLES_ENABLE_F2003")
    force_variable(EXAMPLES_ENABLE_F2003 BOOL "${DOCSTR}" OFF)
  endif()

  # Fortran 2003 examples only support 64-bit indices
  if(EXAMPLES_ENABLE_F2003 AND (NOT (SUNDIALS_INDEX_SIZE MATCHES "64")))
    print_warning("F2003 examples are not compatible with ${SUNDIALS_INDEX_SIZE}-bit indices"
                  "Disabling EXAMPLES_ENABLE_F2003")
    force_variable(EXAMPLES_ENABLE_F2003 BOOL "${DOCSTR}" OFF)
  endif()

else()

  # set back to OFF (in case was ON)
  if(EXAMPLES_ENABLE_F2003)
    print_warning("EXAMPLES_ENABLE_F2003 is ON but F2003_INTERFACE_ENABLE is OFF"
                  "Disabling EXAMPLES_ENABLE_F2003")
    force_variable(EXAMPLES_ENABLE_F2003 BOOL "${DOCSTR}" OFF)
  endif()
  hide_variable(EXAMPLES_ENABLE_F2003)

endif()

# CUDA examples
set(DOCSTR "Build SUNDIALS CUDA examples")
if(CUDA_ENABLE)
  option(EXAMPLES_ENABLE_CUDA "${DOCSTR}" ON)
else()
  if(EXAMPLES_ENABLE_CUDA)
    print_warning("EXAMPLES_ENABLE_CUDA is ON but CUDA_ENABLE is OFF"
                  "Disabling EXAMPLES_ENABLE_CUDA")
    force_variable(EXAMPLES_ENABLE_CUDA BOOL "${DOCSTR}" OFF)
  endif()
endif()

# If any of the above examples are enabled set EXAMPLES_ENABLED to TRUE
if(EXAMPLES_ENABLE_C OR
    EXAMPLES_ENABLE_F77 OR
    EXAMPLES_ENABLE_CXX OR
    EXAMPLES_ENABLE_F90 OR
    EXAMPLES_ENABLE_F2003 OR
    EXAMPLES_ENABLE_CUDA)
  set(EXAMPLES_ENABLED TRUE)
else()
  set(EXAMPLES_ENABLED FALSE)
endif()

# ---------------------------------------------------------------
# Install examples?
# ---------------------------------------------------------------

# Enable installing examples by default
set(DOCSTR "Install SUNDIALS examples")
if(EXAMPLES_ENABLED)
  option(EXAMPLES_INSTALL "${DOCSTR}" ON)
else()
  force_variable(EXAMPLES_INSTALL BOOL "${DOCSTR}" OFF)
  hide_variable(EXAMPLES_INSTALL)
endif()

# If examples are to be exported, check where we should install them.
if(EXAMPLES_INSTALL)

  show_variable(EXAMPLES_INSTALL_PATH PATH
    "Output directory for installing example files"
    "${CMAKE_INSTALL_PREFIX}/examples")

  if(NOT EXAMPLES_INSTALL_PATH)
    print_warning("The example installation path is empty"
      "Example installation path was reset to its default value")
    set(EXAMPLES_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/examples" CACHE STRING
      "Output directory for installing example files" FORCE)
  endif()

else()

  hide_variable(EXAMPLES_INSTALL_PATH)

endif()


# ==============================================================================
# Advanced (hidden) options
# ==============================================================================

# ------------------------------------------------------------------------------
# Manually specify the Fortran name-mangling scheme
#
# The build system tries to infer the Fortran name-mangling scheme using a
# Fortran compiler and defaults to using lower case and one underscore if the
# scheme can not be determined. If a working Fortran compiler is not available
# or the user needs to override the inferred or default scheme, the following
# options specify the case and number of appended underscores corresponding to
# the Fortran name-mangling scheme of symbol names that do not themselves
# contain underscores. This is all we really need for the FCMIX and LAPACK
# interfaces. A working Fortran compiler is only necessary for building Fortran
# example programs.
# ------------------------------------------------------------------------------

# The case to use in the name-mangling scheme
show_variable(SUNDIALS_F77_FUNC_CASE STRING
  "case of Fortran function names (lower/upper)"
  "")

# The number of underscores of appended in the name-mangling scheme
show_variable(SUNDIALS_F77_FUNC_UNDERSCORES STRING
  "number of underscores appended to Fortran function names (none/one/two)"
  "")

# Hide the name-mangling varibales as advanced options
mark_as_advanced(FORCE SUNDIALS_F77_FUNC_CASE)
mark_as_advanced(FORCE SUNDIALS_F77_FUNC_UNDERSCORES)

# If used, both case and underscores must be set
if((NOT SUNDIALS_F77_FUNC_CASE) AND SUNDIALS_F77_FUNC_UNDERSCORES)
  message(FATAL_ERROR
    "If SUNDIALS_F77_FUNC_UNDERSCORES is set, SUNDIALS_F77_FUNC_CASE must also be set.")
endif()

if(SUNDIALS_F77_FUNC_CASE AND (NOT SUNDIALS_F77_FUNC_UNDERSCORES))
  message(FATAL_ERROR
    "If SUNDIALS_F77_FUNC_CASE is set, SUNDIALS_F77_FUNC_UNDERSCORES must also be set.")
endif()

# ------------------------------------------------------------------------------
# Advanced testing options
#
# NOTE: Development examples are currently used for internal testing and may
# produce erroneous failures when run on different systems as the pass/fail
# status is determined by comparing the output against a saved output file.
# ------------------------------------------------------------------------------

# Include development examples in regression tests
option(SUNDIALS_TEST_DEVTESTS "Include development tests in make test" OFF)
mark_as_advanced(FORCE SUNDIALS_TEST_DEVTESTS)

# Include unit tests in regression tests
option(SUNDIALS_TEST_UNITTESTS "Include unit tests in make test" OFF)
mark_as_advanced(FORCE SUNDIALS_TEST_UNITTESTS)

# ===============================================================
# Add any platform specifc settings
# ===============================================================

# Under Windows, add compiler directive to inhibit warnings
# about use of unsecure functions

if(WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

if(APPLE)
  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
endif(APPLE)

# ===============================================================
# Fortran and C++ settings
# ===============================================================

# ---------------------------------------------------------------
# A Fortran compiler is needed to:
# (a) Determine compiler the name-mangling scheme if FCMIX or
#     LAPACK are enabled
# (b) Compile example programs if F77 or F90 examples are enabled
# ---------------------------------------------------------------

# Do we need a Fortran name-mangling scheme?
if(F77_INTERFACE_ENABLE OR LAPACK_ENABLE)
  set(NEED_FORTRAN_NAME_MANGLING TRUE)
endif()

# Did the user provide a name-mangling scheme?
if(SUNDIALS_F77_FUNC_CASE AND SUNDIALS_F77_FUNC_UNDERSCORES)

  string(TOUPPER ${SUNDIALS_F77_FUNC_CASE} SUNDIALS_F77_FUNC_CASE)
  string(TOUPPER ${SUNDIALS_F77_FUNC_UNDERSCORES} SUNDIALS_F77_FUNC_UNDERSCORES)

  # Based on the given case and number of underscores, set the C preprocessor
  # macro definitions. Since SUNDIALS never uses symbols names containing
  # underscores we set the name-mangling schemes to be the same. In general,
  # names of symbols with and without underscore may be mangled differently
  # (e.g. g77 mangles mysub to mysub_ and my_sub to my_sub__)
  if(SUNDIALS_F77_FUNC_CASE MATCHES "LOWER")
    if(SUNDIALS_F77_FUNC_UNDERSCORES MATCHES "NONE")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name")
    elseif(SUNDIALS_F77_FUNC_UNDERSCORES MATCHES "ONE")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name ## _")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name ## _")
    elseif(SUNDIALS_F77_FUNC_UNDERSCORES MATCHES "TWO")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name ## __")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name ## __")
    else()
      message(FATAL_ERROR "Invalid SUNDIALS_F77_FUNC_UNDERSCORES option.")
    endif()
  elseif(SUNDIALS_F77_FUNC_CASE MATCHES "UPPER")
    if(SUNDIALS_F77_FUNC_UNDERSCORES MATCHES "NONE")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME")
    elseif(SUNDIALS_F77_FUNC_UNDERSCORES MATCHES "ONE")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME ## _")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## _")
    elseif(SUNDIALS_F77_FUNC_UNDERSCORES MATCHES "TWO")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME ## __")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## __")
    else()
      message(FATAL_ERROR "Invalid SUNDIALS_F77_FUNC_UNDERSCORES option.")
    endif()
  else()
    message(FATAL_ERROR "Invalid SUNDIALS_F77_FUNC_CASE option.")
  endif()

  # name-mangling scheme has been manually set
  set(NEED_FORTRAN_NAME_MANGLING FALSE)

endif()

# Do we need a Fortran compiler?
if(F2003_INTERFACE_ENABLE OR EXAMPLES_ENABLE_F77 OR EXAMPLES_ENABLE_F90 OR NEED_FORTRAN_NAME_MANGLING)
  include(SundialsFortran)
endif()

# Ensure that F90 compiler is found if F90 examples are enabled
if (EXAMPLES_ENABLE_F90 AND (NOT F90_FOUND))
  print_warning("Compiler with F90 support not found" "Disabling F90 Examples")
  set(DOCSTR "Build F90 examples")
  force_variable(EXAMPLES_ENABLE_F90 "${DOCSTR}" OFF)
endif()

# Ensure that F90 compiler found if F2003 interface is enabled
if (F2003_INTERFACE_ENABLE AND (NOT F2003_FOUND))
  print_warning("Compiler with F2003 support not found" "Disabling F2003 Interface")
  set(DOCSTR "Enable Fortran 2003 interfaces")
  force_variable(F2003_INTERFACE_ENABLE BOOL "${DOCSTR}" OFF)
endif()


# ---------------------------------------------------------------
# A C++ compiler is needed if:
# (a) C++ examples are enabled
# (b) CUDA is enabled
# (c) RAJA is enabled
# (d) Trilinos is enabled
# (e) SuperLU_DIST is enabled
# ---------------------------------------------------------------

if(EXAMPLES_ENABLE_CXX OR CUDA_ENABLE OR RAJA_ENABLE OR Trilinos_ENABLE OR
    SUPERLUDIST_ENABLE)
  include(SundialsCXX)
endif()

# ---------------------------------------------------------------
# Setup CUDA. Since CUDA is its own language we do this
# separate from the TPLs.
# ---------------------------------------------------------------

if(CUDA_ENABLE)
  # If a user did not provide the host compiler, then we
  # assume that they want to use the CXX compiler that was set.
  if(NOT (DEFINED CMAKE_CUDA_HOST_COMPILER))
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "NVCC host compiler")
    mark_as_advanced(CMAKE_CUDA_HOST_COMPILER)
  endif()
  enable_language(CUDA)
  message(STATUS "The CUDA host compiler is: ${CMAKE_CUDA_HOST_COMPILER}")

  # Need this as long as CUDA libraries like cuSOLVER are not
  # available through some other way.
  find_package(CUDA REQUIRED)

  # Allow the user to tell us what CUDA architecture to use
  if(NOT DEFINED CUDA_ARCH)
    message(STATUS "CUDA compute architecture set to SUNDIALS default sm_30 since it was not specified")
    set(CUDA_ARCH "sm_30" CACHE STRING "Set CUDA_ARCH to SUNDIALS default sm_30" FORCE)
  endif()

  if(CUDA_ARCH)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch ${CUDA_ARCH}")
  endif()
endif(CUDA_ENABLE)

# ---------------------------------------------------------------
# Now that all languages are setup, we can configure them more.
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Decide how to compile MPI codes. We must check for MPI if
# MPI is enabled or if Trilinos is enabled because the Trilinos
# examples may need MPI without us turning on the MPI SUNDIALS
# components.
# ---------------------------------------------------------------

if(MPI_ENABLE OR Trilinos_ENABLE OR SUPERLUDIST_ENABLE)
  include(SundialsMPI)
endif()

if(MPI_ENABLE)
  if(NOT MPI_C_FOUND)
    print_warning("MPI not functional" "Parallel support will not be provided")
  else()
    set(IS_MPI_ENABLED "#ifndef SUNDIALS_MPI_ENABLED\n#define SUNDIALS_MPI_ENABLED 1\n#endif")
  endif()
endif()

# always define FMPI_COMM_F2C in sundials_fconfig.h file
if(MPIC_MPI2)
  set(F77_MPI_COMM_F2C "#define SUNDIALS_MPI_COMM_F2C 1")
  set(FMPI_COMM_F2C ".true.")
else()
  set(F77_MPI_COMM_F2C "#define SUNDIALS_MPI_COMM_F2C 0")
  set(FMPI_COMM_F2C ".false.")
endif()

# -------------------------------------------------------------
# Find OpenMP
# -------------------------------------------------------------

if(OPENMP_ENABLE OR OPENMP_DEVICE_ENABLE OR
    (SUPERLUDIST_ENABLE AND SUPERLUDIST_OpenMP))
  include(SundialsOpenMP)
endif()

# -------------------------------------------------------------
# Find PThreads
# -------------------------------------------------------------

if(PTHREAD_ENABLE)
  find_package(Threads REQUIRED)
  if(CMAKE_USE_PTHREADS_INIT)
    set(PTHREADS_FOUND TRUE)
    message(STATUS "Using Pthreads")
  else()
    set(PTHREADS_FOUND FALSE)
    print_error("Could not determine Pthreads compiler flags")
  endif()
endif()

# ===============================================================
# Find (and test) external packages
# ===============================================================

# ---------------------------------------------------------------
# Find (and test) the Lapack libraries
# ---------------------------------------------------------------

# If LAPACK is needed, first try to find the appropriate
# libraries and linker flags needed to link against them.

if(LAPACK_ENABLE)

  # find LAPACK and BLAS Libraries
  include(SundialsLapack)

  # show after include so FindLapack can locate LAPCK_LIBRARIES if necessary
  show_variable(LAPACK_LIBRARIES STRING "Lapack and Blas libraries" "${LAPACK_LIBRARIES}")

  if(LAPACK_LIBRARIES AND NOT LAPACK_FOUND)
    print_warning("LAPACK not functional"
                  "Blas/Lapack support will not be provided")
  else()
    #set sundials_config.h symbol via sundials_config.in
    set(SUNDIALS_BLAS_LAPACK TRUE)
  endif()

else()

  hide_variable(LAPACK_LIBRARIES)

endif()

# ---------------------------------------------------------------
# Find (and test) the SuperLUDIST libraries
# ---------------------------------------------------------------

if(SUPERLUDIST_ENABLE)
  include(SundialsSuperLUDIST)
else()
  set(SUPERLUDIST_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")
  set(SUNDIALS_SUPERLUDIST FALSE)
endif()

# ---------------------------------------------------------------
# Find (and test) the SUPERLUMT libraries
# ---------------------------------------------------------------

if(SUPERLUMT_ENABLE)
  include(SundialsSuperLUMT)
else()
  set(SUPERLUMT_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")
  set(SUNDIALS_SUPERLUMT FALSE)
endif()

# ---------------------------------------------------------------
# Find (and test) the KLU libraries
# ---------------------------------------------------------------

# If KLU is requested, first try to find the appropriate libraries to
# link against them.

if(KLU_ENABLE)

  show_variable(KLU_INCLUDE_DIR PATH "KLU include directory"
    "${KLU_INCLUDE_DIR}")
  show_variable(KLU_LIBRARY_DIR PATH
    "Klu library directory" "${KLU_LIBRARY_DIR}")

  include(SundialsKLU)

  if(KLU_FOUND)
    # sundials_config.h symbol
    set(SUNDIALS_KLU TRUE)
    include_directories(${KLU_INCLUDE_DIR})
  endif(KLU_FOUND)

  if(KLU_LIBRARIES AND NOT KLU_FOUND)
    print_warning("KLU not functional - support will not be provided"
                  "Double check spelling of include path and specified libraries (search is case sensitive)")
  endif(KLU_LIBRARIES AND NOT KLU_FOUND)

else()

  hide_variable(KLU_LIBRARY_DIR)
  hide_variable(KLU_INCLUDE_DIR)
  set (KLU_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")

endif(KLU_ENABLE)

# ---------------------------------------------------------------
# Find (and test) the hypre libraries
# ---------------------------------------------------------------

# >>>>>>> NOTE: Need to add check for hypre precision and integer type

if(HYPRE_ENABLE)
  show_variable(HYPRE_INCLUDE_DIR PATH "HYPRE include directory"
    "${HYPRE_INCLUDE_DIR}")
  show_variable(HYPRE_LIBRARY_DIR PATH
    "HYPRE library directory" "${HYPRE_LIBRARY_DIR}")

  include(SundialsHypre)

  if(HYPRE_FOUND)
    # sundials_config.h symbol
    set(SUNDIALS_HYPRE TRUE)
    include_directories(${HYPRE_INCLUDE_DIR})
  endif(HYPRE_FOUND)

  if(HYPRE_LIBRARIES AND NOT HYPRE_FOUND)
    print_warning("HYPRE not functional - support will not be provided"
                  "Found hypre library, test code does not work")
  endif(HYPRE_LIBRARIES AND NOT HYPRE_FOUND)

else()

  hide_variable(HYPRE_INCLUDE_DIR)
  hide_variable(HYPRE_LIBRARY_DIR)
  set (HYPRE_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")

endif()

# ---------------------------------------------------------------
# Find (and test) the PETSc libraries
# ---------------------------------------------------------------

if(PETSC_ENABLE)
  include(SundialsPETSC)
else()
  set(PETSC_DISABLED TRUE CACHE INTERNAL "GUI - return when first set")
endif()

# -------------------------------------------------------------
# Find RAJA
# -------------------------------------------------------------

if(RAJA_ENABLE)
  # Look for CMake configuration file in RAJA installation
  find_package(RAJA CONFIG
    PATHS ${RAJA_DIR} ${RAJA_DIR}/share/raja/cmake
    REQUIRED)
endif(RAJA_ENABLE)


# -------------------------------------------------------------
# Find Trilinos
# -------------------------------------------------------------

if(Trilinos_ENABLE)
  include(SundialsTrilinos)
endif(Trilinos_ENABLE)


# ===============================================================
# At this point all the configuration options are set.
# ===============================================================

# ---------------------------------------------------------------
# Configure the header file sundials_config.h
# ---------------------------------------------------------------

# All required substitution variables should be available at this point.
# Generate the header file and place it in the binary dir.
configure_file(
  ${PROJECT_SOURCE_DIR}/include/sundials/sundials_config.in
  ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  )
configure_file(
  ${PROJECT_SOURCE_DIR}/include/sundials/sundials_fconfig.in
  ${PROJECT_BINARY_DIR}/include/sundials/sundials_fconfig.h
  )

# Add the include directory in the source tree and the one in
# the binary tree (for the header file sundials_config.h)
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)

# ---------------------------------------------------------------
# Enable testing and add source and example files to the build.
# ---------------------------------------------------------------

# Enable testing
if(EXAMPLES_ENABLED)
  include(SundialsTesting)
endif()

# Add selected packages and modules to the build
add_subdirectory(src)

# Add selected examples to the build
if(EXAMPLES_ENABLED)
  add_subdirectory(examples)
endif()

# Add unit tests to the build
if(SUNDIALS_TEST_UNITTESTS)
  add_subdirectory(test/unit_tests)
endif()

# ---------------------------------------------------------------
# Install configuration header files and license file
# ---------------------------------------------------------------

# install configured header file
install(
  FILES ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  DESTINATION include/sundials
  )

# install configured header file for Fortran 90
install(
  FILES ${PROJECT_BINARY_DIR}/include/sundials/sundials_fconfig.h
  DESTINATION include/sundials
  )

# install shared Fortran 2003 modules
if(F2003_INTERFACE_ENABLE)
  # While the .mod files get generated for static and shared
  # libraries, they are identical. So only install one set
  # of the .mod files.
  if(BUILD_STATIC_LIBS)
    install(
      DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}_STATIC/
      DESTINATION ${Fortran_INSTALL_MODDIR}
      )
  else()
    install(
      DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}_SHARED/
      DESTINATION ${Fortran_INSTALL_MODDIR}
      )
  endif()
endif()

# install license and notice files
install(
  FILES ${PROJECT_SOURCE_DIR}/LICENSE
  DESTINATION include/sundials
  )
install(
  FILES ${PROJECT_SOURCE_DIR}/NOTICE
  DESTINATION include/sundials
  )
