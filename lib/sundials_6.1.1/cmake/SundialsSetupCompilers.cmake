# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# Module that sets up compilers for SUNDIALS.
# ---------------------------------------------------------------

# ===============================================================
# Determine the index type for the compiler
# ===============================================================

include(SundialsIndexSize)

# ===============================================================
# Platform specifc settings
# ===============================================================

if(WIN32)
  # Under Windows, add compiler directive to inhibit warnings
  # about use of unsecure functions.
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

if(APPLE)
  # Allow undefined symbols that will be resolved by a user program.
  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
endif()

# ===============================================================
# RPath settings
# ===============================================================

# only apply rpath settings for builds using shared libs
if(BUILD_SHARED_LIBS)
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  endif()
endif()

# ===============================================================
# C settings
# ===============================================================

set(DOCSTR "The C standard to use (90, 99, 11, 17)")
sundials_option(CMAKE_C_STANDARD STRING "${DOCSTR}" "99"
                OPTIONS "90;99;11;17")
message(STATUS "C standard set to ${CMAKE_C_STANDARD}")

set(DOCSTR "Enable C compiler specific extensions")
sundials_option(CMAKE_C_EXTENSIONS BOOL "${DOCSTR}" OFF)
message(STATUS "C extensions set to ${CMAKE_C_EXTENSIONS}")

# Profiling generally requires ISO C99 or newer for __func__ though some
# compilers define __func__ even with ISO C90.
if(SUNDIALS_BUILD_WITH_PROFILING AND (CMAKE_C_STANDARD STREQUAL "90"))
  message(WARNING "SUNDIALS_BUILD_WITH_PROFILING=ON requires __func__, compilation may fail with CMAKE_C_STANDARD=90")
endif()

# ---------------------------------------------------------------
# Check for POSIX timers
#
# 199309L is the minimum POSIX version needed for struct timespec
# and clock_monotonic()
# ---------------------------------------------------------------
include(SundialsPOSIXTimers)

if(SUNDIALS_POSIX_TIMERS AND POSIX_TIMERS_NEED_POSIX_C_SOURCE)
  set(DOCSTR "Value of _POSIX_C_SOURCE")
  sundials_option(SUNDIALS_POSIX_C_SOURCE STRING "${DOCSTR}" "199309L"
                  ADVANCED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=${SUNDIALS_POSIX_C_SOURCE}")
endif()

# Check if profiling is being built with no timers.
if(SUNDIALS_BUILD_WITH_PROFILING AND
   (NOT ENABLE_CALIPER) AND
   (NOT ENABLE_MPI) AND
   (NOT SUNDIALS_POSIX_TIMERS))
  message(SEND_ERROR "The SUNDIALS native profiler requires POSIX timers or MPI_Wtime, but neither were found.")
endif()

# ---------------------------------------------------------------
# Check for deprecated attribute with message
# ---------------------------------------------------------------
if(WIN32)
  set(COMPILER_DEPRECATED_MSG_ATTRIBUTE "__declspec(deprecated(msg))" CACHE INTERNAL "")
else()
  set(COMPILER_DEPRECATED_MSG_ATTRIBUTE "__attribute__ ((__deprecated__(msg)))" CACHE INTERNAL "")
endif()
check_c_source_compiles("
  #define msg \"test\"
  ${COMPILER_DEPRECATED_MSG_ATTRIBUTE} int somefunc() { return 0; }
  int main() { return somefunc();}" COMPILER_HAS_DEPRECATED_MSG
)

# ===============================================================
# Fortran settings
# ===============================================================

# When LAPACK is enabled we will need a Fortran compiler to infer the
# name-mangling scheme if it is not set by the user
if(ENABLE_LAPACK)
  set(NEED_FORTRAN_NAME_MANGLING TRUE)
endif()

# ------------------------------------------------------------------------------
# Allow the user to manually specify the Fortran name-mangling scheme
#
# The build system tries to infer the Fortran name-mangling scheme using a
# Fortran compiler and defaults to using lower case and one underscore if the
# scheme can not be determined. If a working Fortran compiler is not available
# or the user needs to override the inferred or default scheme, the following
# options specify the case and number of appended underscores corresponding to
# the Fortran name-mangling scheme of symbol names that do not themselves
# contain underscores. This is all we really need for the LAPACK interfaces. A
# working Fortran compiler is only necessary for building Fortran example
# programs.
# ------------------------------------------------------------------------------

# The case to use in the name-mangling scheme
sundials_option(SUNDIALS_F77_FUNC_CASE STRING
                "case of Fortran function names (lower/upper)"
                ""
                ADVANCED)

# The number of underscores of appended in the name-mangling scheme
sundials_option(SUNDIALS_F77_FUNC_UNDERSCORES STRING
                "number of underscores appended to Fortran function names (none/one/two)"
                ""
                ADVANCED)

# If used, both case and underscores must be set
if((NOT SUNDIALS_F77_FUNC_CASE) AND SUNDIALS_F77_FUNC_UNDERSCORES)
  print_error("If SUNDIALS_F77_FUNC_UNDERSCORES is set, "
                      "SUNDIALS_F77_FUNC_CASE must also be set.")
endif()
if(SUNDIALS_F77_FUNC_CASE AND (NOT SUNDIALS_F77_FUNC_UNDERSCORES))
  print_error("If SUNDIALS_F77_FUNC_CASE is set, "
                      "SUNDIALS_F77_FUNC_UNDERSCORES must also be set.")
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
      print_error("Invalid SUNDIALS_F77_FUNC_UNDERSCORES option.")
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
      print_error("Invalid SUNDIALS_F77_FUNC_UNDERSCORES option.")
    endif()
  else()
    print_error("Invalid SUNDIALS_F77_FUNC_CASE option.")
  endif()

  # name-mangling scheme has been manually set
  set(NEED_FORTRAN_NAME_MANGLING FALSE)

endif()

# Do we need a Fortran compiler?
if(BUILD_FORTRAN_MODULE_INTERFACE OR
    NEED_FORTRAN_NAME_MANGLING)
  include(SundialsSetupFortran)
endif()

# ===============================================================
# C++ settings
# ===============================================================

if(BUILD_BENCHMARKS OR EXAMPLES_ENABLE_CXX OR
    ENABLE_CUDA OR
    ENABLE_HIP OR
    ENABLE_SYCL OR
    ENABLE_RAJA OR
    ENABLE_TRILINOS OR
    ENABLE_SUPERLUDIST OR
    ENABLE_MAGMA)
  include(SundialsSetupCXX)
endif()

# ===============================================================
# CUDA settings
# ===============================================================

if(ENABLE_CUDA)
  include(SundialsSetupCuda)
  # we treat CUDA as both a TPL and a language
  list(APPEND SUNDIALS_TPL_LIST "CUDA")
endif()

# ===============================================================
# HIP settings
# ===============================================================

if(ENABLE_HIP)
  include(SundialsSetupHIP)
  # we treat HIP as both a TPL and a language
  list(APPEND SUNDIALS_TPL_LIST "HIP")
endif()

# ===============================================================
# Configure presentation of language options
# ===============================================================

# List of enabled languages
set(_SUNDIALS_ENABLED_LANGS "C")
if(CXX_FOUND)
  list(APPEND _SUNDIALS_ENABLED_LANGS "CXX")
endif()
if(Fortran_FOUND)
  list(APPEND _SUNDIALS_ENABLED_LANGS "Fortran")
endif()
if(CUDA_FOUND)
  list(APPEND _SUNDIALS_ENABLED_LANGS "CUDA")
endif()

# Upper case version of build type
string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)

# Make build type specific flag options ADVANCED,
# except for the one corresponding to the current build type
foreach(lang ${_SUNDIALS_ENABLED_LANGS})
  foreach(build_type DEBUG;RELEASE;RELWITHDEBINFO;MINSIZEREL)
    if("${_cmake_build_type}" STREQUAL "${build_type}")
      message(STATUS "Appending ${lang} ${build_type} flags")
      mark_as_advanced(CLEAR CMAKE_${lang}_FLAGS_${build_type})
    else()
      mark_as_advanced(FORCE CMAKE_${lang}_FLAGS_${build_type})
    endif()
  endforeach()
  # show the language compiler and flags
  mark_as_advanced(CLEAR CMAKE_${lang}_COMPILER CMAKE_${lang}_FLAGS)
endforeach()

# ===============================================================
# Additional compiler flags
#
# TODO(DJG): Set flags based on CMAKE_<language>_COMPILER_ID
# ===============================================================

if(ENABLE_ALL_WARNINGS)
  message(STATUS "Enabling all compiler warnings")

  set(CMAKE_C_FLAGS "-Wall -Wpedantic -Wextra -Wno-unused-parameter -Wno-deprecated-declarations -Wno-unused-function ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -Wno-unused-parameter -Wno-deprecated-declarations -Wno-unused-function ${CMAKE_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "-Wall -Wpedantic -Wno-unused-dummy-argument -Wno-c-binding-type -ffpe-summary=none ${CMAKE_Fortran_FLAGS}")
endif()

if(ENABLE_WARNINGS_AS_ERRORS)
  message(STATUS "Enabling compiler warnings as errors")

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror")
endif()

if(ENABLE_ADDRESS_SANITIZER)
  message(STATUS "Enabling address sanitizer")

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
endif()

if(SUNDIALS_DEBUG)
  message(STATUS "Adding debugging preprocessor directives")

  foreach(debug ${_SUNDIALS_DEBUG_OPTIONS})
    if (${${debug}})
      add_compile_definitions(${debug})
    endif()
  endforeach()
endif()

# ===============================================================
# Configure compilers for installed examples
# ===============================================================

foreach(lang ${_SUNDIALS_ENABLED_LANGS})
  if(SUNDIALS_BUILD_WITH_PROFILING AND ENABLE_MPI)
    if(DEFINED MPI_${lang}_COMPILER)
      set(_EXAMPLES_${lang}_COMPILER "${MPI_${lang}_COMPILER}" CACHE INTERNAL "${lang} compiler for installed examples")
    endif()
  else()
    set(_EXAMPLES_${lang}_COMPILER "${CMAKE_${lang}_COMPILER}" CACHE INTERNAL "${lang} compiler for installed examples")
  endif()
endforeach()
