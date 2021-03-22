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
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
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
# Fortran settings
# ===============================================================

# ---------------------------------------------------------------
# A Fortran compiler is needed to:
# (a) Determine compiler the name-mangling scheme if the F77
#     interfaces or LAPACK are enabled
# (b) Compile example programs if F77 or F90 examples are enabled
# ---------------------------------------------------------------

# Do we need a Fortran name-mangling scheme?
if(BUILD_FORTRAN77_INTERFACE OR ENABLE_LAPACK)
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
# contain underscores. This is all we really need for the FCMIX and LAPACK
# interfaces. A working Fortran compiler is only necessary for building Fortran
# example programs.
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
    EXAMPLES_ENABLE_F77 OR
    EXAMPLES_ENABLE_F90 OR
    NEED_FORTRAN_NAME_MANGLING)
  include(SundialsSetupFortran)
endif()

# ===============================================================
# C++ settings
# ===============================================================

# ---------------------------------------------------------------
# A C++ compiler is only needed if:
# (a) C++ examples are enabled
# (b) CUDA is enabled
# (c) HIP is enabled
# (d) SYCL is enabled
# (e) RAJA is enabled
# (f) Trilinos is enabled
# (g) SuperLU_DIST is enabled
# (e) MAGMA is enabled
# ---------------------------------------------------------------

if(EXAMPLES_ENABLE_CXX OR
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
  list(APPEND SUNDIALS_TPL_LIST "HIP")
endif()

# ===============================================================
# Configure presentation of language options
# ===============================================================

set(build_types DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
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

# Make build type specific flag options ADVANCED,
# except for the one corresponding to the current build type
foreach(lang ${_SUNDIALS_ENABLED_LANGS})
  foreach(build_type ${build_types})
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    if(${_cmake_build_type} MATCHES "${build_type}")
      message("Appending ${lang} ${build_type} flags")
      mark_as_advanced(CLEAR CMAKE_${lang}_FLAGS_${build_type})
    else()
      mark_as_advanced(FORCE CMAKE_${lang}_FLAGS_${build_type})
    endif()
  endforeach()
  # show the language compiler and flags
  mark_as_advanced(CLEAR CMAKE_${lang}_COMPILER CMAKE_${lang}_FLAGS)
endforeach()
