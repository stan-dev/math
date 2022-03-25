# ---------------------------------------------------------------
# Programmer(s): Steven Smith and Cody J. Balos @ LLNL
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
# KLU find module that creates an imported target for KLU.
# The target is SUNDIALS::KLU.
#
# The variable KLU_LIBRARY_DIR can be used to control
# where the module looks for the library.
#
# The variable KLU_INCLUDE_DIR can be used to set the
# include path for the library.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   KLU_FOUND       - system has KLU library
#   KLU_LIBRARY     - the KLU library
#   KLU_INCLUDE_DIR - the KLU include path
#   KLU_LIBRARIES   - all of the libraries needed for KLU
# ---------------------------------------------------------------

# Set library prefixes for Windows
if(WIN32)
  set(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})
  set(CMAKE_FIND_LIBRARY_SUFFIXES d.lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
elseif(APPLE)
  set(CMAKE_FIND_LIBRARY_SUFFIXES d.a ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

### Find include dir
find_path(temp_KLU_INCLUDE_DIR klu.h ${KLU_INCLUDE_DIR})
if (temp_KLU_INCLUDE_DIR)
    set(KLU_INCLUDE_DIR ${temp_KLU_INCLUDE_DIR})
endif()
unset(temp_KLU_INCLUDE_DIR CACHE)

if (KLU_LIBRARY)
    # We have (or were given) KLU_LIBRARY - get path to use for other Suitesparse libs
    get_filename_component(KLU_LIBRARY_DIR ${KLU_LIBRARY} PATH)

    # force CACHE update to show user DIR that will be used
    set(KLU_LIBRARY_DIR ${KLU_LIBRARY_DIR} CACHE PATH "" FORCE)

else ()
    # find library with user provided directory path
    set(KLU_LIBRARY_NAME klu)
    find_library(KLU_LIBRARY ${KLU_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
endif ()
mark_as_advanced(KLU_LIBRARY)

if (NOT AMD_LIBRARY)
    set(AMD_LIBRARY_NAME amd)
    find_library(AMD_LIBRARY ${AMD_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(AMD_LIBRARY)
endif ()

if (NOT COLAMD_LIBRARY)
    set(COLAMD_LIBRARY_NAME colamd)
    find_library(COLAMD_LIBRARY ${COLAMD_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(COLAMD_LIBRARY)
endif ()

if (NOT BTF_LIBRARY)
    set(BTF_LIBRARY_NAME btf)
    find_library( BTF_LIBRARY ${BTF_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(BTF_LIBRARY)
endif ()

if (NOT SUITESPARSECONFIG_LIBRARY)
    set(SUITESPARSECONFIG_LIBRARY_NAME suitesparseconfig)
    # NOTE: no prefix for this library on windows
    if(WIN32 AND NOT MSYS)
        set(CMAKE_FIND_LIBRARY_PREFIXES "")
    endif()
    find_library( SUITESPARSECONFIG_LIBRARY ${SUITESPARSECONFIG_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(SUITESPARSECONFIG_LIBRARY)
endif ()

set(KLU_LIBRARIES ${KLU_LIBRARY} ${AMD_LIBRARY} ${COLAMD_LIBRARY} ${BTF_LIBRARY} ${SUITESPARSECONFIG_LIBRARY})

# set package variables including KLU_FOUND
find_package_handle_standard_args(KLU
  REQUIRED_VARS
    KLU_LIBRARY
    KLU_LIBRARIES
    KLU_INCLUDE_DIR)

# Create target for KLU
if(KLU_FOUND)

  if(NOT TARGET SUNDIALS::KLU)
    add_library(SUNDIALS::KLU UNKNOWN IMPORTED)
  endif()

  set_target_properties(SUNDIALS::KLU PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${KLU_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${KLU_LIBRARIES}"
    IMPORTED_LOCATION "${KLU_LIBRARY}")

endif()
