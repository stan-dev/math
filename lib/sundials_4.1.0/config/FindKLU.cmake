# ---------------------------------------------------------------
# Programmer:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Find KLU library.
# 

# Set library prefixes for Windows
if(WIN32)
  set(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})
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
    FIND_LIBRARY(AMD_LIBRARY ${AMD_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(AMD_LIBRARY)
endif ()

if (NOT COLAMD_LIBRARY)
    set(COLAMD_LIBRARY_NAME colamd)
    FIND_LIBRARY(COLAMD_LIBRARY ${COLAMD_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(COLAMD_LIBRARY)
endif ()

if (NOT BTF_LIBRARY)
    set(BTF_LIBRARY_NAME btf)
    FIND_LIBRARY( BTF_LIBRARY ${BTF_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(BTF_LIBRARY)
endif ()

if (NOT SUITESPARSECONFIG_LIBRARY)
    set(SUITESPARSECONFIG_LIBRARY_NAME suitesparseconfig)
    # NOTE: no prefix for this library on windows
    if(WIN32 AND NOT MSYS)
        set(CMAKE_FIND_LIBRARY_PREFIXES "")
    endif()
    FIND_LIBRARY( SUITESPARSECONFIG_LIBRARY ${SUITESPARSECONFIG_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(SUITESPARSECONFIG_LIBRARY)
endif ()

set(KLU_LIBRARIES ${KLU_LIBRARY} ${AMD_LIBRARY} ${COLAMD_LIBRARY} ${BTF_LIBRARY} ${SUITESPARSECONFIG_LIBRARY})
