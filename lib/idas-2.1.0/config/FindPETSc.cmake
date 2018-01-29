# ---------------------------------------------------------------
# Programmer:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# Find PETSC library.
# 

IF(WIN32)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")
endif(WIN32)

### Find include dir
find_path(temp_PETSC_INCLUDE_DIR petsc.h ${PETSC_INCLUDE_DIR})
if (temp_PETSC_INCLUDE_DIR)
    set(PETSC_INCLUDE_DIR ${temp_PETSC_INCLUDE_DIR})
endif()
unset(temp_PETSC_INCLUDE_DIR CACHE)

if (PETSC_LIBRARY)
    # We have (or were given) PETSC_LIBRARY - get path to use for any related libs
    get_filename_component(PETSC_LIBRARY_DIR ${PETSC_LIBRARY} PATH)
    
    # force CACHE update to show user DIR that will be used
    set(PETSC_LIBRARY_DIR ${PETSC_LIBRARY_DIR} CACHE PATH "" FORCE)
    
else ()
    # find library with user provided directory path
    set(PETSC_LIBRARY_NAMES petsc PETSC)
    find_library(PETSC_LIBRARY 
      NAMES ${PETSC_LIBRARY_NAMES}
      PATHS ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH
      )
endif ()
mark_as_advanced(PETSC_LIBRARY)

set(PETSC_LIBRARIES ${PETSC_LIBRARY})
