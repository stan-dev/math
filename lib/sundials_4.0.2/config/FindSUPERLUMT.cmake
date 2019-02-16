# ---------------------------------------------------------------
# Programmer:  Eddy Banks @ LLNL
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
# SUPERLUMT tests for SUNDIALS CMake-based configuration.
# 

# make sure valid thread type - if not, then warn and return
STRING(TOUPPER "${SUPERLUMT_THREAD_TYPE}" SUPERLUMT_THREAD_TYPE_UPPER)
If(SUPERLUMT_THREAD_TYPE AND
    NOT SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "OPENMP" AND
    NOT SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
    PRINT_WARNING("Unknown thread type: ${SUPERLUMT_THREAD_TYPE}" "Please enter Pthread or OpenMP")
ENDIF()

# find the SUPERLUMT include directory path
IF(SUPERLUMT_THREAD_TYPE)
  # if have user input for thread type - set postfix of library name
  set(POST ${SUPERLUMT_THREAD_TYPE_UPPER})

  ### Find include dir
  find_path(temp_SUPERLUMT_INCLUDE_DIR slu_mt_ddefs.h ${SUPERLUMT_INCLUDE_DIR})
  if (temp_SUPERLUMT_INCLUDE_DIR)
    set(SUPERLUMT_INCLUDE_DIR ${temp_SUPERLUMT_INCLUDE_DIR})
  endif()
  unset(temp_SUPERLUMT_INCLUDE_DIR CACHE)
ENDIF()

IF(MSVC)
  SET(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})
ENDIF()

if(SUPERLUMT_LIBRARY)
  get_filename_component(SUPERLUMT_LIBRARY_DIR ${SUPERLUMT_LIBRARY} PATH)
  set(SUPERLUMT_LIBRARY_DIR ${SUPERLUMT_LIBRARY_DIR} CACHE PATH "" FORCE)
else()
  # find library with user provided directory path
  set(SUPERLUMT_LIBRARY_NAME superlu_mt_${POST})
  find_library(SUPERLUMT_LIBRARY ${SUPERLUMT_LIBRARY_NAME} ${SUPERLUMT_LIBRARY_DIR} NO_DEFAULT_PATH)
endif()
mark_as_advanced(SUPERLUMT_LIBRARY)

# add threading library (pthread or openmp)
If(SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
  # add pthread to libraries
  find_library(SUPERLUMT_THREAD_LIBRARY
      NAMES pthread
      PATHS /usr/lib /usr/local/lib
      "$ENV{ProgramFiles}/SUPERLUMT/Lib"
      )
ELSE(SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "OPENMP")
  # add openmp to libraries
  if(NOT OPENMP_FOUND)
    find_package( OpenMP REQUIRED)
  endif()
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()  
ENDIF()
mark_as_advanced(SUPERLUMT_THREAD_LIBRARY)

# add to SUPERLUMT_LIBRARIES (Note: will be 'not found' if either are not found 
set(SUPERLUMT_LIBRARIES ${SUPERLUMT_LIBRARY} ${SUPERLUMT_THREAD_LIBRARY})

# If LAPACK/BLAS not enabled - find BLAS with SUPERLUMT
if(NOT BLAS_ENABLE AND NOT LAPACK_ENABLE)
  set(SUPERLUMT_BLAS_LIBRARY_NAME blas_${POST})

  #unset(SUPERLUMT_BLAS_LIBRARIES CACHE)
  FIND_LIBRARY(SUPERLUMT_BLAS_LIBRARIES ${SUPERLUMT_BLAS_LIBRARY_NAME} ${SUPERLUMT_LIBRARY_DIR} NO_DEFAULT_PATH)

  if (NOT SUPERLUMT_BLAS_LIBRARIES)
    PRINT_WARNING("Can't find SUPERLUMT_BLAS_LIBRARY, support will not be provided."
                  "Try setting BLAS_ENABLE or LAPACK_ENABLE to ON")
    SET(SUPERLUMT_FOUND FALSE)
  else ()
    set(SUPERLUMT_LIBRARIES ${SUPERLUMT_LIBRARIES} ${SUPERLUMT_BLAS_LIBRARIES})
  endif ()
  mark_as_advanced(SUPERLUMT_BLAS_LIBRARIES)
endif(NOT BLAS_ENABLE AND NOT LAPACK_ENABLE)
