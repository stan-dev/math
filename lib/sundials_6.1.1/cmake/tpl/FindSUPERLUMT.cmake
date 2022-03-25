# ---------------------------------------------------------------
# Programmer(s): Eddy Banks and David J. Gardner @ LLNL
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
# SuperLUMT find module that creates an imported target for
# SuperLU_MT. The target is SUNDIALS::SUPERLUMT.
#
# The variable SUPERLUMT_LIBRARY_DIR can be used to control
# where the module looks for the library.
#
# The variable SUPERLUMT_INCLUDE_DIR can be used to set the
# include path for the library.
#
# Additional libraries can be passed in SUPERLUMT_LIBRARIES.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   SUPERLUMT_FOUND       - system has SuperLU_MT library
#   SUPERLUMT_LIBRARY     - the SuperLU_MT library
#   SUPERLUMT_THREAD_TYPE - the SuperLU_MT threading
# ---------------------------------------------------------------

# check for valid thread type
string(TOUPPER ${SUPERLUMT_THREAD_TYPE} _upper_SUPERLUMT_THREAD_TYPE)
force_variable(SUPERLUMT_THREAD_TYPE STRING "SuperLU_MT threading type: OPENMP or PTHREAD" ${_upper_SUPERLUMT_THREAD_TYPE})

if(SUPERLUMT_THREAD_TYPE AND
    NOT SUPERLUMT_THREAD_TYPE STREQUAL "OPENMP" AND
    NOT SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
  print_error("Unknown thread type: ${SUPERLUMT_THREAD_TYPE}" "Please enter PTHREAD or OPENMP")
endif()

# check if the threading library has been found
if(SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
  # find pthread libraries
  if(NOT PTHREADS_FOUND)
    find_package(Threads REQUIRED)
    if(CMAKE_USE_PTHREADS_INIT)
      set(PTHREADS_FOUND TRUE)
      message(STATUS "Using Pthreads")
    else()
      set(PTHREADS_FOUND FALSE)
      print_error("Could not determine Pthreads compiler flags")
    endif()
  endif()
else(SUPERLUMT_THREAD_TYPE STREQUAL "OPENMP")
  # find openmp libraries
  if(NOT OPENMP_FOUND)
    find_package(OpenMP REQUIRED)
  endif()
endif()

# Set SuperLU_MT library name with thread type postfix
set(SUPERLUMT_LIBRARY_NAME superlu_mt_${SUPERLUMT_THREAD_TYPE})

if(MSVC)
  set(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})
endif()

# Check if SUPERLUMT_LIBRARIES contains the superlu_mt
# library as well as TPLs. If so, extract it into the
# SUPERLUMT_LIBRARY variable.
if(SUPERLUMT_LIBRARIES MATCHES "${SUPERLUMT_LIBRARY_NAME}")
  foreach(lib ${SUPERLUMT_LIBRARIES})
    if(lib MATCHES "${SUPERLUMT_LIBRARY_NAME}")
      set(SUPERLUMT_LIBRARY ${lib})
    endif()
  endforeach()
endif()

# find library
if(NOT SUPERLUMT_LIBRARY)
  # search user provided directory path
  find_library(SUPERLUMT_LIBRARY ${SUPERLUMT_LIBRARY_NAME}
    PATHS ${SUPERLUMT_LIBRARY_DIR} NO_DEFAULT_PATH)
  # if user didn't provide a path, search anywhere
  if(NOT (SUPERLUMT_LIBRARY_DIR OR SUPERLUMT_LIBRARY))
    find_library(SUPERLUMT_LIBRARY ${SUPERLUMT_LIBRARY_NAME})
  endif()
  mark_as_advanced(SUPERLUMT_LIBRARY)
endif()

# set the libraries, stripping out 'NOTFOUND' from previous attempts
string(REPLACE "SUPERLUMT_LIBRARY-NOTFOUND" "" SUPERLUMT_LIBRARIES "${SUPERLUMT_LIBRARIES}")
set(SUPERLUMT_LIBRARIES "${SUPERLUMT_LIBRARY};${SUPERLUMT_LIBRARIES}" CACHE STRING "" FORCE)

# set the library dir option if it wasn't preset
if(SUPERLUMT_LIBRARY AND (NOT SUPERLUMT_LIBRARY_DIR))
  get_filename_component(SUPERLUMT_LIBRARY_DIR ${SUPERLUMT_LIBRARY} DIRECTORY)
  set(SUPERLUMT_LIBRARY_DIR ${SUPERLUMT_LIBRARY_DIR} CACHE PATH "" FORCE)
endif()

# set the include dir option if it wasn't preset
if(SUPERLUMT_LIBRARY AND (NOT SUPERLUMT_INCLUDE_DIR))
  get_filename_component(SUPERLUMT_INCLUDE_DIR ${SUPERLUMT_LIBRARY_DIR} DIRECTORY)
  set(SUPERLUMT_INCLUDE_DIR "${SUPERLUMT_INCLUDE_DIR}/include" CACHE PATH "" FORCE)
endif()

# set a more informative error message in case the library was not found
set(SUPERLUMT_NOT_FOUND_MESSAGE "\
************************************************************************\n\
ERROR: Could not find SuperLU_MT. Please check the variables:\n\
       SUPERLUMT_INCLUDE_DIR and SUPERLUMT_LIBRARY_DIR\n\
************************************************************************")

# set package variables including SUPERLUMT_FOUND
find_package_handle_standard_args(SUPERLUMT
  REQUIRED_VARS
    SUPERLUMT_LIBRARY
    SUPERLUMT_LIBRARIES
    SUPERLUMT_INCLUDE_DIR
    SUPERLUMT_THREAD_TYPE
  FAIL_MESSAGE
    "${SUPERLUMT_NOT_FOUND_MESSAGE}"
  )

# Create target for SuperLU_MT
if(SUPERLUMT_FOUND)

  if(NOT TARGET SUNDIALS::SUPERLUMT)
    add_library(SUNDIALS::SUPERLUMT UNKNOWN IMPORTED)
  endif()

  set_target_properties(SUNDIALS::SUPERLUMT PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SUPERLUMT_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${SUPERLUMT_LIBRARIES}"
    IMPORTED_LOCATION "${SUPERLUMT_LIBRARY}")

endif()
