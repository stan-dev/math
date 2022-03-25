# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# XBRAID find module that creates an imported target for XBRAID.
# The target is SUNDIALS::XBRAID.
#
# The variable XBRAID_DIR can be used to control where the module
# looks for the library.
#
#   XBRAID_LIBRARIES - (advanced) the libraries to link against
#   XBRAID_INCLUDES  - (advanced) the directories to include
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   XBRAID_FOUND - system has the XBRAID library
# ------------------------------------------------------------------------------

# Check if we are locating XBraid using the root install directory or a list of
# include directories and link libraries
if (XBRAID_INCLUDES OR XBRAID_LIBRARIES)

  if (XBRAID_INCLUDES AND XBRAID_LIBRARIES)

    set(XBRAID_DIR "" CACHE PATH "Path to the root of XBraid installation" FORCE)

  else ()

    string(CONCAT msg
      "Both XBRAID_INCLUDES and XBRAID_LIBRARIES must be provided:\n"
      "  XBRAID_INCLUDES=${XBRAID_INCLUDES}\n"
      "  XBRAID_LIBRARIES=${XBRAID_LIBRARIES}")
    message(FATAL_ERROR ${msg})

  endif ()

else ()

  set(XBRAID_INCLUDES "" CACHE STRING "Semi-colon separated list of XBraid include directories" FORCE)
  set(XBRAID_LIBRARIES "" CACHE STRING "Semi-colon separated list of XBraid link libraries" FORCE)

endif ()

# unset cache values for multiple passes
unset(XBRAID_INCLUDE_DIR CACHE)
unset(XBRAID_LIBRARY CACHE)

unset(XBRAID_INCS CACHE)
unset(XBRAID_LIBS CACHE)

if (XBRAID_INCLUDES AND XBRAID_LIBRARIES)

  message(STATUS "Finding XBraid using XBRAID_INCLUDES and XBRAID_LIBRARIES")

  # extract path from XBRAID_INCLUDES
  foreach (include_dir ${XBRAID_INCLUDES})
    if (EXISTS "${include_dir}/braid.h")
      set(XBRAID_INCLUDE_DIR "${include_dir}" CACHE "XBraid include directory")
      break()
    endif ()
  endforeach ()

  # check if the include directory was found
  if (NOT XBRAID_INCLUDE_DIR)
    string(CONCAT msg
      "Could not determine XBraid include directory from XBRAID_INCLUDES:\n"
      "  XBRAID_INCLUDES=${XBRAID_INCLUDES}\n")
    message(FATAL_ERROR ${msg})
  endif ()

  # extract library from XBRAID_LIBRARIES
  foreach (library_path ${XBRAID_LIBRARIES})
    get_filename_component(library_name "${library_path}" NAME)
    if (library_name MATCHES "braid")
      set(XBRAID_LIBRARY "${library_path}" CACHE "XBraid library")
      break()
    endif ()
  endforeach ()

  # check if the library directory was found
  if (NOT XBRAID_LIBRARY)
    string(CONCAT msg
      "Could not determine XBraid library from XBRAID_LIBRARIES:\n"
      "  XBRAID_LIBRARIES=${XBRAID_LIBRARIES}")
    message(FATAL_ERROR ${msg})
  endif ()

else ()

  message(STATUS "Finding XBraid using XBRAID_DIR")

  # find XBRAID_DIR
  if (NOT XBRAID_DIR)

    message(STATUS "Looking for XBraid in common install locations")
    find_path(XBRAID_DIR include/braid.h braid/braid.h)

  endif ()

  # check if XBRAID_DIR was set/found
  if (NOT XBRAID_DIR)

    string(CONCAT msg
      "Could not locate XBraid install directory please set:\n"
      "  - XBRAID_DIR\n"
      "or used the advanced options\n"
      "  - XBRAID_INCLUDES and XBRAID_LIBRARIES.")
    message(FATAL_ERROR ${msg})

  endif ()

  # Find the include dir
  find_path(XBRAID_INCLUDE_DIR braid.h
    PATHS
    ${XBRAID_DIR}
    PATH_SUFFIXES
    include braid
    DOC
    "XBraid include directory"
    NO_DEFAULT_PATH)

  # check if the include directory was found
  if (NOT XBRAID_INCLUDE_DIR)
    string(CONCAT msg
      "Could not determine XBraid include directory from XBRAID_DIR:\n"
      "  XBRAID_DIR=${XBRAID_DIR}\n")
    message(FATAL_ERROR ${msg})
  endif ()

  # Find the library
  find_library(XBRAID_LIBRARY braid
    PATHS
    ${XBRAID_DIR}
    PATH_SUFFIXES
    lib braid
    DOC
    "XBraid library"
    NO_DEFAULT_PATH)

  # check if the library was found
  if (NOT XBRAID_LIBRARY)
    string(CONCAT msg
      "Could not determine XBraid library from XBRAID_DIR:\n"
      "  XBRAID_DIR=${XBRAID_DIR}\n")
    message(FATAL_ERROR ${msg})
  endif ()

endif ()

# set package variables including XBRAID_FOUND
find_package_handle_standard_args(XBRAID
  REQUIRED_VARS
  XBRAID_INCLUDE_DIR
  XBRAID_LIBRARY
  )

# XBraid target
if (XBRAID_FOUND)

  # create target if necessary
  if (NOT TARGET SUNDIALS::XBRAID)
    add_library(SUNDIALS::XBRAID UNKNOWN IMPORTED)
  endif ()

  # update target properties (for multiple passes)
  set_target_properties(SUNDIALS::XBRAID PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${XBRAID_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${XBRAID_LIBRARIES}"
    IMPORTED_LOCATION "${XBRAID_LIBRARY}")

  # set variables for output message, compile tests, and
  # CMake/Makefile templates
  if (XBRAID_INCLUDES AND XBRAID_LIBRARIES)
    set(XBRAID_INCS "${XBRAID_INCLUDES}" CACHE INTERNAL
      "Internal XBraid includes")
    set(XBRAID_LIBS "${XBRAID_LIBRARIES}" CACHE INTERNAL
      "Internal XBraid libraries")
  else ()
    set(XBRAID_INCS "${XBRAID_INCLUDE_DIR}" CACHE INTERNAL
      "Internal XBraid includes")
    set(XBRAID_LIBS "${XBRAID_LIBRARY}" CACHE INTERNAL
      "Internal XBraid libraries")
  endif ()

endif ()
