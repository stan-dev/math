# -----------------------------------------------------------------------------
# Programmer: Cody J. Balos and Slaven Peles @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# - Find Trilinos
#   This will create a Trilinos::Trilinos target which can be used in
#   target_link_libraries to get all the necessary libraries and includes
#   without using the Trilinos_LIBRARIES or Trilinos_INCLUDE_DIRS variables.
#   However, the Trilinos_CXX_COMPILER_FLAGS sill has to be used. Also creates
#   the variable Trilinos_MPI which indicates if Trilinos was built with MPI.
# -----------------------------------------------------------------------------

set(Trilinos_FOUND FALSE)

# First try and find Trilinos using Trilinos_DIR only.
find_package(Trilinos
  NAMES Trilinos TRILINOS
  PATHS
    ${Trilinos_DIR}/lib/cmake/Trilinos
    ${Trilinos_DIR}
  NO_DEFAULT_PATH)

# If Trilinos_DIR was not provided, try and find Trilinos
# somewhere else unless using in the xSDK mode.
if (NOT (Trilinos_FOUND OR USE_XSDK_DEFAULTS))
  find_package(Trilinos
    NAMES Trilinos TRILINOS
    PATHS
      ${Trilinos_DIR}/lib/cmake/Trilinos
      ${Trilinos_DIR})
endif()

if(Trilinos_FOUND)
  message(STATUS "Looking for Trilinos... success")

  # check if Trilinos was built with MPI
  if(";${Trilinos_TPL_LIST};" MATCHES ";MPI;")
    set(Trilinos_MPI TRUE)
  else()
    set(Trilinos_MPI FALSE)
  endif()
else()
  message(STATUS "Looking for Trilinos... failed")
endif()

if(Trilinos_FOUND AND NOT TARGET Trilinos::Trilinos)
  add_library(Trilinos::Trilinos IMPORTED INTERFACE)
  set_target_properties(Trilinos::Trilinos PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Trilinos_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${Trilinos_LIBRARIES}")
endif()

