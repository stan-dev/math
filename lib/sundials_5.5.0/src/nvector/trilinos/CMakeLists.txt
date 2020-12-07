# ---------------------------------------------------------------
# Programmer(s): Slaven Peles, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the Trilinos NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_Trilinos\n\")")

# Set Trilinos compilers/flags
set(CMAKE_CXX_COMPILER ${Trilinos_INTERFACE_CXX_COMPILER})
set(CMAKE_C_COMPILER   ${Trilinos_INTERFACE_C_COMPILER})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_INTERFACE_CXX_COMPILER_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Trilinos_INTERFACE_C_COMPILER_FLAGS}")

# Add variable nvectrilinos_HEADERS with the exported NVECTRILINOS header files
set(nvectrilinos_HEADERS
  ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_trilinos.h
  )

set(nvectrilinos_CXXHEADERS
  ${SUNDIALS_SOURCE_DIR}/include/nvector/trilinos/SundialsTpetraVectorInterface.hpp
  ${SUNDIALS_SOURCE_DIR}/include/nvector/trilinos/SundialsTpetraVectorKernels.hpp
  )

# Create the library
sundials_add_library(sundials_nvectrilinos
  SOURCES
    nvector_trilinos.cpp
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC SUNDIALS::TRILINOS
  OUTPUT_NAME
    sundials_nvectrilinos
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

# Install the NVECTRILINOS header files
install(FILES ${nvectrilinos_HEADERS} DESTINATION include/nvector)
install(FILES ${nvectrilinos_CXXHEADERS} DESTINATION include/nvector/trilinos)

#
message(STATUS "Added NVECTOR_Trilinos module")
