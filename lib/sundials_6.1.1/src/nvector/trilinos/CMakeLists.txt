# ---------------------------------------------------------------
# Programmer(s): Slaven Peles, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the Trilinos NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_Trilinos\n\")")

# Set Trilinos compilers/flags
set(CMAKE_CXX_COMPILER ${Trilinos_INTERFACE_CXX_COMPILER})
set(CMAKE_C_COMPILER   ${Trilinos_INTERFACE_C_COMPILER})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_INTERFACE_CXX_COMPILER_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Trilinos_INTERFACE_C_COMPILER_FLAGS}")

# Create the library
sundials_add_library(sundials_nvectrilinos
  SOURCES
    nvector_trilinos.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_trilinos.h
  INCLUDE_SUBDIR
    nvector
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

# Install the exported NVECTRILINOS CXX header files in a separate subdirectory
install(FILES
        ${SUNDIALS_SOURCE_DIR}/include/nvector/trilinos/SundialsTpetraVectorInterface.hpp
        ${SUNDIALS_SOURCE_DIR}/include/nvector/trilinos/SundialsTpetraVectorKernels.hpp
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nvector/trilinos")

message(STATUS "Added NVECTOR_Trilinos module")
