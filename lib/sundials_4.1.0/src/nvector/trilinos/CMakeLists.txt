# ---------------------------------------------------------------
# Programmer:  Slaven Peles @ LLNL
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
# CMakeLists.txt file for the Trilinos NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_Trilinos\n\")")

# Set Trilinos compilers/flags
SET(CMAKE_CXX_COMPILER ${Trilinos_INTERFACE_CXX_COMPILER})
SET(CMAKE_C_COMPILER   ${Trilinos_INTERFACE_C_COMPILER})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_INTERFACE_CXX_COMPILER_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Trilinos_INTERFACE_C_COMPILER_FLAGS}")

# Add variable nvectrilinos_SOURCES with the sources for the NVECTRILINOS lib
SET(nvectrilinos_SOURCES nvector_trilinos.cpp)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECTRILINOS library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvectrilinos_HEADERS with the exported NVECTRILINOS header files
SET(nvectrilinos_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_trilinos.h
  )

SET(nvectrilinos_CXXHEADERS
  ${sundials_SOURCE_DIR}/include/nvector/trilinos/SundialsTpetraVectorInterface.hpp
  ${sundials_SOURCE_DIR}/include/nvector/trilinos/SundialsTpetraVectorKernels.hpp
  )

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECTRILINOS library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECTRILINOS library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvectrilinos_static STATIC ${nvectrilinos_SOURCES} ${shared_SOURCES})

  TARGET_LINK_LIBRARIES(sundials_nvectrilinos_static
    PUBLIC Trilinos::Trilinos)

  SET_TARGET_PROPERTIES(sundials_nvectrilinos_static
    PROPERTIES OUTPUT_NAME sundials_nvectrilinos CLEAN_DIRECT_OUTPUT 1)

  INSTALL(TARGETS sundials_nvectrilinos_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECTRILINOS library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECTRILINOS library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvectrilinos_shared SHARED ${nvectrilinos_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvectrilinos_shared
      PUBLIC m)
  ENDIF()

  TARGET_LINK_LIBRARIES(sundials_nvectrilinos_shared
    PUBLIC Trilinos::Trilinos)

  SET_TARGET_PROPERTIES(sundials_nvectrilinos_shared
    PROPERTIES OUTPUT_NAME sundials_nvectrilinos CLEAN_DIRECT_OUTPUT 1)

  SET_TARGET_PROPERTIES(sundials_nvectrilinos_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})

  INSTALL(TARGETS sundials_nvectrilinos_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECTRILINOS header files
INSTALL(FILES ${nvectrilinos_HEADERS} DESTINATION include/nvector)
INSTALL(FILES ${nvectrilinos_CXXHEADERS} DESTINATION include/nvector/trilinos)

#
MESSAGE(STATUS "Added NVECTOR_Trilinos module")
