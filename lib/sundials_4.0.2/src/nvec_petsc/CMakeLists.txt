# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
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
# CMakeLists.txt file for the petsc NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_PETSC\n\")")

IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Add variable nvecpetsc_SOURCES with the sources for the NVECPARHYP lib
SET(nvecpetsc_SOURCES nvector_petsc.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPARHYP library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_mpi.c
  )

# Add variable nvecpetsc_HEADERS with the exported NVECPARHYP header files
SET(nvecpetsc_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_petsc.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR})

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECPARHYP library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECPARHYP library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecpetsc_static STATIC ${nvecpetsc_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecpetsc_static
    PROPERTIES OUTPUT_NAME sundials_nvecpetsc CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecpetsc_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECPARHYP library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECPARHYP library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecpetsc_shared SHARED ${nvecpetsc_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvecpetsc_shared m)
  ENDIF()

  # nvecpetsc depends on PETSc
  TARGET_LINK_LIBRARIES(sundials_nvecpetsc_shared ${PETSC_LIBRARIES})

  SET_TARGET_PROPERTIES(sundials_nvecpetsc_shared
    PROPERTIES OUTPUT_NAME sundials_nvecpetsc CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecpetsc_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecpetsc_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECPARHYP header files
INSTALL(FILES ${nvecpetsc_HEADERS} DESTINATION include/nvector)

#
MESSAGE(STATUS "Added NVECTOR_PETSC module")
