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
# CMakeLists.txt file for the parallel NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_PARALLEL\n\")")

IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Add variable nvecparallel_SOURCES with the sources for the NVECPARALLEL lib
SET(nvecparallel_SOURCES nvector_parallel.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPARALLEL library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_mpi.c
  )

# Add variable nvecparallel_HEADERS with the exported NVECPARALLEL header files
SET(nvecparallel_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_parallel.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECPARALLEL library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECPARALLEL library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecparallel_static STATIC ${nvecparallel_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecparallel_static
    PROPERTIES OUTPUT_NAME sundials_nvecparallel CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecparallel_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECPARALLEL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECPARALLEL library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecparallel_shared SHARED ${nvecparallel_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvecparallel_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_nvecparallel_shared
    PROPERTIES OUTPUT_NAME sundials_nvecparallel CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecparallel_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecparallel_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECPARALLEL header files
INSTALL(FILES ${nvecparallel_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled and MPI works, build and install the FNVECPARALLEL library
IF(F77_INTERFACE_ENABLE AND MPI_Fortran_FOUND)
  SET(fnvecparallel_SOURCES fnvector_parallel.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecparallel_static STATIC ${fnvecparallel_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecparallel_static
      PROPERTIES OUTPUT_NAME sundials_fnvecparallel CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecparallel_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecparallel_shared SHARED ${fnvecparallel_SOURCES})

    # fnvecparallel depends on nvecparallel
    TARGET_LINK_LIBRARIES(sundials_fnvecparallel_shared sundials_nvecparallel_shared)

    SET_TARGET_PROPERTIES(sundials_fnvecparallel_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecparallel CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecparallel_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecparallel_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF()
#
MESSAGE(STATUS "Added NVECTOR_PARALLEL module")
