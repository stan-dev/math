# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
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
# CMakeLists.txt file for the petsc NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_PETSC\n\")")

IF(MPI_MPICC)
  # use MPI_MPICC as the compiler
  SET(CMAKE_C_COMPILER ${MPI_MPICC})
ELSE(MPI_MPICC)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_MPICC)

# Add variable nvecpetsc_SOURCES with the sources for the NVECPARHYP lib
SET(nvecpetsc_SOURCES nvector_petsc.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPARHYP library
SET(shared_SOURCES 
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
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
  INSTALL(TARGETS sundials_nvecpetsc_static DESTINATION lib)
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
  INSTALL(TARGETS sundials_nvecpetsc_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECPARHYP header files
INSTALL(FILES ${nvecpetsc_HEADERS} DESTINATION include/nvector)

#
MESSAGE(STATUS "Added NVECTOR_PETSC module")
