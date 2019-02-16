# ---------------------------------------------------------------
# Programmer:  Steve Smith, and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the openmp NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_OPENMP\n\")")

# Add F90 module if F2003 interface is enabled
IF(F90_FOUND AND F2003_INTERFACE_ENABLE)
  ADD_SUBDIRECTORY(F90)
ENDIF(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable nvecopenmp_SOURCES with the sources for the NVECOPENMP lib
SET(nvecopenmp_SOURCES nvector_openmp.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECOPENMP library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecopenmp_HEADERS with the exported NVECOPENMP header files
SET(nvecopenmp_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_openmp.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# Use C flags for linker as well.
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")


# Rules for building and installing the static library:
#  - Add the build target for the NVECOPENMP library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECOPENMP library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecopenmp_static STATIC ${nvecopenmp_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecopenmp_static
    PROPERTIES OUTPUT_NAME sundials_nvecopenmp CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecopenmp_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECOPENMP library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECOPENMP library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecopenmp_shared SHARED ${nvecopenmp_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvecopenmp_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_nvecopenmp_shared
    PROPERTIES OUTPUT_NAME sundials_nvecopenmp CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecopenmp_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecopenmp_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECOPENMP header files
INSTALL(FILES ${nvecopenmp_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECOPENMP library
IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fnvecopenmp_SOURCES fnvector_openmp.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecopenmp_static STATIC ${fnvecopenmp_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecopenmp_static
      PROPERTIES OUTPUT_NAME sundials_fnvecopenmp CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecopenmp_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecopenmp_shared ${fnvecopenmp_SOURCES})

    # fnvecopenmp depends on nvecopenmp
    TARGET_LINK_LIBRARIES(sundials_fnvecopenmp_shared sundials_nvecopenmp_shared)

    SET_TARGET_PROPERTIES(sundials_fnvecopenmp_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecopenmp CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecopenmp_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecopenmp_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added NVECTOR_OPENMP module")
