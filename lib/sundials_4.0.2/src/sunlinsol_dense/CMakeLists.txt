# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds, Ashley Crawford @ SMU
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
# CMakeLists.txt file for the dense SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_DENSE\n\")")

# Add F90 module if F2003 interface is enabled
IF(F90_FOUND AND F2003_INTERFACE_ENABLE)
  ADD_SUBDIRECTORY(F90)
ENDIF(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable sunlinsoldense_SOURCES with the sources for the SUNLINSOLDENSE lib
SET(sunlinsoldense_SOURCES sunlinsol_dense.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLDENSE library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_dense.c
  )

# Add variable sunlinsoldense_HEADERS with the exported SUNLINSOLDENSE header files
SET(sunlinsoldense_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_dense.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLDENSE library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLDENSE library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsoldense_static STATIC ${sunlinsoldense_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsoldense_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsoldense CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsoldense_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLDENSE library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLDENSE library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsoldense_shared SHARED ${sunlinsoldense_SOURCES} ${shared_SOURCES})

  TARGET_LINK_LIBRARIES(sundials_sunlinsoldense_shared sundials_sunmatrixdense_shared)

  SET_TARGET_PROPERTIES(sundials_sunlinsoldense_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsoldense CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsoldense_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsoldense_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLDENSE header files
INSTALL(FILES ${sunlinsoldense_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLDENSE library
IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fsunlinsoldense_SOURCES fsunlinsol_dense.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsoldense_static STATIC ${fsunlinsoldense_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsoldense_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsoldense CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsoldense_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsoldense_shared ${fsunlinsoldense_SOURCES})

    # fsunlinsoldense depends on fnvecserial, fsunmatrixdense, sunlinsoldense
    TARGET_LINK_LIBRARIES(sundials_fsunlinsoldense_shared
      sundials_fnvecserial_shared
      sundials_fsunmatrixdense_shared
      sundials_sunlinsoldense_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsoldense_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsoldense CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsoldense_shared
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsoldense_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_DENSE module")
