# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the band SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_BAND\n\")")

# Add F90 module if F2003 interface is enabled
IF(F90_FOUND AND F2003_INTERFACE_ENABLE)
  ADD_SUBDIRECTORY(F90)
ENDIF(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable sunlinsolband_SOURCES with the sources for the SUNLINSOLBAND lib
SET(sunlinsolband_SOURCES sunlinsol_band.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLBAND library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_band.c
  )

# Add variable sunlinsolband_HEADERS with the exported SUNLINSOLBAND header files
SET(sunlinsolband_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_band.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLBAND library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLBAND library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolband_static STATIC ${sunlinsolband_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolband_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolband CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolband_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLBAND library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLBAND library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolband_shared SHARED ${sunlinsolband_SOURCES} ${shared_SOURCES})

  # sunlinsolband depends on sunmatrixband
  TARGET_LINK_LIBRARIES(sundials_sunlinsolband_shared sundials_sunmatrixband_shared)

  SET_TARGET_PROPERTIES(sundials_sunlinsolband_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolband CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolband_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolband_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLBAND header files
INSTALL(FILES ${sunlinsolband_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLBAND library
IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fsunlinsolband_SOURCES fsunlinsol_band.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolband_static STATIC ${fsunlinsolband_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolband_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolband CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolband_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolband_shared ${fsunlinsolband_SOURCES})

    # fsunlinsolband depends on fnvecserial, fsunmatrixband, sunlinsolband
    TARGET_LINK_LIBRARIES(sundials_fsunlinsolband_shared
      sundials_fnvecserial_shared
      sundials_fsunmatrixband_shared
      sundials_sunlinsolband_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsolband_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolband CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolband_shared
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolband_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_BAND module")
