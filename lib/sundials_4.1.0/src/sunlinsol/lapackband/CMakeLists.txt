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
# CMakeLists.txt file for the LapackBand SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKBAND\n\")")

# Add variable sunlinsollapackband_SOURCES with the sources for the SUNLINSOLLAPACKBAND lib
SET(sunlinsollapackband_SOURCES sunlinsol_lapackband.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLLAPACKBAND library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  )

# Add variable sunlinsollapackband_HEADERS with the
# exported SUNLINSOLLAPACKBAND header files
SET(sunlinsollapackband_HEADERS 
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_lapackband.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLLAPACKBAND library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLLAPACKBAND library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsollapackband_static STATIC ${sunlinsollapackband_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsollapackband_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsollapackband CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsollapackband_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLLAPACKBAND library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLLAPACKBAND library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsollapackband_shared SHARED ${sunlinsollapackband_SOURCES} ${shared_SOURCES})

  # sunlinsollapackband depends on sunmatrixband and LAPACK
  TARGET_LINK_LIBRARIES(sundials_sunlinsollapackband_shared
    sundials_sunmatrixband_shared
    ${LAPACK_LIBRARIES})

  SET_TARGET_PROPERTIES(sundials_sunlinsollapackband_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsollapackband CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsollapackband_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsollapackband_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLLAPACKBAND header files
INSTALL(FILES ${sunlinsollapackband_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLLAPACKBAND library
IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fsunlinsollapackband_SOURCES fsunlinsol_lapackband.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsollapackband_static STATIC ${fsunlinsollapackband_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackband_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsollapackband CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsollapackband_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsollapackband_shared ${fsunlinsollapackband_SOURCES})

    # fsunlinsollapackband depends on fnvecserial, fsunmatrixband, sunlinsollapackband
    TARGET_LINK_LIBRARIES(sundials_fsunlinsollapackband_shared
      sundials_fnvecserial_shared
      sundials_fsunmatrixband_shared
      sundials_sunlinsollapackband_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackband_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsollapackband CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackband_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsollapackband_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_LAPACKBAND module")
