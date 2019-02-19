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
# CMakeLists.txt file for the SPFGMR SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_SPFGMR\n\")")

# Add F90 module if F2003 interface is enabled
IF(F90_FOUND AND F2003_INTERFACE_ENABLE)
  ADD_SUBDIRECTORY(F90)
ENDIF(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable sunlinsolspfgmr_SOURCES with the sources for the SUNLINSOLSPFGMR lib
SET(sunlinsolspfgmr_SOURCES sunlinsol_spfgmr.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSPFGMR library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  )

# Add variable sunlinsolspfgmr_HEADERS with the exported SUNLINSOLSPFGMR header files
SET(sunlinsolspfgmr_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_spfgmr.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSPFGMR library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSPFGMR library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolspfgmr_static STATIC ${sunlinsolspfgmr_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolspfgmr_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolspfgmr_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLSPFGMR library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLSPFGMR library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolspfgmr_shared SHARED ${sunlinsolspfgmr_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolspfgmr_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunlinsolspfgmr_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolspfgmr_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolspfgmr_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLSPFGMR header files
INSTALL(FILES ${sunlinsolspfgmr_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLSPFGMR library
IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fsunlinsolspfgmr_SOURCES fsunlinsol_spfgmr.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolspfgmr_static STATIC ${fsunlinsolspfgmr_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspfgmr_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolspfgmr_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolspfgmr_shared ${fsunlinsolspfgmr_SOURCES})

    # fsunlinsolspfgmr depends on fnvecserial and sunlinsolspfgmr
    TARGET_LINK_LIBRARIES(sundials_fsunlinsolspfgmr_shared
      sundials_fnvecserial_shared
      sundials_sunlinsolspfgmr_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsolspfgmr_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspfgmr_shared
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolspfgmr_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_SPFGMR module")
