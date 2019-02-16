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
# CMakeLists.txt file for the PCG SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_PCG\n\")")

# Add F90 module if F2003 interface is enabled
IF(F90_FOUND AND F2003_INTERFACE_ENABLE)
  ADD_SUBDIRECTORY(F90)
ENDIF(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable sunlinsolpcg_SOURCES with the sources for the SUNLINSOLPCG lib
SET(sunlinsolpcg_SOURCES sunlinsol_pcg.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLPCG library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  )

# Add variable sunlinsolpcg_HEADERS with the exported SUNLINSOLPCG header files
SET(sunlinsolpcg_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_pcg.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLPCG library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLPCG library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolpcg_static STATIC ${sunlinsolpcg_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolpcg_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolpcg CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolpcg_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLPCG library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLPCG library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolpcg_shared SHARED ${sunlinsolpcg_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolpcg_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunlinsolpcg_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolpcg CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolpcg_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolpcg_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLPCG header files
INSTALL(FILES ${sunlinsolpcg_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLPCG library
IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fsunlinsolpcg_SOURCES fsunlinsol_pcg.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolpcg_static STATIC ${fsunlinsolpcg_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolpcg_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolpcg CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolpcg_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolpcg_shared ${fsunlinsolpcg_SOURCES})

    # fsunlinsolpcg dependes on fnvecserial and sunlinsolpcg_shared
    TARGET_LINK_LIBRARIES(sundials_fsunlinsolpcg_shared
      sundials_fnvecserial_shared
      sundials_sunlinsolpcg_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsolpcg_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolpcg CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolpcg_shared
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolpcg_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_PCG module")
