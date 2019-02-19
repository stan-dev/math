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
# CMakeLists.txt file for the serial NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_SERIAL\n\")")

# Add F90 module if F2003 interface is enabled
IF(F90_FOUND AND F2003_INTERFACE_ENABLE)
  ADD_SUBDIRECTORY(F90)
ENDIF(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable nvecserial_SOURCES with the sources for the NVECSERIAL lib
SET(nvecserial_SOURCES nvector_serial.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECSERIAL library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecserial_HEADERS with the exported NVECSERIAL header files
SET(nvecserial_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_serial.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECSERIAL library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECSERIAL library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecserial_static STATIC ${nvecserial_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecserial_static
    PROPERTIES OUTPUT_NAME sundials_nvecserial CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecserial_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECSERIAL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECSERIAL library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecserial_shared SHARED ${nvecserial_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvecserial_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_nvecserial_shared
    PROPERTIES OUTPUT_NAME sundials_nvecserial CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecserial_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecserial_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECSERIAL header files
INSTALL(FILES ${nvecserial_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECSERIAL library
IF(F77_INTERFACE_ENABLE AND F77_FOUND)
  SET(fnvecserial_SOURCES fnvector_serial.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecserial_static STATIC ${fnvecserial_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecserial_static
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecserial_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecserial_shared ${fnvecserial_SOURCES})

    # fnvecserial depends on nvecserial
    TARGET_LINK_LIBRARIES(sundials_fnvecserial_shared sundials_nvecserial_shared)

    SET_TARGET_PROPERTIES(sundials_fnvecserial_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecserial_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecserial_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added NVECTOR_SERIAL module")
