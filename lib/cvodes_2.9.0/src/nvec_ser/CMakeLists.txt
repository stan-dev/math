# ---------------------------------------------------------------
# $Revision: 4075 $
# $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
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
# CMakeLists.txt file for the serial NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_SERIAL\n\")")

# Add variable nvecserial_SOURCES with the sources for the NVECSERIAL lib
SET(nvecserial_SOURCES nvector_serial.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECSERIAL library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvecserial_HEADERS with the exported NVECSERIAL header files
SET(nvecserial_HEADERS nvector_serial.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecserial_HEADERS)

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
  INSTALL(TARGETS sundials_nvecserial_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECSERIAL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECSERIAL library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecserial_shared SHARED ${nvecserial_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecserial_shared
    PROPERTIES OUTPUT_NAME sundials_nvecserial CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecserial_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecserial_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECSERIAL header files
INSTALL(FILES ${nvecserial_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECSERIAL library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fnvecserial_SOURCES fnvector_serial.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecserial_static STATIC ${fnvecserial_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecserial_static
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecserial_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecserial_shared ${fnvecserial_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecserial_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecserial_shared 
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecserial_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added NVECTOR_SERIAL module")
