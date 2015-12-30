# ---------------------------------------------------------------
# $Revision: 4075 $
# $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
# ---------------------------------------------------------------
# Programmer:  Steve Smith @ LLNL
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
# CMakeLists.txt file for the pthreads NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_PTHREADS\n\")")

# Add variable nvecpthreads_SOURCES with the sources for the NVECPTHREADS lib
SET(nvecpthreads_SOURCES nvector_pthreads.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPTHREADS library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvecpthreads_HEADERS with the exported NVECPTHREADS header files
SET(nvecpthreads_HEADERS nvector_pthreads.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecpthreads_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECPTHREADS library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECPTHREADS library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecpthreads_static STATIC ${nvecpthreads_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecpthreads_static
    PROPERTIES OUTPUT_NAME sundials_nvecpthreads CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecpthreads_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECPTHREADS library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECPTHREADS library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecpthreads_shared SHARED ${nvecpthreads_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecpthreads_shared
    PROPERTIES OUTPUT_NAME sundials_nvecpthreads CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecpthreads_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecpthreads_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECPTHREADS header files
INSTALL(FILES ${nvecpthreads_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECPTHREADS library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fnvecpthreads_SOURCES fnvector_pthreads.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecpthreads_static STATIC ${fnvecpthreads_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecpthreads_static
      PROPERTIES OUTPUT_NAME sundials_fnvecpthreads CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecpthreads_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecpthreads_shared ${fnvecpthreads_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecpthreads_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecpthreads CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecpthreads_shared 
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecpthreads_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added NVECTOR_PTHREADS module")
