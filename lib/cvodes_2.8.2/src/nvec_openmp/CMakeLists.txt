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
# CMakeLists.txt file for the openmp NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_OPENMP\n\")")

# Add variable nvecopenmp_SOURCES with the sources for the NVECOPENMP lib
SET(nvecopenmp_SOURCES nvector_openmp.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECOPENMP library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvecopenmp_HEADERS with the exported NVECOPENMP header files
SET(nvecopenmp_HEADERS nvector_openmp.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecopenmp_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  # Use C flags for linker as well.
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
ENDIF(OPENMP_FOUND)

# Rules for building and installing the static library:
#  - Add the build target for the NVECOPENMP library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECOPENMP library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecopenmp_static STATIC ${nvecopenmp_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecopenmp_static
    PROPERTIES OUTPUT_NAME sundials_nvecopenmp CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecopenmp_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECOPENMP library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECOPENMP library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecopenmp_shared SHARED ${nvecopenmp_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecopenmp_shared
    PROPERTIES OUTPUT_NAME sundials_nvecopenmp CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecopenmp_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecopenmp_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECOPENMP header files
INSTALL(FILES ${nvecopenmp_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECOPENMP library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fnvecopenmp_SOURCES fnvector_openmp.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecopenmp_static STATIC ${fnvecopenmp_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecopenmp_static
      PROPERTIES OUTPUT_NAME sundials_fnvecopenmp CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecopenmp_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecopenmp_shared ${fnvecopenmp_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecopenmp_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecopenmp CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecopenmp_shared 
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecopenmp_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added NVECTOR_OPENMP module")
