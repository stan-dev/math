# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2017, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the SPTFQMR SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_SPTFQMR\n\")")

# Add variable sunlinsolsptfqmr_SOURCES with the sources for the SUNLINSOLSPTFQMR lib
SET(sunlinsolsptfqmr_SOURCES sunlinsol_sptfqmr.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSPTFQMR library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  )

# Add variable sunlinsolsptfqmr_HEADERS with the exported SUNLINSOLSPTFQMR header files
SET(sunlinsolsptfqmr_HEADERS 
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_sptfqmr.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSPTFQMR library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSPTFQMR library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolsptfqmr_static STATIC ${sunlinsolsptfqmr_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolsptfqmr_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolsptfqmr CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolsptfqmr_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLSPTFQMR library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLSPTFQMR library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolsptfqmr_shared SHARED ${sunlinsolsptfqmr_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolsptfqmr_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunlinsolsptfqmr_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolsptfqmr CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolsptfqmr_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolsptfqmr_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLSPTFQMR header files
INSTALL(FILES ${sunlinsolsptfqmr_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLSPTFQMR library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsolsptfqmr_SOURCES fsunlinsol_sptfqmr.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolsptfqmr_static STATIC ${fsunlinsolsptfqmr_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolsptfqmr_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolsptfqmr CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolsptfqmr_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolsptfqmr_shared ${fsunlinsolsptfqmr_SOURCES})

    # fsunlinsolsptfqmr depends on fnvecserial and sunlinsolsptfqmr
    TARGET_LINK_LIBRARIES(sundials_fsunlinsolsptfqmr_shared
      sundials_fnvecserial_shared
      sundials_sunlinsolsptfqmr_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsolsptfqmr_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolsptfqmr CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolsptfqmr_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolsptfqmr_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_SPTFQMR module")
