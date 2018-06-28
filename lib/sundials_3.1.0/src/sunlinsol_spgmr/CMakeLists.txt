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
# CMakeLists.txt file for the SPGMR SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_SPGMR\n\")")

# Add variable sunlinsolspgmr_SOURCES with the sources for the SUNLINSOLSPGMR lib
SET(sunlinsolspgmr_SOURCES sunlinsol_spgmr.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSPGMR library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  )

# Add variable sunlinsolspgmr_HEADERS with the exported SUNLINSOLSPGMR header files
SET(sunlinsolspgmr_HEADERS 
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_spgmr.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSPGMR library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSPGMR library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolspgmr_static STATIC ${sunlinsolspgmr_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolspgmr_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolspgmr CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolspgmr_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLSPGMR library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLSPGMR library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolspgmr_shared SHARED ${sunlinsolspgmr_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolspgmr_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunlinsolspgmr_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolspgmr CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolspgmr_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolspgmr_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLSPGMR header files
INSTALL(FILES ${sunlinsolspgmr_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLSPGMR library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsolspgmr_SOURCES fsunlinsol_spgmr.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolspgmr_static STATIC ${fsunlinsolspgmr_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspgmr_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolspgmr CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolspgmr_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolspgmr_shared ${fsunlinsolspgmr_SOURCES})

    # fsunlinsolspgmr depends on fnvecserial and sunlinsolspgmr
    TARGET_LINK_LIBRARIES(sundials_fsunlinsolspgmr_shared
      sundials_fnvecserial_shared
      sundials_sunlinsolspgmr_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsolspgmr_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolspgmr CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspgmr_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolspgmr_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_SPGMR module")
