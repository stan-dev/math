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
# CMakeLists.txt file for the SPBCGS SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_SPBCGS\n\")")

# Add variable sunlinsolspbcgs_SOURCES with the sources for the SUNLINSOLSPBCGS lib
SET(sunlinsolspbcgs_SOURCES sunlinsol_spbcgs.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSPBCGS library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  )

# Add variable sunlinsolspbcgs_HEADERS with the exported SUNLINSOLSPBCGS header files
SET(sunlinsolspbcgs_HEADERS 
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_spbcgs.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSPBCGS library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSPBCGS library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolspbcgs_static STATIC ${sunlinsolspbcgs_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolspbcgs_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolspbcgs CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolspbcgs_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLSPBCGS library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLSPBCGS library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolspbcgs_shared SHARED ${sunlinsolspbcgs_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolspbcgs_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunlinsolspbcgs_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolspbcgs CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolspbcgs_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolspbcgs_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLSPBCGS header files
INSTALL(FILES ${sunlinsolspbcgs_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLSPBCGS library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsolspbcgs_SOURCES fsunlinsol_spbcgs.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolspbcgs_static STATIC ${fsunlinsolspbcgs_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspbcgs_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolspbcgs CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolspbcgs_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolspbcgs_shared ${fsunlinsolspbcgs_SOURCES})

    # fsunlinsolspbcgs depdens on fnvecserial and sunlinsolspbcgs
    TARGET_LINK_LIBRARIES(sundials_fsunlinsolspbcgs_shared
      sundials_fnvecserial_shared
      sundials_sunlinsolspbcgs_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsolspbcgs_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolspbcgs CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspbcgs_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolspbcgs_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_SPBCGS module")
