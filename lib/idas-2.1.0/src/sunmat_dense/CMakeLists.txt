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
# CMakeLists.txt file for the dense SUNMatrix library

INSTALL(CODE "MESSAGE(\"\nInstall SUNMATRIX_DENSE\n\")")

# Add variable sunmatrixdense_SOURCES with the sources for the SUNMATRIXDENSE lib
SET(sunmatrixdense_SOURCES sunmatrix_dense.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNMATRIXDENSE library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable sunmatrixdense_HEADERS with the exported SUNMATRIXDENSE header files
SET(sunmatrixdense_HEADERS 
  ${sundials_SOURCE_DIR}/include/sunmatrix/sunmatrix_dense.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNMATRIXDENSE library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNMATRIXDENSE library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunmatrixdense_static STATIC ${sunmatrixdense_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunmatrixdense_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixdense CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunmatrixdense_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNMATRIXDENSE library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNMATRIXDENSE library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunmatrixdense_shared SHARED ${sunmatrixdense_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunmatrixdense_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunmatrixdense_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixdense CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunmatrixdense_shared
    PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
  INSTALL(TARGETS sundials_sunmatrixdense_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNMATRIXDENSE header files
INSTALL(FILES ${sunmatrixdense_HEADERS} DESTINATION include/sunmatrix)

# If FCMIX is enabled, build and install the FSUNMATRIXDENSE library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunmatrixdense_SOURCES fsunmatrix_dense.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunmatrixdense_static STATIC ${fsunmatrixdense_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunmatrixdense_static
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixdense CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunmatrixdense_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunmatrixdense_shared ${fsunmatrixdense_SOURCES})

    # fsunmatrixdense depends on sunmatrixdense
    TARGET_LINK_LIBRARIES(sundials_fsunmatrixdense_shared sundials_sunmatrixdense_shared)

    SET_TARGET_PROPERTIES(sundials_fsunmatrixdense_shared
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixdense CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunmatrixdense_shared 
      PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
    INSTALL(TARGETS sundials_fsunmatrixdense_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNMATRIX_DENSE module")
