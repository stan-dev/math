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
# CMakeLists.txt file for the superlumt SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLUMT\n\")")

# Add variable sunlinsolsuperlumt_SOURCES with the sources for the SUNLINSOLSUPERLUMT lib
SET(sunlinsolsuperlumt_SOURCES sunlinsol_superlumt.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSUPERLUMT library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  )

# Add variable sunlinsolsuperlumt_HEADERS with the exported SUNLINSOLSUPERLUMT header files
SET(sunlinsolsuperlumt_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_superlumt.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSUPERLUMT library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSUPERLUMT library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolsuperlumt_static STATIC ${sunlinsolsuperlumt_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolsuperlumt_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolsuperlumt CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolsuperlumt_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLSUPERLUMT library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLSUPERLUMT library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolsuperlumt_shared SHARED ${sunlinsolsuperlumt_SOURCES} ${shared_SOURCES})

  # sunlinsolsuperlumt depends on sunmatrixsparse and SuperLU_MT
  TARGET_LINK_LIBRARIES(sundials_sunlinsolsuperlumt_shared
    sundials_sunmatrixsparse_shared
    ${SUPERLUMT_LIBRARIES})
  
  # if a blas library is provided, assume SuperLU_MT depends on it
  IF(BLAS_FOUND)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolsuperlumt_shared ${BLAS_LIBRARIES})
  ENDIF()

  # if a blas library was not provided but a lapack library was set, assume 
  # SuperLU_MT depends depends on a blas library provided in LAPACK_LIBRARIES
  IF((NOT BLAS_FOUND) AND LAPACK_FOUND)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolsuperlumt_shared ${LAPACK_LIBRARIES})
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunlinsolsuperlumt_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolsuperlumt CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolsuperlumt_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolsuperlumt_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLSUPERLUMT header files
INSTALL(FILES ${sunlinsolsuperlumt_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled and building with 64-bit indices, then build and install 
# the FSUNLINSOLSUPERLUMT library
IF(SUNDIALS_INDEX_TYPE MATCHES "INT32_T")

  PRINT_WARNING("The Fortran interface to SuperLU_MT is not compatible with ${SUNDIALS_INDEX_TYPE} integers"
                "Fortran interfaces to SuperLU_MT will not be built")

ELSE(SUNDIALS_INDEX_TYPE MATCHES "INT32_T")

  IF(FCMIX_ENABLE AND F77_FOUND)
    SET(fsunlinsolsuperlumt_SOURCES fsunlinsol_superlumt.c)

    IF(BUILD_STATIC_LIBS)
      ADD_LIBRARY(sundials_fsunlinsolsuperlumt_static STATIC ${fsunlinsolsuperlumt_SOURCES})
      SET_TARGET_PROPERTIES(sundials_fsunlinsolsuperlumt_static
        PROPERTIES OUTPUT_NAME sundials_fsunlinsolsuperlumt CLEAN_DIRECT_OUTPUT 1)
      INSTALL(TARGETS sundials_fsunlinsolsuperlumt_static DESTINATION lib)
    ENDIF(BUILD_STATIC_LIBS)

    IF(BUILD_SHARED_LIBS)
      ADD_LIBRARY(sundials_fsunlinsolsuperlumt_shared ${fsunlinsolsuperlumt_SOURCES})

      # fsunlinsolsuperlumt dependes on fnvecserial, fsunmatrixsparse, sunlinsolsuperlumt
      TARGET_LINK_LIBRARIES(sundials_fsunlinsolsuperlumt_shared
        sundials_fnvecserial_shared
        sundials_fsunmatrixsparse_shared
        sundials_sunlinsolsuperlumt_shared)

      SET_TARGET_PROPERTIES(sundials_fsunlinsolsuperlumt_shared
        PROPERTIES OUTPUT_NAME sundials_fsunlinsolsuperlumt CLEAN_DIRECT_OUTPUT 1)
      SET_TARGET_PROPERTIES(sundials_fsunlinsolsuperlumt_shared 
        PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
      INSTALL(TARGETS sundials_fsunlinsolsuperlumt_shared DESTINATION lib)
    ENDIF(BUILD_SHARED_LIBS)

  ENDIF(FCMIX_ENABLE AND F77_FOUND)

ENDIF(SUNDIALS_INDEX_TYPE MATCHES "INT32_T")

#
MESSAGE(STATUS "Added SUNLINSOL_SUPERLUMT module")
