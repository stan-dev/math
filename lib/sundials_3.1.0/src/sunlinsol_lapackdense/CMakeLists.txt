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
# CMakeLists.txt file for the LapackDense SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKDENSE\n\")")

# Add variable sunlinsollapackdense_SOURCES with the sources for
# the SUNLINSOLLAPACKDENSE lib
SET(sunlinsollapackdense_SOURCES sunlinsol_lapackdense.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLLAPACKDENSE library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  )

# Add variable sunlinsollapackdense_HEADERS with the
# exported SUNLINSOLLAPACKDENSE header files
SET(sunlinsollapackdense_HEADERS 
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_lapackdense.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLLAPACKDENSE library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLLAPACKDENSE library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsollapackdense_static STATIC ${sunlinsollapackdense_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsollapackdense_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsollapackdense CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsollapackdense_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLLAPACKDENSE library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLLAPACKDENSE library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsollapackdense_shared SHARED ${sunlinsollapackdense_SOURCES} ${shared_SOURCES})

  # sunlinsollapackdense depends on sunmatrixdense and LAPACK
  TARGET_LINK_LIBRARIES(sundials_sunlinsollapackdense_shared
    sundials_sunmatrixdense_shared
    ${LAPACK_LIBRARIES})
  
  #	if a blas library was provided, assume LAPACK depends on it
  IF(BLAS_FOUND)
    TARGET_LINK_LIBRARIES(sundials_sunlinsollapackdense_shared ${BLAS_LIBRARIES})
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_sunlinsollapackdense_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsollapackdense CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsollapackdense_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsollapackdense_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLLAPACKDENSE header files
INSTALL(FILES ${sunlinsollapackdense_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLLAPACKDENSE library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsollapackdense_SOURCES fsunlinsol_lapackdense.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsollapackdense_static STATIC ${fsunlinsollapackdense_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackdense_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsollapackdense CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsollapackdense_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsollapackdense_shared ${fsunlinsollapackdense_SOURCES})

    # fsunlinsollapack depends on fnvecserial, fsunmatrixdense, sunlinsollapackdense
    TARGET_LINK_LIBRARIES(sundials_fsunlinsollapackdense_shared
      sundials_fnvecserial_shared
      sundials_fsunmatrixdense_shared
      sundials_sunlinsollapackdense_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackdense_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsollapackdense CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackdense_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsollapackdense_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_LAPACKDENSE module")
