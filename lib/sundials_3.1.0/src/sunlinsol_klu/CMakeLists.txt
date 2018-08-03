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
# CMakeLists.txt file for the klu SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_KLU\n\")")

# Add variable sunlinsolklu_SOURCES with the sources for the SUNLINSOLKLU lib
SET(sunlinsolklu_SOURCES sunlinsol_klu.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLKLU library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  )

# Add variable sunlinsolklu_HEADERS with the exported SUNLINSOLKLU header files
SET(sunlinsolklu_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_klu.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLKLU library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLKLU library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolklu_static STATIC ${sunlinsolklu_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolklu_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolklu CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolklu_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLKLU library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLKLU library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolklu_shared SHARED ${sunlinsolklu_SOURCES} ${shared_SOURCES})

  # sunlinsolklu depends on sunmatrixsparse and KLU
  TARGET_LINK_LIBRARIES(sundials_sunlinsolklu_shared
    sundials_sunmatrixsparse_shared
    ${KLU_LIBRARIES})

  SET_TARGET_PROPERTIES(sundials_sunlinsolklu_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolklu CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolklu_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolklu_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLKLU header files
INSTALL(FILES ${sunlinsolklu_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled and building with 64-bit indices, then build and install 
# the FSUNLINSOLKLU library
IF(SUNDIALS_INDEX_TYPE MATCHES "INT32_T")

  PRINT_WARNING("The Fortran interface to KLU is not compatible with ${SUNDIALS_INDEX_TYPE} integers"
                "Fortran interfaces to KLU will not be built")

ELSE(SUNDIALS_INDEX_TYPE MATCHES "INT32_T")

  IF(FCMIX_ENABLE AND F77_FOUND)
    SET(fsunlinsolklu_SOURCES fsunlinsol_klu.c)

    IF(BUILD_STATIC_LIBS)
      ADD_LIBRARY(sundials_fsunlinsolklu_static STATIC ${fsunlinsolklu_SOURCES})
      SET_TARGET_PROPERTIES(sundials_fsunlinsolklu_static
        PROPERTIES OUTPUT_NAME sundials_fsunlinsolklu CLEAN_DIRECT_OUTPUT 1)
      INSTALL(TARGETS sundials_fsunlinsolklu_static DESTINATION lib)
    ENDIF(BUILD_STATIC_LIBS)

    IF(BUILD_SHARED_LIBS)
      ADD_LIBRARY(sundials_fsunlinsolklu_shared ${fsunlinsolklu_SOURCES})

      # fsunlinsolklu depends on fnvecserial, fsunmatrixsparse, sunlinsolklu
      TARGET_LINK_LIBRARIES(sundials_fsunlinsolklu_shared
        sundials_fnvecserial_shared
        sundials_fsunmatrixsparse_shared
        sundials_sunlinsolklu_shared)

      SET_TARGET_PROPERTIES(sundials_fsunlinsolklu_shared
        PROPERTIES OUTPUT_NAME sundials_fsunlinsolklu CLEAN_DIRECT_OUTPUT 1)
      SET_TARGET_PROPERTIES(sundials_fsunlinsolklu_shared 
        PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
      INSTALL(TARGETS sundials_fsunlinsolklu_shared DESTINATION lib)
    ENDIF(BUILD_SHARED_LIBS)

  ENDIF(FCMIX_ENABLE AND F77_FOUND)

ENDIF(SUNDIALS_INDEX_TYPE MATCHES "INT32_T")

#
MESSAGE(STATUS "Added SUNLINSOL_KLU module")
