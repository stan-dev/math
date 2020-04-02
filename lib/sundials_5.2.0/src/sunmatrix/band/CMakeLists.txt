# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the band SUNMatrix library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_BAND\n\")")

# Add variable sunmatrixband_SOURCES with the sources for the SUNMATRIXBAND lib
set(sunmatrixband_SOURCES sunmatrix_band.c)

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNMATRIXBAND library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable sunmatrixband_HEADERS with the exported SUNMATRIXBAND header files
set(sunmatrixband_HEADERS
  ${sundials_SOURCE_DIR}/include/sunmatrix/sunmatrix_band.h
  )

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNMATRIXBAND library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNMATRIXBAND library
if(BUILD_STATIC_LIBS)
  add_library(sundials_sunmatrixband_static STATIC ${sunmatrixband_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_sunmatrixband_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixband CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_sunmatrixband_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNMATRIXBAND library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNMATRIXBAND library
if(BUILD_SHARED_LIBS)
  add_library(sundials_sunmatrixband_shared SHARED ${sunmatrixband_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_sunmatrixband_shared m)
  endif()

  set_target_properties(sundials_sunmatrixband_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixband CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_sunmatrixband_shared
    PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
  install(TARGETS sundials_sunmatrixband_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the SUNMATRIXBAND header files
install(FILES ${sunmatrixband_HEADERS} DESTINATION include/sunmatrix)

# If FCMIX is enabled, build and install the FSUNMATRIXBAND library
if(F77_INTERFACE_ENABLE AND F77_FOUND)
  set(fsunmatrixband_SOURCES fsunmatrix_band.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunmatrixband_static STATIC ${fsunmatrixband_SOURCES})
    set_target_properties(sundials_fsunmatrixband_static
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixband CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fsunmatrixband_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunmatrixband_shared ${fsunmatrixband_SOURCES})

    # fsunmatrixband depends on sunmatrixband
    target_link_libraries(sundials_fsunmatrixband_shared sundials_sunmatrixband_shared)

    set_target_properties(sundials_fsunmatrixband_shared
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixband CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fsunmatrixband_shared
      PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
    install(TARGETS sundials_fsunmatrixband_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

#
message(STATUS "Added SUNMATRIX_BAND module")
