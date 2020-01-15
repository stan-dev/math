# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the dense SUNMatrix library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_DENSE\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable sunmatrixdense_SOURCES with the sources for the SUNMATRIXDENSE lib
set(sunmatrixdense_SOURCES sunmatrix_dense.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNMATRIXDENSE library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable sunmatrixdense_HEADERS with the exported SUNMATRIXDENSE header files
set(sunmatrixdense_HEADERS
  ${sundials_SOURCE_DIR}/include/sunmatrix/sunmatrix_dense.h
  )

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNMATRIXDENSE library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNMATRIXDENSE library
if(BUILD_STATIC_LIBS)
  add_library(sundials_sunmatrixdense_static STATIC ${sunmatrixdense_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_sunmatrixdense_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixdense CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_sunmatrixdense_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNMATRIXDENSE library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNMATRIXDENSE library
if(BUILD_SHARED_LIBS)
  add_library(sundials_sunmatrixdense_shared SHARED ${sunmatrixdense_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_sunmatrixdense_shared m)
  endif()

  set_target_properties(sundials_sunmatrixdense_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixdense CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_sunmatrixdense_shared
    PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
  install(TARGETS sundials_sunmatrixdense_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the SUNMATRIXDENSE header files
install(FILES ${sunmatrixdense_HEADERS} DESTINATION include/sunmatrix)

# If FCMIX is enabled, build and install the FSUNMATRIXDENSE library
if(F77_INTERFACE_ENABLE AND F77_FOUND)
  set(fsunmatrixdense_SOURCES fsunmatrix_dense.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunmatrixdense_static STATIC ${fsunmatrixdense_SOURCES})
    set_target_properties(sundials_fsunmatrixdense_static
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixdense CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fsunmatrixdense_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunmatrixdense_shared ${fsunmatrixdense_SOURCES})

    # fsunmatrixdense depends on sunmatrixdense
    target_link_libraries(sundials_fsunmatrixdense_shared sundials_sunmatrixdense_shared)

    set_target_properties(sundials_fsunmatrixdense_shared
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixdense CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fsunmatrixdense_shared
      PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
    install(TARGETS sundials_fsunmatrixdense_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

#
message(STATUS "Added SUNMATRIX_DENSE module")
