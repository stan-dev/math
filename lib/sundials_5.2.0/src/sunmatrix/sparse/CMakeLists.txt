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
# CMakeLists.txt file for the sparse SUNMatrix library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_SPARSE\n\")")

# Add variable sunmatrixsparse_SOURCES with the sources for the SUNMATRIXSPARSE lib
set(sunmatrixsparse_SOURCES sunmatrix_sparse.c)

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNMATRIXSPARSE library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable sunmatrixsparse_HEADERS with the exported SUNMATRIXSPARSE header files
set(sunmatrixsparse_HEADERS
  ${sundials_SOURCE_DIR}/include/sunmatrix/sunmatrix_sparse.h
  )

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNMATRIXSPARSE library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNMATRIXSPARSE library
if(BUILD_STATIC_LIBS)
  add_library(sundials_sunmatrixsparse_static STATIC ${sunmatrixsparse_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_sunmatrixsparse_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixsparse CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_sunmatrixsparse_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNMATRIXSPARSE library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNMATRIXSPARSE library
if(BUILD_SHARED_LIBS)
  add_library(sundials_sunmatrixsparse_shared SHARED ${sunmatrixsparse_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_sunmatrixsparse_shared m)
  endif()

  set_target_properties(sundials_sunmatrixsparse_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixsparse CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_sunmatrixsparse_shared
    PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
  install(TARGETS sundials_sunmatrixsparse_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the SUNMATRIXSPARSE header files
install(FILES ${sunmatrixsparse_HEADERS} DESTINATION include/sunmatrix)

# If FCMIX is enabled, build and install the FSUNMATRIXSPARSE library
if(F77_INTERFACE_ENABLE AND F77_FOUND)
  set(fsunmatrixsparse_SOURCES fsunmatrix_sparse.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunmatrixsparse_static STATIC ${fsunmatrixsparse_SOURCES})
    set_target_properties(sundials_fsunmatrixsparse_static
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixsparse CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fsunmatrixsparse_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunmatrixsparse_shared ${fsunmatrixsparse_SOURCES})

    # fsunmatrixsparse depends on sunmatrixsparse
    target_link_libraries(sundials_fsunmatrixsparse_shared sundials_sunmatrixsparse_shared)

    set_target_properties(sundials_fsunmatrixsparse_shared
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixsparse CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fsunmatrixsparse_shared
      PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
    install(TARGETS sundials_fsunmatrixsparse_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

#
message(STATUS "Added SUNMATRIX_SPARSE module")
