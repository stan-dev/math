# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# CMakeLists.txt file for the SuperLU SLU_NR_LOC SUNMatrix library
# -----------------------------------------------------------------------------

install(CODE "message(\"\nInstall SUNMATRIX_SLUNRLOC\n\")")

# Add variable sunmatrixslunrloc_SOURCES with the sources for the
# sunmatrixslunrloc lib
set(sunmatrixslunrloc_SOURCES sunmatrix_slunrloc.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the sunmatrixslunrloc library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable sunmatrixslunrloc_HEADERS with the exported sunmatrixslunrloc
# header files
set(sunmatrixslunrloc_HEADERS
  ${sundials_SOURCE_DIR}/include/sunmatrix/sunmatrix_slunrloc.h
  )

# Use MPI
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Rules for building and installing the static library:
#  - Add the build target for the sunmatrixslunrloc library
#  - Set the library name and make sure it is not deleted
#  - Install the sunmatrixslunrloc library
if(BUILD_STATIC_LIBS)
  add_library(sundials_sunmatrixslunrloc_static STATIC
    ${sunmatrixslunrloc_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunmatrixslunrloc_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixslunrloc CLEAN_DIRECT_OUTPUT 1)

  target_link_libraries(sundials_sunmatrixslunrloc_static
    PUBLIC SuperLU_DIST::SuperLU_DIST PRIVATE m)

  if(SUPERLUDIST_OpenMP AND OPENMP_ENABLE AND OPENMP_FOUND)
    target_compile_options(sundials_sunmatrixslunrloc_static
      PUBLIC ${OpenMP_C_FLAGS})
  endif()

  target_compile_definitions(sundials_sunmatrixslunrloc_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunmatrixslunrloc_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the sunmatrixslunrloc library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the sunmatrixslunrloc library
if(BUILD_SHARED_LIBS)
  add_library(sundials_sunmatrixslunrloc_shared SHARED
    ${sunmatrixslunrloc_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunmatrixslunrloc_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixslunrloc CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})

  target_link_libraries(sundials_sunmatrixslunrloc_shared
    PUBLIC SuperLU_DIST::SuperLU_DIST PRIVATE m)

  if(SUPERLUDIST_OpenMP AND OPENMP_ENABLE AND OPENMP_FOUND)
    target_compile_options(sundials_sunmatrixslunrloc_shared
      PUBLIC ${OpenMP_C_FLAGS})
  endif()

  target_compile_definitions(sundials_sunmatrixslunrloc_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunmatrixslunrloc_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the sunmatrixslunrloc header files
install(FILES ${sunmatrixslunrloc_HEADERS} DESTINATION include/sunmatrix)

message(STATUS "Added SUNMATRIX_SLUNRLOC module")
