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
# CMakeLists.txt file for the superlu-dist SUNLinearSolver library
# -----------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLUDIST\n\")")

# Source files for the library
set(sunlinsolsuperludist_SOURCES sunlinsol_superludist.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c)

# Exported header files
set(sunlinsolsuperludist_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_superludist.h)

# Use MPI
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsolsuperludist_static
    STATIC ${sunlinsolsuperludist_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolsuperludist_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolsuperludist
    CLEAN_DIRECT_OUTPUT 1)

  target_link_libraries(sundials_sunlinsolsuperludist_static
    PUBLIC sundials_sunmatrixslunrloc_static SuperLU_DIST::SuperLU_DIST
    PRIVATE m)

  if(SUPERLUDIST_OpenMP)
    target_compile_options(sundials_sunlinsolsuperludist_static
      PUBLIC ${OpenMP_C_FLAGS})
  endif()

  target_compile_definitions(sundials_sunlinsolsuperludist_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolsuperludist_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsolsuperludist_shared
    SHARED ${sunlinsolsuperludist_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolsuperludist_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolsuperludist
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  target_link_libraries(sundials_sunlinsolsuperludist_shared
    PUBLIC sundials_sunmatrixslunrloc_shared SuperLU_DIST::SuperLU_DIST
    PRIVATE m)

  if(SUPERLUDIST_OpenMP)
    target_compile_options(sundials_sunlinsolsuperludist_shared
      PUBLIC ${OpenMP_C_FLAGS})
  endif()

  target_compile_definitions(sundials_sunlinsolsuperludist_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolsuperludist_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsolsuperludist_HEADERS} DESTINATION include/sunlinsol)

message(STATUS "Added SUNLINSOL_SUPERLUDIST module")
