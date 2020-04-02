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
# CMakeLists.txt file for the superlumt SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLUMT\n\")")

# Source files for the library
set(sunlinsolsuperlumt_SOURCES sunlinsol_superlumt.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c)

# Exported header files
set(sunlinsolsuperlumt_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_superlumt.h)

# Include OpenMP flags if SuperLU_MT is using OpenMP
if(SUPERLUMT_THREAD_TYPE STREQUAL "OPENMP")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsolsuperlumt_static
    STATIC ${sunlinsolsuperlumt_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolsuperlumt_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolsuperlumt
    CLEAN_DIRECT_OUTPUT 1)

  # depends on sunmatrixsparse and SuperLU_MT
  target_link_libraries(sundials_sunlinsolsuperlumt_static
    PUBLIC sundials_sunmatrixsparse_static SuperLU_MT::SuperLU_MT)

  # link to Pthread library if SuperLU_MT is using Pthreads
  if(SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
    target_link_libraries(sundials_sunlinsolsuperlumt_static
      PUBLIC Threads::Threads)
  endif()

  target_compile_definitions(sundials_sunlinsolsuperlumt_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolsuperlumt_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsolsuperlumt_shared
    SHARED ${sunlinsolsuperlumt_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolsuperlumt_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolsuperlumt
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  # depends on sunmatrixsparse and SuperLU_MT
  target_link_libraries(sundials_sunlinsolsuperlumt_shared
    PUBLIC sundials_sunmatrixsparse_shared SuperLU_MT::SuperLU_MT)

  # link to Pthread library if SuperLU_MT is using Pthreads
  if(SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
    target_link_libraries(sundials_sunlinsolsuperlumt_shared
      PUBLIC Threads::Threads)
  endif()

  target_compile_definitions(sundials_sunlinsolsuperlumt_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolsuperlumt_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsolsuperlumt_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_SUPERLUMT module")

# If FCMIX is enabled, build and install the Fortran77 library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(fsunlinsolsuperlumt_SOURCES fsunlinsol_superlumt.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsolsuperlumt_static
      STATIC ${fsunlinsolsuperlumt_SOURCES})
    set_target_properties(sundials_fsunlinsolsuperlumt_static
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolsuperlumt
      CLEAN_DIRECT_OUTPUT 1)
    # depends on fnvecserial, fsunmatrixsparse, sunlinsolsuperlumt
    target_link_libraries(sundials_fsunlinsolsuperlumt_static
      PUBLIC
      sundials_fnvecserial_static
      sundials_fsunmatrixsparse_static
      sundials_sunlinsolsuperlumt_static)
    install(TARGETS sundials_fsunlinsolsuperlumt_static
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsolsuperlumt_shared
      SHARED ${fsunlinsolsuperlumt_SOURCES})
    set_target_properties(sundials_fsunlinsolsuperlumt_shared
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolsuperlumt
      CLEAN_DIRECT_OUTPUT 1
      VERSION ${sunlinsollib_VERSION}
      SOVERSION ${sunlinsollib_SOVERSION})
    # depends on fnvecserial, fsunmatrixsparse, sunlinsolsuperlumt
    target_link_libraries(sundials_fsunlinsolsuperlumt_shared
      PUBLIC
      sundials_fnvecserial_shared
      sundials_fsunmatrixsparse_shared
      sundials_sunlinsolsuperlumt_shared)
    install(TARGETS sundials_fsunlinsolsuperlumt_shared
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

  #
  message(STATUS "Added SUNLINSOL_SUPERLUMT F77 interface")

endif(F77_INTERFACE_ENABLE AND F77_FOUND)
