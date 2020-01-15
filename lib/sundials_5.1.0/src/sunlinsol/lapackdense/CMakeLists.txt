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
# CMakeLists.txt file for the LapackDense SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKDENSE\n\")")

# Source files for the library
set(sunlinsollapackdense_SOURCES sunlinsol_lapackdense.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c)

# Exported header files
set(sunlinsollapackdense_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_lapackdense.h)

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsollapackdense_static
    STATIC ${sunlinsollapackdense_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsollapackdense_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsollapackdense
    CLEAN_DIRECT_OUTPUT 1)

  # depends on sunmatrixdense and LAPACK
  target_link_libraries(sundials_sunlinsollapackdense_static
    PUBLIC sundials_sunmatrixdense_static ${LAPACK_LIBRARIES})

  target_compile_definitions(sundials_sunlinsollapackdense_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsollapackdense_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsollapackdense_shared
    SHARED ${sunlinsollapackdense_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsollapackdense_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsollapackdense
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  # depends on sunmatrixdense and LAPACK
  target_link_libraries(sundials_sunlinsollapackdense_shared
    PUBLIC sundials_sunmatrixdense_shared ${LAPACK_LIBRARIES})

  target_compile_definitions(sundials_sunlinsollapackdense_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsollapackdense_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsollapackdense_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_LAPACKDENSE module")

# If FCMIX is enabled, build and install the Fortran77 library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(fsunlinsollapackdense_SOURCES fsunlinsol_lapackdense.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsollapackdense_static
      STATIC ${fsunlinsollapackdense_SOURCES})
    set_target_properties(sundials_fsunlinsollapackdense_static
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsollapackdense
      CLEAN_DIRECT_OUTPUT 1)
    # depends on fnvecserial, fsunmatrixdense, sunlinsollapackdense
    target_link_libraries(sundials_fsunlinsollapackdense_static
      PUBLIC
      sundials_fnvecserial_static
      sundials_fsunmatrixdense_static
      sundials_sunlinsollapackdense_static)
    install(TARGETS sundials_fsunlinsollapackdense_static
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsollapackdense_shared
      SHARED ${fsunlinsollapackdense_SOURCES})
    set_target_properties(sundials_fsunlinsollapackdense_shared
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsollapackdense
      CLEAN_DIRECT_OUTPUT 1
      VERSION ${sunlinsollib_VERSION}
      SOVERSION ${sunlinsollib_SOVERSION})
    # depends on fnvecserial, fsunmatrixdense, sunlinsollapackdense
    target_link_libraries(sundials_fsunlinsollapackdense_shared
      PUBLIC
      sundials_fnvecserial_shared
      sundials_fsunmatrixdense_shared
      sundials_sunlinsollapackdense_shared)
    install(TARGETS sundials_fsunlinsollapackdense_shared
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

  #
  message(STATUS "Added SUNLINSOL_LAPACKDENSE F77 interface")

endif(F77_INTERFACE_ENABLE AND F77_FOUND)
