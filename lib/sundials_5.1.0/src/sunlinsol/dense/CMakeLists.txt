# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds and Ashley Crawford @ SMU
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
# CMakeLists.txt file for the dense SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_DENSE\n\")")

# Source files for the library
set(sunlinsoldense_SOURCES sunlinsol_dense.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_dense.c)

# Exported header files
set(sunlinsoldense_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_dense.h)

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsoldense_static
    STATIC ${sunlinsoldense_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsoldense_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsoldense
    CLEAN_DIRECT_OUTPUT 1)

  # sunlinsoldense depends on sunmatrixdense
  target_link_libraries(sundials_sunlinsoldense_static
    PUBLIC sundials_sunmatrixdense_static)

  target_compile_definitions(sundials_sunlinsoldense_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsoldense_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsoldense_shared
    SHARED ${sunlinsoldense_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsoldense_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsoldense
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  # sunlinsoldense depends on sunmatrixdense
  target_link_libraries(sundials_sunlinsoldense_shared
    PUBLIC sundials_sunmatrixdense_shared)

  target_compile_definitions(sundials_sunlinsoldense_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsoldense_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsoldense_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_DENSE module")

# If FCMIX is enabled, build and install the Fortran77 library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(fsunlinsoldense_SOURCES fsunlinsol_dense.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsoldense_static
      STATIC ${fsunlinsoldense_SOURCES})
    set_target_properties(sundials_fsunlinsoldense_static
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsoldense
      CLEAN_DIRECT_OUTPUT 1)
    # depends on fnvecserial, fsunmatrixdense, sunlinsoldense
    target_link_libraries(sundials_fsunlinsoldense_static
      PUBLIC
      sundials_fnvecserial_static
      sundials_fsunmatrixdense_static
      sundials_sunlinsoldense_static)
    install(TARGETS sundials_fsunlinsoldense_static
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsoldense_shared
      SHARED ${fsunlinsoldense_SOURCES})
    set_target_properties(sundials_fsunlinsoldense_shared
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsoldense
      CLEAN_DIRECT_OUTPUT 1
      VERSION ${sunlinsollib_VERSION}
      SOVERSION ${sunlinsollib_SOVERSION})
    # depends on fnvecserial, fsunmatrixdense, sunlinsoldense
    target_link_libraries(sundials_fsunlinsoldense_shared
      PUBLIC
      sundials_fnvecserial_shared
      sundials_fsunmatrixdense_shared
      sundials_sunlinsoldense_shared)
    install(TARGETS sundials_fsunlinsoldense_shared
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

  #
  message(STATUS "Added SUNLINSOL_DENSE F77 interface")

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)
