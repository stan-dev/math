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
# CMakeLists.txt file for the band SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_BAND\n\")")

# Source files for the library
set(sunlinsolband_SOURCES sunlinsol_band.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_band.c)

# Exported header files
set(sunlinsolband_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_band.h)

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsolband_static
    STATIC ${sunlinsolband_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolband_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolband
    CLEAN_DIRECT_OUTPUT 1)

  # depends on sunmatrixband
  target_link_libraries(sundials_sunlinsolband_static
    PUBLIC sundials_sunmatrixband_static)

  target_compile_definitions(sundials_sunlinsolband_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolband_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsolband_shared
    SHARED ${sunlinsolband_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolband_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolband
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  # depends on sunmatrixband
  target_link_libraries(sundials_sunlinsolband_shared
    PUBLIC sundials_sunmatrixband_shared)

  target_compile_definitions(sundials_sunlinsolband_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolband_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsolband_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_BAND module")

# If FCMIX is enabled, build and install the Fortran77 library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(fsunlinsolband_SOURCES fsunlinsol_band.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsolband_static
      STATIC ${fsunlinsolband_SOURCES})
    set_target_properties(sundials_fsunlinsolband_static
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolband
      CLEAN_DIRECT_OUTPUT 1)
    # depends on fnvecserial, fsunmatrixband, sunlinsolband
    target_link_libraries(sundials_fsunlinsolband_static
      PUBLIC
      sundials_fnvecserial_static
      sundials_fsunmatrixband_static
      sundials_sunlinsolband_static)
    install(TARGETS sundials_fsunlinsolband_static
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsolband_shared
      SHARED ${fsunlinsolband_SOURCES})
    set_target_properties(sundials_fsunlinsolband_shared
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolband
      CLEAN_DIRECT_OUTPUT 1
      VERSION ${sunlinsollib_VERSION}
      SOVERSION ${sunlinsollib_SOVERSION})
    # depends on fnvecserial, fsunmatrixband, sunlinsolband
    target_link_libraries(sundials_fsunlinsolband_shared
      PUBLIC
      sundials_fnvecserial_shared
      sundials_fsunmatrixband_shared
      sundials_sunlinsolband_shared)
    install(TARGETS sundials_fsunlinsolband_shared
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

  #
  message(STATUS "Added SUNLINSOL_BAND F77 interface")

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)
