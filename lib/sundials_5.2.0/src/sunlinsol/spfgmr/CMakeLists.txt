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
# CMakeLists.txt file for the SPFGMR SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SPFGMR\n\")")

# Source files for the library
set(sunlinsolspfgmr_SOURCES sunlinsol_spfgmr.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c)

# Exported header files
set(sunlinsolspfgmr_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_spfgmr.h)

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsolspfgmr_static
    STATIC ${sunlinsolspfgmr_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolspfgmr_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolspfgmr
    CLEAN_DIRECT_OUTPUT 1)

  if(UNIX)
    target_link_libraries(sundials_sunlinsolspfgmr_static
      PRIVATE m)
  endif()

  target_compile_definitions(sundials_sunlinsolspfgmr_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolspfgmr_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsolspfgmr_shared
    SHARED ${sunlinsolspfgmr_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolspfgmr_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolspfgmr
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  if(UNIX)
    target_link_libraries(sundials_sunlinsolspfgmr_shared
      PRIVATE m)
  endif()

  target_compile_definitions(sundials_sunlinsolspfgmr_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolspfgmr_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsolspfgmr_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_SPFGMR module")

# If FCMIX is enabled, build and install the Fortran77 library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(fsunlinsolspfgmr_SOURCES fsunlinsol_spfgmr.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsolspfgmr_static
      STATIC ${fsunlinsolspfgmr_SOURCES})
    set_target_properties(sundials_fsunlinsolspfgmr_static
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolspfgmr
      CLEAN_DIRECT_OUTPUT 1)
    # depends on fnvecserial and sunlinsolspfgmr
    target_link_libraries(sundials_fsunlinsolspfgmr_static
      PUBLIC
      sundials_fnvecserial_static
      sundials_sunlinsolspfgmr_static)
    install(TARGETS sundials_fsunlinsolspfgmr_static
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsolspfgmr_shared
      SHARED ${fsunlinsolspfgmr_SOURCES})
    set_target_properties(sundials_fsunlinsolspfgmr_shared
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolspfgmr
      CLEAN_DIRECT_OUTPUT 1
      VERSION ${sunlinsollib_VERSION}
      SOVERSION ${sunlinsollib_SOVERSION})
    # depends on fnvecserial and sunlinsolspfgmr
    target_link_libraries(sundials_fsunlinsolspfgmr_shared
      PUBLIC
      sundials_fnvecserial_shared
      sundials_sunlinsolspfgmr_shared)
    install(TARGETS sundials_fsunlinsolspfgmr_shared
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

  #
  message(STATUS "Added SUNLINSOL_SPFGMR F77 interface")

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)
