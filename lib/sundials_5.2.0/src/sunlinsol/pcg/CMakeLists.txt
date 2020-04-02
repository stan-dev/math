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
# CMakeLists.txt file for the PCG SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_PCG\n\")")

# Source files for the library
set(sunlinsolpcg_SOURCES sunlinsol_pcg.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c)

# Exported header files
set(sunlinsolpcg_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_pcg.h)

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsolpcg_static
    STATIC ${sunlinsolpcg_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolpcg_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolpcg
    CLEAN_DIRECT_OUTPUT 1)

  if(UNIX)
    target_link_libraries(sundials_sunlinsolpcg_static
      PRIVATE m)
  endif()

  target_compile_definitions(sundials_sunlinsolpcg_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolpcg_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsolpcg_shared
    SHARED ${sunlinsolpcg_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolpcg_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolpcg
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  if(UNIX)
    target_link_libraries(sundials_sunlinsolpcg_shared
      PRIVATE m)
  endif()

  target_compile_definitions(sundials_sunlinsolpcg_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolpcg_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsolpcg_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_PCG module")

# If FCMIX is enabled, build and install the Fortran77 library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(fsunlinsolpcg_SOURCES fsunlinsol_pcg.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsolpcg_static
      STATIC ${fsunlinsolpcg_SOURCES})
    set_target_properties(sundials_fsunlinsolpcg_static
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolpcg
      CLEAN_DIRECT_OUTPUT 1)
    # depends on fnvecserial and sunlinsolpcg
    target_link_libraries(sundials_fsunlinsolpcg_static
      PUBLIC
      sundials_fnvecserial_static
      sundials_sunlinsolpcg_static)
    install(TARGETS sundials_fsunlinsolpcg_static
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsolpcg_shared
      SHARED ${fsunlinsolpcg_SOURCES})
    set_target_properties(sundials_fsunlinsolpcg_shared
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolpcg
      CLEAN_DIRECT_OUTPUT 1
      VERSION ${sunlinsollib_VERSION}
      SOVERSION ${sunlinsollib_SOVERSION})
    # depends on fnvecserial and sunlinsolpcg
    target_link_libraries(sundials_fsunlinsolpcg_shared
      PUBLIC
      sundials_fnvecserial_shared
      sundials_sunlinsolpcg_shared)
    install(TARGETS sundials_fsunlinsolpcg_shared
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

  #
  message(STATUS "Added SUNLINSOL_PCG F77 interface")

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)
