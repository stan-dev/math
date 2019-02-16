# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# CMakeLists.txt file for the Newton SUNNonlinearSolver library
# ------------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNNONLINSOL_NEWTON\n\")")

# Add F90 module if F2003 interface is enabled
if(F90_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(F90)
endif(F90_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable with the sources for the library
set(lib_SOURCES sunnonlinsol_newton.c)

# Add variable with the common SUNDIALS sources to be included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector_senswrapper.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nonlinearsolver.c
  )

# Add variable with the exported header files
set(lib_HEADERS
  ${sundials_SOURCE_DIR}/include/sunnonlinsol/sunnonlinsol_newton.h
  )

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)
  add_library(sundials_sunnonlinsolnewton_static STATIC ${lib_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_sunnonlinsolnewton_static
    PROPERTIES OUTPUT_NAME sundials_sunnonlinsolnewton CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_sunnonlinsolnewton_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)
  add_library(sundials_sunnonlinsolnewton_shared SHARED ${lib_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_sunnonlinsolnewton_shared m)
  endif()

  set_target_properties(sundials_sunnonlinsolnewton_shared
    PROPERTIES OUTPUT_NAME sundials_sunnonlinsolnewton CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_sunnonlinsolnewton_shared
    PROPERTIES VERSION ${sunnonlinsollib_VERSION} SOVERSION ${sunnonlinsollib_SOVERSION})
  install(TARGETS sundials_sunnonlinsolnewton_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Install the header files
INSTALL(FILES ${lib_HEADERS} DESTINATION include/sunnonlinsol)

# If FCMIX is enabled, build and install the F library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(flib_SOURCES fsunnonlinsol_newton.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunnonlinsolnewton_static STATIC ${flib_SOURCES})
    set_target_properties(sundials_fsunnonlinsolnewton_static
      PROPERTIES OUTPUT_NAME sundials_fsunnonlinsolnewton CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fsunnonlinsolnewton_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunnonlinsolnewton_shared ${flib_SOURCES})

    # depends on fnvecserial and sunnonlinsolnewton
    target_link_libraries(sundials_fsunnonlinsolnewton_shared
      sundials_fnvecserial_shared
      sundials_sunnonlinsolnewton_shared)

    set_target_properties(sundials_fsunnonlinsolnewton_shared
      PROPERTIES OUTPUT_NAME sundials_fsunnonlinsolnewton CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fsunnonlinsolnewton_shared
      PROPERTIES VERSION ${sunnonlinsollib_VERSION} SOVERSION ${sunnonlinsollib_SOVERSION})
    install(TARGETS sundials_fsunnonlinsolnewton_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()

endif()

#
message(STATUS "Added SUNNONLINSOL_NEWTON module")
