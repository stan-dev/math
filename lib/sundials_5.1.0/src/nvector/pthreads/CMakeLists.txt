# ---------------------------------------------------------------
# Programmer(s): Steve Smith and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the pthreads NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_PTHREADS\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable nvecpthreads_SOURCES with the sources for the NVECPTHREADS lib
set(nvecpthreads_SOURCES nvector_pthreads.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPTHREADS library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecpthreads_HEADERS with the exported NVECPTHREADS header files
set(nvecpthreads_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_pthreads.h
  )

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECPTHREADS library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECPTHREADS library
if(BUILD_STATIC_LIBS)
  add_library(sundials_nvecpthreads_static STATIC ${nvecpthreads_SOURCES} ${shared_SOURCES})

  target_link_libraries(sundials_nvecpthreads_static Threads::Threads)
  if(UNIX)
    target_link_libraries(sundials_nvecpthreads_static m)
  endif()

  set_target_properties(sundials_nvecpthreads_static
    PROPERTIES OUTPUT_NAME sundials_nvecpthreads CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nvecpthreads_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECPTHREADS library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECPTHREADS library
if(BUILD_SHARED_LIBS)
  add_library(sundials_nvecpthreads_shared SHARED ${nvecpthreads_SOURCES} ${shared_SOURCES})

  target_link_libraries(sundials_nvecpthreads_shared Threads::Threads)
  if(UNIX)
    target_link_libraries(sundials_nvecpthreads_shared m)
  endif()

  set_target_properties(sundials_nvecpthreads_shared
    PROPERTIES OUTPUT_NAME sundials_nvecpthreads CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_nvecpthreads_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_nvecpthreads_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the NVECPTHREADS header files
install(FILES ${nvecpthreads_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECPTHREADS library
if(F77_INTERFACE_ENABLE AND F77_FOUND)
  set(fnvecpthreads_SOURCES fnvector_pthreads.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fnvecpthreads_static STATIC ${fnvecpthreads_SOURCES})
    set_target_properties(sundials_fnvecpthreads_static
      PROPERTIES OUTPUT_NAME sundials_fnvecpthreads CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fnvecpthreads_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fnvecpthreads_shared ${fnvecpthreads_SOURCES})

    # fnvecpthreads depends on nvecpthreads
    target_link_libraries(sundials_fnvecpthreads_shared sundials_nvecpthreads_shared)

    set_target_properties(sundials_fnvecpthreads_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecpthreads CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fnvecpthreads_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    install(TARGETS sundials_fnvecpthreads_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)
endif(F77_INTERFACE_ENABLE AND F77_FOUND)

#
message(STATUS "Added NVECTOR_PTHREADS module")
