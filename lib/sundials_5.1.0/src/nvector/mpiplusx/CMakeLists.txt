# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the MPIPlusX NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_MPIPLUSX\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

if(MPI_C_FOUND AND MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
elseif(MPI_C_FOUND)
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Add variable nvecmpiplusx_SOURCES with the sources for the NVECMPIPLUSX lib
set(nvecmpiplusx_SOURCES
  nvector_mpiplusx.c
  ${sundials_SOURCE_DIR}/src/nvector/manyvector/nvector_manyvector.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECMPIPLUSX library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecmpiplusx_HEADERS with the exported NVECMPIPLUSX header files
set(nvecmpiplusx_HEADERS ${sundials_SOURCE_DIR}/include/nvector/nvector_manyvector.h)
if(MPI_ENABLE AND MPI_C_FOUND)
  list(APPEND nvecmpiplusx_HEADERS
      ${sundials_SOURCE_DIR}/include/nvector/nvector_mpimanyvector.h
      ${sundials_SOURCE_DIR}/include/nvector/nvector_mpiplusx.h)
endif()

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECMPIPLUSX library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECMPIPLUSX library
if(BUILD_STATIC_LIBS)
  add_library(sundials_nvecmpiplusx_static STATIC ${nvecmpiplusx_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_nvecmpiplusx_static PROPERTIES
                        OUTPUT_NAME sundials_nvecmpiplusx
                        CLEAN_DIRECT_OUTPUT 1)

  install(TARGETS sundials_nvecmpiplusx_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECMPIPLUSX library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECMPIPLUSX library
if(BUILD_SHARED_LIBS)
  add_library(sundials_nvecmpiplusx_shared SHARED ${nvecmpiplusx_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_nvecmpiplusx_shared m)
  endif()

  set_target_properties(sundials_nvecmpiplusx_shared PROPERTIES
                        OUTPUT_NAME sundials_nvecmpiplusx
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION}
                        CLEAN_DIRECT_OUTPUT 1)

  install(TARGETS sundials_nvecmpiplusx_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the NVECMPIPLUSX header files
install(FILES ${nvecmpiplusx_HEADERS} DESTINATION include/nvector)

message(STATUS "Added NVECTOR_MPIPLUSX module")
