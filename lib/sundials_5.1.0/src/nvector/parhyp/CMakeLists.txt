# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
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
# CMakeLists.txt file for the parhyp NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_PARHYP\n\")")

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Add variable nvecparhyp_SOURCES with the sources for the NVECPARHYP lib
set(nvecparhyp_SOURCES nvector_parhyp.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPARHYP library
set(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
  )
add_prefix(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvecparhyp_HEADERS with the exported NVECPARHYP header files
set(nvecparhyp_HEADERS nvector_parhyp.h)
add_prefix(${sundials_SOURCE_DIR}/include/nvector/ nvecparhyp_HEADERS)

# Add source directory to include directories
include_directories(.)
include_directories(${HYPRE_INCLUDE_DIR})

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECPARHYP library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECPARHYP library
if(BUILD_STATIC_LIBS)
  add_library(sundials_nvecparhyp_static STATIC ${nvecparhyp_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_nvecparhyp_static
    PROPERTIES OUTPUT_NAME sundials_nvecparhyp CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nvecparhyp_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECPARHYP library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECPARHYP library
if(BUILD_SHARED_LIBS)
  add_library(sundials_nvecparhyp_shared SHARED ${nvecparhyp_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_nvecparhyp_shared m)
  endif()

  # nvecparhyp depends on hypre
  target_link_libraries(sundials_nvecparhyp_shared ${HYPRE_LIBRARIES})

  set_target_properties(sundials_nvecparhyp_shared
    PROPERTIES OUTPUT_NAME sundials_nvecparhyp CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_nvecparhyp_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_nvecparhyp_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the NVECPARHYP header files
install(FILES ${nvecparhyp_HEADERS} DESTINATION include/nvector)

#
message(STATUS "Added NVECTOR_PARHYP module")
