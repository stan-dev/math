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
# CMakeLists.txt file for the petsc NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_PETSC\n\")")

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Add variable nvecpetsc_SOURCES with the sources for the NVECPARHYP lib
set(nvecpetsc_SOURCES nvector_petsc.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPARHYP library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecpetsc_HEADERS with the exported NVECPARHYP header files
set(nvecpetsc_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_petsc.h
  )

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

#
if(BUILD_STATIC_LIBS)
  add_library(sundials_nvecpetsc_static STATIC ${nvecpetsc_SOURCES} ${shared_SOURCES})

  target_link_libraries(sundials_nvecpetsc_static PUBLIC PETSC::VEC)

  set_target_properties(sundials_nvecpetsc_static PROPERTIES
                        OUTPUT_NAME sundials_nvecpetsc
                        CLEAN_DIRECT_OUTPUT 1)

  install(TARGETS sundials_nvecpetsc_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

#
if(BUILD_SHARED_LIBS)
  add_library(sundials_nvecpetsc_shared SHARED ${nvecpetsc_SOURCES} ${shared_SOURCES})

  target_link_libraries(sundials_nvecpetsc_shared PUBLIC PETSC::VEC)

  if(UNIX)
    target_link_libraries(sundials_nvecpetsc_shared PRIVATE m)
  endif()

  set_target_properties(sundials_nvecpetsc_shared PROPERTIES
                        OUTPUT_NAME sundials_nvecpetsc
                        CLEAN_DIRECT_OUTPUT 1
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION})

  install(TARGETS sundials_nvecpetsc_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the NVECPARHYP header files
install(FILES ${nvecpetsc_HEADERS} DESTINATION include/nvector)

#
message(STATUS "Added NVECTOR_PETSC module")
