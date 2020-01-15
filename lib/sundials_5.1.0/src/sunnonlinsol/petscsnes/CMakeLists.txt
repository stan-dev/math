# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# CMakeLists.txt file for the PETSc SUNNonlinearSolver library
# ------------------------------------------------------------------------------
install(CODE "MESSAGE(\"\nInstall SUNNONLINSOL_PETSCSNES\n\")")

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Add variable with the sources for the library
set(lib_SOURCES sunnonlinsol_petscsnes.c)

# Add variable with the common SUNDIALS sources to be included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nonlinearsolver.c
  )

# Add variable with the exported header files
set(lib_HEADERS
  ${sundials_SOURCE_DIR}/include/sunnonlinsol/sunnonlinsol_petscsnes.h
  )

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

if(BUILD_STATIC_LIBS)
  add_library(sundials_sunnonlinsolpetscsnes_static STATIC ${lib_SOURCES} ${shared_SOURCES})

  target_link_libraries(sundials_sunnonlinsolpetscsnes_static
    PUBLIC sundials_nvecpetsc_static PETSC::SNES)

  set_target_properties(sundials_sunnonlinsolpetscsnes_static PROPERTIES
                        OUTPUT_NAME sundials_sunnonlinsolpetscsnes
                        CLEAN_DIRECT_OUTPUT 1)

  install(TARGETS sundials_sunnonlinsolpetscsnes_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(BUILD_SHARED_LIBS)
  add_library(sundials_sunnonlinsolpetscsnes_shared SHARED ${lib_SOURCES} ${shared_SOURCES})

  target_link_libraries(sundials_sunnonlinsolpetscsnes_shared
    PUBLIC sundials_nvecpetsc_shared PETSC::SNES)

  if(UNIX)
    target_link_libraries(sundials_sunnonlinsolpetscsnes_shared PRIVATE m)
  endif()

  set_target_properties(sundials_sunnonlinsolpetscsnes_shared PROPERTIES
                        OUTPUT_NAME sundials_sunnonlinsolpetscsnes
                        CLEAN_DIRECT_OUTPUT 1
                        VERSION ${sunnonlinsollib_VERSION}
                        SOVERSION ${sunnonlinsollib_SOVERSION})

  install(TARGETS sundials_sunnonlinsolpetscsnes_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Install the header files
install(FILES ${lib_HEADERS} DESTINATION include/sunnonlinsol)

message(STATUS "Added SUNNONLINSOL_PETSCSNES module")
