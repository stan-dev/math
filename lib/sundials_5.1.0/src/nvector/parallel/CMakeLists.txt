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
# CMakeLists.txt file for the parallel NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_PARALLEL\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Add variable nvecparallel_SOURCES with the sources for the NVECPARALLEL lib
set(nvecparallel_SOURCES nvector_parallel.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECPARALLEL library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecparallel_HEADERS with the exported NVECPARALLEL header files
set(nvecparallel_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_parallel.h
  )

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECPARALLEL library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECPARALLEL library
if(BUILD_STATIC_LIBS)
  add_library(sundials_nvecparallel_static STATIC ${nvecparallel_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_nvecparallel_static
    PROPERTIES OUTPUT_NAME sundials_nvecparallel CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nvecparallel_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECPARALLEL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECPARALLEL library
if(BUILD_SHARED_LIBS)
  add_library(sundials_nvecparallel_shared SHARED ${nvecparallel_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_nvecparallel_shared m)
  endif()

  set_target_properties(sundials_nvecparallel_shared
    PROPERTIES OUTPUT_NAME sundials_nvecparallel CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_nvecparallel_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_nvecparallel_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the NVECPARALLEL header files
install(FILES ${nvecparallel_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled and MPI works, build and install the FNVECPARALLEL library
if(F77_INTERFACE_ENABLE AND MPI_Fortran_FOUND)
  set(fnvecparallel_SOURCES fnvector_parallel.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fnvecparallel_static STATIC ${fnvecparallel_SOURCES})
    set_target_properties(sundials_fnvecparallel_static
      PROPERTIES OUTPUT_NAME sundials_fnvecparallel CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fnvecparallel_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fnvecparallel_shared SHARED ${fnvecparallel_SOURCES})

    # fnvecparallel depends on nvecparallel
    target_link_libraries(sundials_fnvecparallel_shared sundials_nvecparallel_shared)

    set_target_properties(sundials_fnvecparallel_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecparallel CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fnvecparallel_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    install(TARGETS sundials_fnvecparallel_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

endif()
#
message(STATUS "Added NVECTOR_PARALLEL module")
