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
# CMakeLists.txt file for the ManyVector NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_MANYVECTOR\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable nvecmanyvector_SOURCES with the sources for the NVECMANYVECTOR lib
set(nvecmanyvector_SOURCES nvector_manyvector.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECMANYVECTOR library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecmanyvector_HEADERS with the exported NVECMANYVECTOR header files
set(nvecmanyvector_HEADERS ${sundials_SOURCE_DIR}/include/nvector/nvector_manyvector.h)

if(MPI_ENABLE AND MPI_C_FOUND)
  install(CODE "MESSAGE(\"\nInstall NVECTOR_MPIMANYVECTOR\n\")")
  list(APPEND nvecmanyvector_HEADERS ${sundials_SOURCE_DIR}/include/nvector/nvector_mpimanyvector.h)
  if(MPI_C_COMPILER)
    # use MPI wrapper as the compiler
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  elseif()
    # add MPI_INCLUDE_PATH to include directories
    include_directories(${MPI_INCLUDE_PATH})
  endif()
endif()

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECMANYVECTOR library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECMANYVECTOR library
if(BUILD_STATIC_LIBS)

  # ----------------------- No MPI
  add_library(sundials_nvecmanyvector_static STATIC ${nvecmanyvector_SOURCES} ${shared_SOURCES})

  target_compile_definitions(sundials_nvecmanyvector_static PRIVATE -DSUNDIALS_MPI_ENABLED=0)

  set_target_properties(sundials_nvecmanyvector_static PROPERTIES
                        OUTPUT_NAME sundials_nvecmanyvector
                        CLEAN_DIRECT_OUTPUT 1)

  install(TARGETS sundials_nvecmanyvector_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- With MPI
  if(MPI_ENABLE AND MPI_C_FOUND)
    add_library(sundials_nvecmpimanyvector_static STATIC ${nvecmanyvector_SOURCES} ${shared_SOURCES})

    target_compile_definitions(sundials_nvecmpimanyvector_static PRIVATE -DSUNDIALS_MPI_ENABLED=1)

    set_target_properties(sundials_nvecmpimanyvector_static PROPERTIES
                          OUTPUT_NAME sundials_nvecmpimanyvector
                          CLEAN_DIRECT_OUTPUT 1)

    install(TARGETS sundials_nvecmpimanyvector_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECMANYVECTOR library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECMANYVECTOR library
if(BUILD_SHARED_LIBS)

  # ----------------------- No MPI
  add_library(sundials_nvecmanyvector_shared SHARED ${nvecmanyvector_SOURCES} ${shared_SOURCES})

  target_compile_definitions(sundials_nvecmanyvector_shared PRIVATE -DSUNDIALS_MPI_ENABLED=0)

  if(UNIX)
    target_link_libraries(sundials_nvecmanyvector_shared PRIVATE m)
  endif()

  set_target_properties(sundials_nvecmanyvector_shared PROPERTIES
                        OUTPUT_NAME sundials_nvecmanyvector
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION}
                        CLEAN_DIRECT_OUTPUT 1)

  install(TARGETS sundials_nvecmanyvector_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- With MPI
  if(MPI_ENABLE AND MPI_C_FOUND)
    add_library(sundials_nvecmpimanyvector_shared SHARED ${nvecmanyvector_SOURCES} ${shared_SOURCES})

    target_compile_definitions(sundials_nvecmpimanyvector_shared PRIVATE -DSUNDIALS_MPI_ENABLED=1)

    if(UNIX)
      target_link_libraries(sundials_nvecmpimanyvector_shared PRIVATE m)
    endif()

    set_target_properties(sundials_nvecmpimanyvector_shared PROPERTIES
                          OUTPUT_NAME sundials_nvecmpimanyvector
                          VERSION ${nveclib_VERSION}
                          SOVERSION ${nveclib_SOVERSION}
                          CLEAN_DIRECT_OUTPUT 1)

    install(TARGETS sundials_nvecmpimanyvector_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()

endif(BUILD_SHARED_LIBS)

# Install the NVECMANYVECTOR header files
install(FILES ${nvecmanyvector_HEADERS} DESTINATION include/nvector)

message(STATUS "Added NVECTOR_MANYVECTOR module")

if(MPI_ENABLE AND MPI_C_FOUND)
  message(STATUS "Added NVECTOR_MPIMANYVECTOR module")
endif()
