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
# CMakeLists.txt file for the openmp NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_OPENMP\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable nvecopenmp_SOURCES with the sources for the NVECOPENMP lib
set(nvecopenmp_SOURCES nvector_openmp.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECOPENMP library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecopenmp_HEADERS with the exported NVECOPENMP header files
set(nvecopenmp_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_openmp.h
  )

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# Use C flags for linker as well.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")


# Rules for building and installing the static library:
#  - Add the build target for the NVECOPENMP library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECOPENMP library
if(BUILD_STATIC_LIBS)
  add_library(sundials_nvecopenmp_static STATIC ${nvecopenmp_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_nvecopenmp_static
    PROPERTIES OUTPUT_NAME sundials_nvecopenmp CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nvecopenmp_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECOPENMP library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECOPENMP library
if(BUILD_SHARED_LIBS)
  add_library(sundials_nvecopenmp_shared SHARED ${nvecopenmp_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_nvecopenmp_shared m)
  endif()

  set_target_properties(sundials_nvecopenmp_shared
    PROPERTIES OUTPUT_NAME sundials_nvecopenmp CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_nvecopenmp_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_nvecopenmp_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the NVECOPENMP header files
install(FILES ${nvecopenmp_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECOPENMP library
if(F77_INTERFACE_ENABLE AND F77_FOUND)
  set(fnvecopenmp_SOURCES fnvector_openmp.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fnvecopenmp_static STATIC ${fnvecopenmp_SOURCES})
    set_target_properties(sundials_fnvecopenmp_static
      PROPERTIES OUTPUT_NAME sundials_fnvecopenmp CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fnvecopenmp_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fnvecopenmp_shared ${fnvecopenmp_SOURCES})

    # fnvecopenmp depends on nvecopenmp
    target_link_libraries(sundials_fnvecopenmp_shared sundials_nvecopenmp_shared)

    set_target_properties(sundials_fnvecopenmp_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecopenmp CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fnvecopenmp_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    install(TARGETS sundials_fnvecopenmp_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

#
message(STATUS "Added NVECTOR_OPENMP module")
