# ---------------------------------------------------------------
# Programmer: David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the openmpdev NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_OPENMPDEV\n\")")

# Add variable nvecopenmpdev_SOURCES with the sources for the NVECOPENMPDEV lib
set(nvecopenmpdev_SOURCES nvector_openmpdev.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECOPENMPDEV library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecopenmpdev_HEADERS with the exported NVECOPENMPDEV header files
set(nvecopenmpdev_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_openmpdev.h
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
#  - Add the build target for the NVECOPENMPDEV library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECOPENMPDEV library
if(BUILD_STATIC_LIBS)
  add_library(sundials_nvecopenmpdev_static STATIC ${nvecopenmpdev_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_nvecopenmpdev_static
    PROPERTIES OUTPUT_NAME sundials_nvecopenmpdev CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nvecopenmpdev_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Rules for building and installing the shared library:
#  - Add the build target for the NVECOPENMPDEV library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECOPENMPDEV library
if(BUILD_SHARED_LIBS)
  add_library(sundials_nvecopenmpdev_shared SHARED ${nvecopenmpdev_SOURCES} ${shared_SOURCES})

  if(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvecopenmpdev_shared m)
  endif()

  set_target_properties(sundials_nvecopenmpdev_shared
    PROPERTIES OUTPUT_NAME sundials_nvecopenmpdev CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_nvecopenmpdev_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_nvecopenmpdev_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif()

# Install the NVECOPENMPDEV header files
INSTALL(FILES ${nvecopenmpdev_HEADERS} DESTINATION include/nvector)

#
MESSAGE(STATUS "Added NVECTOR_OPENMPDEV module")
