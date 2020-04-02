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
# CMakeLists.txt file for the serial NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_SERIAL\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable nvecserial_SOURCES with the sources for the NVECSERIAL lib
set(nvecserial_SOURCES nvector_serial.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECSERIAL library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecserial_HEADERS with the exported NVECSERIAL header files
set(nvecserial_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_serial.h
  )

# Add source directory to include directories
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECSERIAL library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECSERIAL library
if(BUILD_STATIC_LIBS)
  add_library(sundials_nvecserial_static STATIC ${nvecserial_SOURCES} ${shared_SOURCES})
  set_target_properties(sundials_nvecserial_static
    PROPERTIES OUTPUT_NAME sundials_nvecserial CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nvecserial_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECSERIAL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECSERIAL library
if(BUILD_SHARED_LIBS)
  add_library(sundials_nvecserial_shared SHARED ${nvecserial_SOURCES} ${shared_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_nvecserial_shared m)
  endif()

  set_target_properties(sundials_nvecserial_shared
    PROPERTIES OUTPUT_NAME sundials_nvecserial CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_nvecserial_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_nvecserial_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the NVECSERIAL header files
install(FILES ${nvecserial_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled, build and install the FNVECSERIAL library
if(F77_INTERFACE_ENABLE AND F77_FOUND)
  set(fnvecserial_SOURCES fnvector_serial.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fnvecserial_static STATIC ${fnvecserial_SOURCES})
    set_target_properties(sundials_fnvecserial_static
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fnvecserial_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fnvecserial_shared ${fnvecserial_SOURCES})

    # fnvecserial depends on nvecserial
    target_link_libraries(sundials_fnvecserial_shared sundials_nvecserial_shared)

    set_target_properties(sundials_fnvecserial_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fnvecserial_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    install(TARGETS sundials_fnvecserial_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

#
message(STATUS "Added NVECTOR_SERIAL module")
