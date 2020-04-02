# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the cuda NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_CUDA\n\")")

# Add variable nveccuda_SOURCES with the sources for the NVECUDA lib
set(nveccuda_SOURCES nvector_cuda.cu)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECCUDA library
set(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
)
add_prefix(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nveccuda_HEADERS with the exported NVECUDA header files
set(nveccuda_HEADERS nvector_cuda.h)
add_prefix(${sundials_SOURCE_DIR}/include/nvector/ nveccuda_HEADERS)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECCUDA library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECUDA library
if(BUILD_STATIC_LIBS)
  add_library(sundials_nveccuda_static STATIC ${nveccuda_SOURCES} ${shared_SOURCES})
  target_compile_features(sundials_nveccuda_static PUBLIC cxx_std_11)
  set_target_properties(sundials_nveccuda_static PROPERTIES
                        OUTPUT_NAME sundials_nveccuda
                        CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nveccuda_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECUDA library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECUDA library
if(BUILD_SHARED_LIBS)
  add_library(sundials_nveccuda_shared SHARED ${nveccuda_SOURCES} ${shared_SOURCES})
  target_compile_features(sundials_nveccuda_shared PUBLIC cxx_std_11)
  set_target_properties(sundials_nveccuda_shared PROPERTIES
                        OUTPUT_NAME sundials_nveccuda
                        CLEAN_DIRECT_OUTPUT 1
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_nveccuda_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the CUDA NVector header files
install(FILES ${nveccuda_HEADERS} DESTINATION include/nvector)
install(DIRECTORY ${sundials_SOURCE_DIR}/include/nvector/cuda DESTINATION include/nvector)

message(STATUS "Added NVECTOR_CUDA module")
