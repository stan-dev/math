# ---------------------------------------------------------------
# Programmer:  Slaven Peles, and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the cuda NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_CUDA\n\")")

# Add variable nveccuda_SOURCES with the sources for the NVECSERIAL lib
SET(nveccuda_SOURCES nvector_cuda.cu)

# Tell compiler it is a CUDA source
set_source_files_properties(${nveccuda_SOURCES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECCUDA library
SET(shared_SOURCES
  sundials_math.c
  sundials_mpi.c
)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nveccuda_HEADERS with the exported NVECSERIAL header files
SET(nveccuda_HEADERS nvector_cuda.h)
IF (MPI_ENABLE)
  LIST(APPEND nveccuda_HEADERS nvector_mpicuda.h)
ENDIF()
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nveccuda_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ${MPI_CXX_INCLUDE_PATH})

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECCUDA library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECSERIAL library
IF(BUILD_STATIC_LIBS)
  # ----------------------- CUDA only
  # The FindCUDA module does not properly forward compile options using target_* commands.
  # So as long as we use it (required in CMake < 3.8), we have to manually add compile
  # options in the CUDA_ADD_* commands.
  CUDA_ADD_LIBRARY(sundials_nveccuda_static STATIC ${nveccuda_SOURCES} ${shared_SOURCES}
                   OPTIONS -DSUNDIALS_MPI_ENABLED=0)
  TARGET_COMPILE_DEFINITIONS(sundials_nveccuda_static PUBLIC -DSUNDIALS_MPI_ENABLED=0)

  SET_TARGET_PROPERTIES(sundials_nveccuda_static
    PROPERTIES OUTPUT_NAME sundials_nveccuda CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nveccuda_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- MPI+CUDA
  IF(MPI_ENABLE)
    CUDA_ADD_LIBRARY(sundials_nvecmpicuda_static STATIC ${nveccuda_SOURCES} ${shared_SOURCES})
    SET_TARGET_PROPERTIES(sundials_nvecmpicuda_static
      PROPERTIES OUTPUT_NAME sundials_nvecmpicuda CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_nvecmpicuda_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF()
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECSERIAL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECSERIAL library
IF(BUILD_SHARED_LIBS)
  # ----------------------- CUDA only
  CUDA_ADD_LIBRARY(sundials_nveccuda_shared SHARED ${nveccuda_SOURCES} ${shared_SOURCES}
                   OPTIONS -DSUNDIALS_MPI_ENABLED=0)
  TARGET_COMPILE_DEFINITIONS(sundials_nveccuda_shared PUBLIC -DSUNDIALS_MPI_ENABLED=0)
  SET_TARGET_PROPERTIES(sundials_nveccuda_shared
    PROPERTIES OUTPUT_NAME sundials_nveccuda CLEAN_DIRECT_OUTPUT 1
               VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nveccuda_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- MPI+CUDA
  IF(MPI_ENABLE)
    CUDA_ADD_LIBRARY(sundials_nvecmpicuda_shared SHARED ${nveccuda_SOURCES} ${shared_SOURCES})
    SET_TARGET_PROPERTIES(sundials_nvecmpicuda_shared
      PROPERTIES OUTPUT_NAME sundials_nvecmpicuda CLEAN_DIRECT_OUTPUT 1
                 VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_nvecmpicuda_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF()
ENDIF(BUILD_SHARED_LIBS)

# Install the CUDA NVector header files
INSTALL(FILES ${nveccuda_HEADERS} DESTINATION include/nvector)
INSTALL(DIRECTORY ${sundials_SOURCE_DIR}/include/nvector/cuda DESTINATION include/nvector)

MESSAGE(STATUS "Added NVECTOR_CUDA module")
