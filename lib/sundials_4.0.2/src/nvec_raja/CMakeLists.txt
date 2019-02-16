# ---------------------------------------------------------------
# Programmer:  Slaven Peles @ LLNL
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
# CMakeLists.txt file for the RAJA NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_RAJA\n\")")

# Add variable nvecraja_SOURCES with the sources for the NVECRAJA lib
SET(nvecraja_SOURCES
  nvector_raja.cu
)

# Tell compiler it is a CUDA source
set_source_files_properties(${nvecraja_SOURCES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

# Set C++ compiler flags to include RAJA flags
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${RAJA_COMPILE_FLAGS})

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECRAJA library
SET(shared_SOURCES
  sundials_math.c
  sundials_mpi.c
)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nveccuda_HEADERS with the exported NVECRAJA header files
SET(nvecraja_HEADERS nvector_raja.h)
IF (MPI_ENABLE)
  LIST(APPEND nvecraja_HEADERS nvector_mpiraja.h)
ENDIF()
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecraja_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ${MPI_CXX_INCLUDE_PATH})

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECRAJA library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECRAJA library
IF(BUILD_STATIC_LIBS)
  # ----------------------- CUDA only
  # The FindCUDA module does not properly forward compile options using target_* commands.
  # So as long as we use it (required in CMake < 3.8), we have to manually add compile
  # options in the CUDA_ADD_* commands.
  CUDA_ADD_LIBRARY(sundials_nveccudaraja_static STATIC
    ${nvecraja_SOURCES} ${shared_SOURCES}
    OPTIONS -DSUNDIALS_MPI_ENABLED=0
  )
  TARGET_COMPILE_DEFINITIONS(sundials_nveccudaraja_static PUBLIC -DSUNDIALS_MPI_ENABLED=0)
  SET_TARGET_PROPERTIES(sundials_nveccudaraja_static
    PROPERTIES OUTPUT_NAME sundials_nveccudaraja CLEAN_DIRECT_OUTPUT 1
  )
  INSTALL(TARGETS sundials_nveccudaraja_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  # MPI+CUDA
  IF(MPI_ENABLE)
    CUDA_ADD_LIBRARY(sundials_nvecmpicudaraja_static STATIC ${nvecraja_SOURCES} ${shared_SOURCES})
    SET_TARGET_PROPERTIES(sundials_nvecmpicudaraja_static
      PROPERTIES OUTPUT_NAME sundials_nvecmpicudaraja CLEAN_DIRECT_OUTPUT 1
    )
    INSTALL(TARGETS sundials_nvecmpicudaraja_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF()
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECRAJA library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECRAJA library
IF(BUILD_SHARED_LIBS)
  # ----------------------- CUDA only
  CUDA_ADD_LIBRARY(sundials_nveccudaraja_shared SHARED
    ${nvecraja_SOURCES} ${shared_SOURCES}
    OPTIONS -DSUNDIALS_MPI_ENABLED=0
  )
  TARGET_COMPILE_DEFINITIONS(sundials_nveccudaraja_shared PUBLIC -DSUNDIALS_MPI_ENABLED=0)
  SET_TARGET_PROPERTIES(sundials_nveccudaraja_shared
    PROPERTIES OUTPUT_NAME sundials_nveccudaraja CLEAN_DIRECT_OUTPUT 1
               VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION}
  )
  INSTALL(TARGETS sundials_nveccudaraja_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- MPI+CUDA
  IF(MPI_ENABLE)
    CUDA_ADD_LIBRARY(sundials_nvecmpicudaraja_shared SHARED ${nvecraja_SOURCES} ${shared_SOURCES})
    SET_TARGET_PROPERTIES(sundials_nvecmpicudaraja_shared
      PROPERTIES OUTPUT_NAME sundials_nvecmpicudaraja CLEAN_DIRECT_OUTPUT 1
                 VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION}
    )
    INSTALL(TARGETS sundials_nvecmpicudaraja_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF()
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECRAJA header files
INSTALL(FILES ${nvecraja_HEADERS} DESTINATION include/nvector)
INSTALL(DIRECTORY ${sundials_SOURCE_DIR}/include/nvector/raja DESTINATION include/nvector)


MESSAGE(STATUS "Added NVECTOR_RAJA module")
