# ---------------------------------------------------------------
# Programmer:  Slaven Peles @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the cuda NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_RAJA\n\")")

# Add variable nveccuda_SOURCES with the sources for the NVECRAJA lib
SET(nvecraja_SOURCES 
  nvector_raja.cu
)

# Tell compiler it is a CUDA source
set_source_files_properties(${nvecraja_SOURCES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECRAJA library
SET(shared_SOURCES 
  sundials_math.c
)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nveccuda_HEADERS with the exported NVECRAJA header files
SET(nvecraja_HEADERS
  nvector_raja.h 
#  raja/Vector.hpp 
#  raja/VectorKernels.cuh
)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecraja_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECRAJA library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECRAJA library
IF(BUILD_STATIC_LIBS)
  CUDA_ADD_LIBRARY(sundials_nvecraja_static STATIC ${nvecraja_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecraja_static
    PROPERTIES OUTPUT_NAME sundials_nvecraja CLEAN_DIRECT_OUTPUT 1
  )
  INSTALL(TARGETS sundials_nvecraja_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECRAJA library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECRAJA library
IF(BUILD_SHARED_LIBS)
  CUDA_ADD_LIBRARY(sundials_nvecraja_shared SHARED ${nvecraja_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecraja_shared
    PROPERTIES OUTPUT_NAME sundials_nvecraja CLEAN_DIRECT_OUTPUT 1
  )
  SET_TARGET_PROPERTIES(sundials_nvecraja_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION}
  )
#  SET_TARGET_PROPERTIES(sundials_nvecraja_shared 
#    PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE
#  )
  INSTALL(TARGETS sundials_nvecraja_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECRAJA header files
INSTALL(FILES ${nvecraja_HEADERS} DESTINATION include/nvector)
INSTALL(DIRECTORY ${sundials_SOURCE_DIR}/include/nvector/raja DESTINATION include/nvector)


MESSAGE(STATUS "Added NVECTOR_RAJA module")
