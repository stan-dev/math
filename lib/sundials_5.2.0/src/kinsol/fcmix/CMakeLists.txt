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
# CMakeLists.txt file for the FKINSOL library
# ---------------------------------------------------------------

# Add variable fcvode_SOURCES with the sources for the FCVODE library
set(fkinsol_SOURCES
  fkinband.c
  fkinbbd.c
  fkindense.c
  fkinjtimes.c
  fkinnullmatrix.c
  fkinnulllinsol.c
  fkinpreco.c
  fkinsol.c
  fkinsparse.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the ARKODE library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/nvector/serial/fnvector_serial.c
  )

# Add variable sunmatrix_SOURCES with the common SUNMatrix sources which will
# also be included in the ARKODE library
set(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/fsunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/fsunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/fsunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the ARKODE library
set(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol/band/fsunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/dense/fsunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spbcgs/fsunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spfgmr/fsunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spgmr/fsunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/sptfqmr/fsunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/pcg/fsunlinsol_pcg.c
  )

# Add source directories to include directories for access to
# implementation only header files (both for fkinsol and kinsol)
include_directories(.)
include_directories(..)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Only build STATIC libraries (we cannot build shared libraries
# for the FCMIX interfaces due to unresolved symbol errors
# coming from inexistent user-provided functions)

# Add the build target for the FKINSOL library
add_library(sundials_fkinsol_static STATIC
  ${fkinsol_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES} )

# Set the library name and make sure it is not deleted
set_target_properties(sundials_fkinsol_static
  PROPERTIES OUTPUT_NAME sundials_fkinsol CLEAN_DIRECT_OUTPUT 1)

# Install the FKINSOL library
install(TARGETS sundials_fkinsol_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

#
message(STATUS "Added KINSOL FCMIX module")
