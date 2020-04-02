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
# CMakeLists.txt file for the FIDA library
# ---------------------------------------------------------------

# Add variable fida_SOURCES with the sources for the FIDA library
set(fida_SOURCES
  fidaband.c
  fidabbd.c
  fida.c
  fidadense.c
  fidaewt.c
  fidajtimes.c
  fidanullmatrix.c
  fidanullnonlinsol.c
  fidapreco.c
  fidaroot.c
  fidasparse.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the IDA library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/nvector/serial/fnvector_serial.c
  )

# Add variable sunmatrix_SOURCES with the common SUNMatrix sources which will
# also be included in the IDA library
set(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/fsunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/fsunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/fsunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the IDA library
set(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol/band/fsunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/dense/fsunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spbcgs/fsunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spfgmr/fsunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spgmr/fsunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/sptfqmr/fsunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/pcg/fsunlinsol_pcg.c
  )

# Add variable sunnonlinsol_SOURCES with the common SUNNonlinearSolver sources
# which will also be included in the IDA library
set(sunnonlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/newton/fsunnonlinsol_newton.c
  )

# Add source directories to include directories for access to
# implementation only header files (both for fida and ida)
include_directories(.)
include_directories(..)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Only build STATIC libraries (we cannot build shared libraries
# for the FCMIX interfaces due to unresolved symbol errors
# coming from inexistent user-provided functions)

# Add the build target for the FIDA library
add_library(sundials_fida_static STATIC
  ${fida_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES}
  ${sunnonlinsol_SOURCES})

# Set the library name and make sure it is not deleted
set_target_properties(sundials_fida_static
  PROPERTIES OUTPUT_NAME sundials_fida CLEAN_DIRECT_OUTPUT 1)

# Install the FIDA library
install(TARGETS sundials_fida_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
#
message(STATUS "Added IDA FCMIX module")
