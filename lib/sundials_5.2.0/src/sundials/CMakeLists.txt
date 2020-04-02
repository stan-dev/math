# ---------------------------------------------------------------
# Programmer(s): Radu Serban and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the generic SUNDIALS modules
# ---------------------------------------------------------------

# From here we only install the generic SUNDIALS headers.
# The implementations themselves are incorporated in the individual SUNDIALS solver libraries.

install(CODE "MESSAGE(\"\nInstall shared components\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable sundials_HEADERS with the exported SUNDIALS header files
set(sundials_HEADERS
  sundials_band.h
  sundials_dense.h
  sundials_direct.h
  sundials_fnvector.h
  sundials_futils.h
  sundials_iterative.h
  sundials_linearsolver.h
  sundials_math.h
  sundials_matrix.h
  sundials_nonlinearsolver.h
  sundials_mpi_types.h
  sundials_nvector.h
  sundials_types.h
  sundials_version.h
  )

# Add prefix with complete path to the SUNDIALS header files
add_prefix(${sundials_SOURCE_DIR}/include/sundials/ sundials_HEADERS)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Create an object library out of the generic sundials modules
if(BUILD_STATIC_LIBS)
  add_library(sundials_generic_static_obj OBJECT
    sundials_band.c
    sundials_dense.c
    sundials_direct.c
    sundials_futils.c
    sundials_iterative.c
    sundials_linearsolver.c
    sundials_math.c
    sundials_matrix.c
    sundials_nonlinearsolver.c
    sundials_nvector.c
    sundials_nvector_senswrapper.c
    sundials_version.c)
endif()

if(BUILD_SHARED_LIBS)
  add_library(sundials_generic_shared_obj OBJECT
      sundials_band.c
      sundials_dense.c
      sundials_direct.c
      sundials_futils.c
      sundials_iterative.c
      sundials_linearsolver.c
      sundials_math.c
      sundials_matrix.c
      sundials_nonlinearsolver.c
      sundials_nvector.c
      sundials_nvector_senswrapper.c
      sundials_version.c)
  set_target_properties(sundials_generic_shared_obj PROPERTIES
                        POSITION_INDEPENDENT_CODE TRUE)
endif()

# Install the SUNDIALS header files
install(FILES ${sundials_HEADERS} DESTINATION include/sundials)

# If Blas/Lapack support was enabled, install the Lapack interface headers
if(LAPACK_FOUND)
  set(sundials_BL_HEADERS sundials_lapack.h)
  add_prefix(${sundials_SOURCE_DIR}/include/sundials/ sundials_BL_HEADERS)
  install(FILES ${sundials_BL_HEADERS} DESTINATION include/sundials)
endif(LAPACK_FOUND)
