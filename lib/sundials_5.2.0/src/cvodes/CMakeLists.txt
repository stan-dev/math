# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Radu Serban @ LLNL
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
# CMakeLists.txt file for the CVODES library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall CVODES\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable cvodes_SOURCES with the sources for the CVODES library
set(cvodes_SOURCES
  cvodea.c
  cvodea_io.c
  cvodes.c
  cvodes_bandpre.c
  cvodes_bbdpre.c
  cvodes_diag.c
  cvodes_direct.c
  cvodes_io.c
  cvodes_ls.c
  cvodes_nls.c
  cvodes_nls_sim.c
  cvodes_nls_stg.c
  cvodes_nls_stg1.c
  cvodes_spils.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CVODES library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nonlinearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_band.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_dense.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_direct.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_version.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector_senswrapper.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_futils.c
  ${sundials_SOURCE_DIR}/src/nvector/serial/nvector_serial.c
  )

# Add variable sunmatrix_SOURCES with the common SUNMatrix sources which will
# also be included in the CVODES library
set(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/sunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/sunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/sunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the CVODES library
set(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol/band/sunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/dense/sunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spbcgs/sunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spfgmr/sunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spgmr/sunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/pcg/sunlinsol_pcg.c
  )

# Add variable sunnonlinsol_SOURCES with the common SUNNonlinearSolver sources
# which will also be included in the CVODES library
set(sunnonlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/newton/sunnonlinsol_newton.c
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c
  )

# Add variable cvodes_HEADERS with the exported CVODES header files
set(cvodes_HEADERS
  cvodes.h
  cvodes_bandpre.h
  cvodes_bbdpre.h
  cvodes_diag.h
  cvodes_direct.h
  cvodes_ls.h
  cvodes_spils.h
  )

# Add prefix with complete path to the CVODES header files
add_prefix(${sundials_SOURCE_DIR}/include/cvodes/ cvodes_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
include_directories(.)
include_directories(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
if(BUILD_STATIC_LIBS)

  # Add the build target for the static CVODES library
  add_library(sundials_cvodes_static STATIC
    ${cvodes_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES}
    ${sunnonlinsol_SOURCES})

  # Set the library name and make sure it is not deleted
  set_target_properties(sundials_cvodes_static
    PROPERTIES OUTPUT_NAME sundials_cvodes CLEAN_DIRECT_OUTPUT 1)

  # Install the CVODES library
  install(TARGETS sundials_cvodes_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Build the shared library
if(BUILD_SHARED_LIBS)

  # Add the build target for the CVODES library
  add_library(sundials_cvodes_shared SHARED
    ${cvodes_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES}
    ${sunnonlinsol_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_cvodes_shared m)
  endif()

  # Set the library name and make sure it is not deleted
  set_target_properties(sundials_cvodes_shared
    PROPERTIES OUTPUT_NAME sundials_cvodes CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  set_target_properties(sundials_cvodes_shared
    PROPERTIES VERSION ${cvodeslib_VERSION} SOVERSION ${cvodeslib_SOVERSION})

  # Install the CVODES library
  install(TARGETS sundials_cvodes_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the CVODES header files
install(FILES ${cvodes_HEADERS} DESTINATION include/cvodes)

#
message(STATUS "Added CVODES module")
