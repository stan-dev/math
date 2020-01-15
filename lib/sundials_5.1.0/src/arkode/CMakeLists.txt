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
# CMakeLists.txt file for the ARKODE library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall ARKODE\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable arkode_SOURCES with the sources for the ARKODE library
set(arkode_SOURCES
  arkode.c
  arkode_adapt.c
  arkode_arkstep.c
  arkode_arkstep_io.c
  arkode_arkstep_nls.c
  arkode_bandpre.c
  arkode_bbdpre.c
  arkode_butcher.c
  arkode_butcher_dirk.c
  arkode_butcher_erk.c
  arkode_erkstep.c
  arkode_erkstep_io.c
  arkode_interp.c
  arkode_io.c
  arkode_ls.c
  arkode_mristep.c
  arkode_mristep_io.c
  arkode_root.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the ARKODE library
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
# also be included in the ARKODE library
set(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/sunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/sunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/sunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the ARKODE library
set(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol/band/sunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/dense/sunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spbcgs/sunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spfgmr/sunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spgmr/sunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/pcg/sunlinsol_pcg.c
  )

# Add variable sunnonlinsol_SOURCES with the common SUNNonlinearSolver sources which will
# also be included in the ARKODE library
set(sunnonlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/newton/sunnonlinsol_newton.c
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c
  )

# Add variable arkode_HEADERS with the exported ARKODE header files
set(arkode_HEADERS
  arkode.h
  arkode_arkstep.h
  arkode_bandpre.h
  arkode_bbdpre.h
  arkode_butcher.h
  arkode_butcher_erk.h
  arkode_butcher_dirk.h
  arkode_erkstep.h
  arkode_ls.h
  arkode_mristep.h
  )

# Add prefix with complete path to the ARKODE header files
add_prefix(${sundials_SOURCE_DIR}/include/arkode/ arkode_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
include_directories(.)
include_directories(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
if(BUILD_STATIC_LIBS)

  # Add the build target for the static ARKODE library
  add_library(sundials_arkode_static STATIC
    ${arkode_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES} ${sunnonlinsol_SOURCES})

  # Set the library name and make sure it is not deleted
  set_target_properties(sundials_arkode_static
    PROPERTIES OUTPUT_NAME sundials_arkode CLEAN_DIRECT_OUTPUT 1)

  # Install the ARKODE library
  install(TARGETS sundials_arkode_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Build the shared library
if(BUILD_SHARED_LIBS)

  # Add the build target for the ARKODE library
  add_library(sundials_arkode_shared SHARED
    ${arkode_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES} ${sunnonlinsol_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_arkode_shared m)
  endif()

  # Set the library name and make sure it is not deleted
  set_target_properties(sundials_arkode_shared
    PROPERTIES OUTPUT_NAME sundials_arkode CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  set_target_properties(sundials_arkode_shared
    PROPERTIES VERSION ${arkodelib_VERSION} SOVERSION ${arkodelib_SOVERSION})

  # Install the ARKODE library
  install(TARGETS sundials_arkode_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the ARKODE header files
install(FILES ${arkode_HEADERS} DESTINATION include/arkode)

#
message(STATUS "Added ARKODE module")
