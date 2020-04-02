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
# CMakeLists.txt file for the IDA library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall IDA\n\")")

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)

# Add variable ida_SOURCES with the sources for the IDA library
set(ida_SOURCES
  ida.c
  ida_bbdpre.c
  ida_direct.c
  ida_ic.c
  ida_io.c
  ida_ls.c
  ida_nls.c
  ida_spils.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the IDA library
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
# also be included in the IDA library
set(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/sunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/sunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/sunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the IDA library
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
# which will also be included in the IDA library
set(sunnonlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/newton/sunnonlinsol_newton.c
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c
  )

# Add variable ida_HEADERS with the exported IDA header files
set(ida_HEADERS
  ida.h
  ida_bbdpre.h
  ida_direct.h
  ida_ls.h
  ida_spils.h
  )

# Add prefix with complete path to the IDA header files
add_prefix(${sundials_SOURCE_DIR}/include/ida/ ida_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
include_directories(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

# Build the static library
if(BUILD_STATIC_LIBS)

  # Add the build target for the static IDA library
  add_library(sundials_ida_static STATIC
    ${ida_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES}
    ${sunnonlinsol_SOURCES})

  # Set the library name and make sure it is not deleted
  set_target_properties(sundials_ida_static
    PROPERTIES OUTPUT_NAME sundials_ida CLEAN_DIRECT_OUTPUT 1)

  # Install the IDA library
  install(TARGETS sundials_ida_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Build the shared library
if(BUILD_SHARED_LIBS)

  # Add the build target for the IDA library
  add_library(sundials_ida_shared SHARED
    ${ida_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES}
    ${sunnonlinsol_SOURCES})

  if(UNIX)
    target_link_libraries(sundials_ida_shared m)
  endif()

  # Set the library name and make sure it is not deleted
  set_target_properties(sundials_ida_shared
    PROPERTIES OUTPUT_NAME sundials_ida CLEAN_DIRECT_OUTPUT 1)

  # Set VERSION and SOVERSION for shared libraries
  set_target_properties(sundials_ida_shared
    PROPERTIES VERSION ${idalib_VERSION} SOVERSION ${idalib_SOVERSION})

  # Install the IDA library
  install(TARGETS sundials_ida_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the IDA header files
install(FILES ${ida_HEADERS} DESTINATION include/ida)

#
message(STATUS "Added IDA module")
