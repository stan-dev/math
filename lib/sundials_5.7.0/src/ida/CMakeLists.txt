# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Radu Serban, Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
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

# Add variable ida_HEADERS with the exported IDA header files
set(ida_HEADERS
  ida.h
  ida_bbdpre.h
  ida_direct.h
  ida_ls.h
  ida_spils.h
  )

# Add prefix with complete path to the IDA header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/ida/ ida_HEADERS)

# Create the library
sundials_add_library(sundials_ida
  SOURCES
    ${ida_SOURCES}
  HEADERS
    ${ida_HEADERS}
  INCLUDE_SUBDIR
    ida
  OBJECT_LIBRARIES
    sundials_generic_obj
    sundials_nvecserial_obj
    sundials_sunmatrixband_obj
    sundials_sunmatrixdense_obj
    sundials_sunmatrixsparse_obj
    sundials_sunlinsolband_obj
    sundials_sunlinsoldense_obj
    sundials_sunlinsolspbcgs_obj
    sundials_sunlinsolspfgmr_obj
    sundials_sunlinsolspgmr_obj
    sundials_sunlinsolsptfqmr_obj
    sundials_sunlinsolpcg_obj
    sundials_sunnonlinsolnewton_obj
    sundials_sunnonlinsolfixedpoint_obj
  OUTPUT_NAME
    sundials_ida
  VERSION
    ${idalib_VERSION}
  SOVERSION
    ${idalib_SOVERSION}
)

# Finished IDA
message(STATUS "Added IDA module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

if(BUILD_FORTRAN77_INTERFACE)
  add_subdirectory(fcmix)
endif()
