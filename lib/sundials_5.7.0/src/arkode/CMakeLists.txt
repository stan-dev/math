# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the ARKODE library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall ARKODE\n\")")

# Add variable arkode_SOURCES with the sources for the ARKODE library
set(arkode_SOURCES
  arkode_adapt.c
  arkode_arkstep_io.c
  arkode_arkstep_nls.c
  arkode_arkstep.c
  arkode_bandpre.c
  arkode_bbdpre.c
  arkode_butcher_dirk.c
  arkode_butcher_erk.c
  arkode_butcher.c
  arkode_erkstep_io.c
  arkode_erkstep.c
  arkode_interp.c
  arkode_io.c
  arkode_ls.c
  arkode_mri_tables.c
  arkode_mristep_io.c
  arkode_mristep_nls.c
  arkode_mristep.c
  arkode_root.c
  arkode.c
)

# Add variable arkode_HEADERS with the exported ARKODE header files
set(arkode_HEADERS
  arkode_arkstep.h
  arkode_bandpre.h
  arkode_bbdpre.h
  arkode_butcher_dirk.h
  arkode_butcher_erk.h
  arkode_butcher.h
  arkode_erkstep.h
  arkode_ls.h
  arkode_mristep.h
  arkode.h
)

# Add prefix with complete path to the ARKODE header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/arkode/ arkode_HEADERS)

# Create the sundials_arkode library
sundials_add_library(sundials_arkode
  SOURCES
    ${arkode_SOURCES}
  HEADERS
    ${arkode_HEADERS}
  INCLUDE_SUBDIR
    arkode
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
    sundials_arkode
  VERSION
    ${arkodelib_VERSION}
  SOVERSION
    ${arkodelib_SOVERSION}
)

# Finished ARKODE
message(STATUS "Added ARKODE module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()

# Add FCMIX code if the F77 interface is enabled
if(F77_FOUND AND BUILD_FORTRAN77_INTERFACE)
  add_subdirectory(fcmix)
endif()

# Add ARKODE XBraid interface
if(ENABLE_XBRAID)
  add_subdirectory(xbraid)
endif()
