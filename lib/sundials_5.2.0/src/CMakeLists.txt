# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# src level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Add all of the shared SUNDIALS components
add_subdirectory(sundials)
add_subdirectory(nvector)
add_subdirectory(sunmatrix)
add_subdirectory(sunlinsol)
add_subdirectory(sunnonlinsol)

# ARKODE library
if(BUILD_ARKODE)
  add_subdirectory(arkode)
  if(F77_INTERFACE_ENABLE AND F77_FOUND)
    add_subdirectory(arkode/fcmix)
  endif(F77_INTERFACE_ENABLE AND F77_FOUND)
endif(BUILD_ARKODE)

# CVODE library
if(BUILD_CVODE)
  add_subdirectory(cvode)
  if(F77_INTERFACE_ENABLE AND F77_FOUND)
    add_subdirectory(cvode/fcmix)
  endif(F77_INTERFACE_ENABLE AND F77_FOUND)
endif(BUILD_CVODE)

# CVODES library
if(BUILD_CVODES)
  add_subdirectory(cvodes)
endif(BUILD_CVODES)

# IDA library
if(BUILD_IDA)
  add_subdirectory(ida)
  if(F77_INTERFACE_ENABLE AND F77_FOUND)
    add_subdirectory(ida/fcmix)
  endif(F77_INTERFACE_ENABLE AND F77_FOUND)
endif(BUILD_IDA)

# IDAS library
if(BUILD_IDAS)
  add_subdirectory(idas)
endif(BUILD_IDAS)

# KINSOL library
if(BUILD_KINSOL)
  add_subdirectory(kinsol)
  if(F77_INTERFACE_ENABLE AND F77_FOUND)
    add_subdirectory(kinsol/fcmix)
  endif(F77_INTERFACE_ENABLE AND F77_FOUND)
endif(BUILD_KINSOL)

# CPODES library
if(BUILD_CPODES)
  add_subdirectory(cpodes)
endif(BUILD_CPODES)
