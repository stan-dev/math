# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
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
add_subdirectory(sunmemory)

# ARKODE library
if(BUILD_ARKODE)
  add_subdirectory(arkode)
endif()

# CVODE library
if(BUILD_CVODE)
  add_subdirectory(cvode)
endif()

# CVODES library
if(BUILD_CVODES)
  add_subdirectory(cvodes)
endif()

# IDA library
if(BUILD_IDA)
  add_subdirectory(ida)
endif()

# IDAS library
if(BUILD_IDAS)
  add_subdirectory(idas)
endif()

# KINSOL library
if(BUILD_KINSOL)
  add_subdirectory(kinsol)
endif()

# CPODES library
if(BUILD_CPODES)
  add_subdirectory(cpodes)
endif()
