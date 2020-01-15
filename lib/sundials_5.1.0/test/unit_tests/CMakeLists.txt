# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# unit_tests level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# The unit tests will be linked with the library corresponding to the build
# type. Whenever building shared libraries, use them to link the tests.
if(BUILD_SHARED_LIBS)
  set(LINK_LIBRARY_TYPE "shared")
else()
  set(LINK_LIBRARY_TYPE "static")
endif()

# Add ARKode unit tests
if(BUILD_ARKODE)
  add_subdirectory(arkode)
endif()
