# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Cody J. Balos @ LLNL
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
# C++-related tests for SUNDIALS CMake-based configuration.
# ---------------------------------------------------------------

enable_language(CXX)
set(CXX_FOUND TRUE)

# ---------------------------------------------------------------
# Option to specify the C++ standard SUNDIALS will use. Defined
# here so it is set in the same configuration pass as the C++
# compiler and related options.
# ---------------------------------------------------------------

set(DOCSTR "The C++ standard to use if C++ is enabled (98, 11, 14, 17, 20)")
sundials_option(CMAKE_CXX_STANDARD STRING "${DOCSTR}" "11"
                OPTIONS "98;11;14;17;20")

message(STATUS "CXX standard set to ${CMAKE_CXX_STANDARD}")
