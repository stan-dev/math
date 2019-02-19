# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# examples/sunnonlinsol level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# Always add the Newton and fixed point examples
add_subdirectory(newton)
add_subdirectory(fixedpoint)
