# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# matrix level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# Always add SUNDIALS provided matrix modules
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(sparse)

