# ---------------------------------------------------------------
# Programmer:  David Gardner and Slaven Peles @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# examples/sunmatrix level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# Always add the serial sunmatrix dense/band/sparse examples
ADD_SUBDIRECTORY(dense)
ADD_SUBDIRECTORY(band)
ADD_SUBDIRECTORY(sparse)

