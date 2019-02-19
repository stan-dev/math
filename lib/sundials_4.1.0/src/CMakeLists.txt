# ---------------------------------------------------------------
# Programmer:  David Gardner, Slaven Peles, and Cody Balos @ LLNL
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
# src level CMakeLists.txt for SUNDIALS (for cmake build system)
# ---------------------------------------------------------------

# Add all of the shared SUNDIALS components
ADD_SUBDIRECTORY(sundials)
ADD_SUBDIRECTORY(nvector)
ADD_SUBDIRECTORY(sunmatrix)
ADD_SUBDIRECTORY(sunlinsol)
ADD_SUBDIRECTORY(sunnonlinsol)

# ARKODE library
IF(BUILD_ARKODE)
  ADD_SUBDIRECTORY(arkode)
  IF(F77_INTERFACE_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(arkode/fcmix)
  ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)
ENDIF(BUILD_ARKODE)

# CVODE library
IF(BUILD_CVODE)
  ADD_SUBDIRECTORY(cvode)
  IF(F77_INTERFACE_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(cvode/fcmix)
  ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)
ENDIF(BUILD_CVODE)

# CVODES library
IF(BUILD_CVODES)
  ADD_SUBDIRECTORY(cvodes)
ENDIF(BUILD_CVODES)

# IDA library
IF(BUILD_IDA)
  ADD_SUBDIRECTORY(ida)
  IF(F77_INTERFACE_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(ida/fcmix)
  ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)
ENDIF(BUILD_IDA)

# IDAS library
IF(BUILD_IDAS)
  ADD_SUBDIRECTORY(idas)
ENDIF(BUILD_IDAS)

# KINSOL library
IF(BUILD_KINSOL)
  ADD_SUBDIRECTORY(kinsol)
  IF(F77_INTERFACE_ENABLE AND F77_FOUND)
    ADD_SUBDIRECTORY(kinsol/fcmix)
  ENDIF(F77_INTERFACE_ENABLE AND F77_FOUND)
ENDIF(BUILD_KINSOL)

# CPODES library
IF(BUILD_CPODES)
  ADD_SUBDIRECTORY(cpodes)
ENDIF(BUILD_CPODES)

