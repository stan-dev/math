# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Macro that checks if the requested CUDA arch of the format sm_[0-9][0-9]
# is greater or equal to the CUDA_ARCH currently being used.
# 
# Example: sundials_cuda_arch_geq(arch_ok "sm_60")
# -----------------------------------------------------------------------------

macro(sundials_cuda_arch_geq OUTPUT_VAR REQUESTED_ARCH)
  set(options )
  set(oneValueArgs REQUESTED_ARCH)
  set(multiValueArgs )

  cmake_parse_arguments(sundials_cuda_arch_geq "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  string(REGEX REPLACE "sm_([0-9][0-9])" "\\1" _capability_number ${CUDA_ARCH})
  string(REGEX REPLACE "sm_([0-9][0-9])" "\\1" _reqd_capability_number ${REQUESTED_ARCH})
  if (_capability_number GREATER _reqd_capability_number)
    set(${OUTPUT_VAR} TRUE)
  elseif (_capability_number EQUAL _reqd_capability_number)
    set(${OUTPUT_VAR} TRUE)
  else()
    set(${OUTPUT_VAR} FALSE)
  endif()
endmacro()
