# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# This module handles setting the index sizes in SUNDIALS.
# The SUNDIALS_INDEX_SIZE build option can be set to 64 or 32
# for 64-bit or 32-bit signed integer indices. The module will try
# various types for each size until it finds one that matches the
# desired size. If no match is found, a user can provide a type via
# the advanced options SUNDIALS_CINDEX_TYPE and SUNDIALS_FINDEX_TYPE.
# -----------------------------------------------------------------------------

include(CheckTypeSize)

if(SUNDIALS_INDEX_SIZE MATCHES "64")
  set(SUNDIALS_CINDEX_TYPE "")

  # if the user specified an index type use it, otherwise try the standard options
  if (SUNDIALS_INDEX_TYPE)
    set(POSSIBLE_INT64 ${SUNDIALS_INDEX_TYPE})
  else()
    set(POSSIBLE_INT64 int64_t;__int64;long long;long)
  endif()

  foreach(INT64_TYPE ${POSSIBLE_INT64})
    string(REPLACE " " "_" INT64_TYPE_NOSPACE ${INT64_TYPE})
    check_type_size("${INT64_TYPE}" HAS_${INT64_TYPE_NOSPACE})
    if(HAS_${INT64_TYPE_NOSPACE} EQUAL "8")
      set(SUNDIALS_CINDEX_TYPE ${INT64_TYPE})
      message(STATUS "Using ${INT64_TYPE} for indices")
      break()
    endif()
  endforeach()

  if(NOT SUNDIALS_CINDEX_TYPE)
    print_error("No integer type of size 8 was found.\n\
                 Tried ${POSSIBLE_INT64}.\n\
                 Try setting the advanced option SUNDIALS_INDEX_TYPE.")
  endif()

  # set Fortran integer size too
  set(SUNDIALS_FINDEX_TYPE "8")
elseif(SUNDIALS_INDEX_SIZE MATCHES "32")
  set(SUNDIALS_CINDEX_TYPE "")

  # if the user specified an index type use it, otherwise try the standard options
  if (SUNDIALS_INDEX_TYPE)
    set(POSSIBLE_INT32 ${SUNDIALS_INDEX_TYPE})
  else()
    set(POSSIBLE_INT32 int32_t;int;long)
  endif()

  foreach(INT32_TYPE ${POSSIBLE_INT32})
    string(REPLACE " " "_" INT32_TYPE_NOSPACE ${INT32_TYPE})
    check_type_size("${INT32_TYPE}" HAS_${INT32_TYPE_NOSPACE})
    if(HAS_${INT32_TYPE_NOSPACE} EQUAL "4")
      set(SUNDIALS_CINDEX_TYPE ${INT32_TYPE})
      message(STATUS "Using ${INT32_TYPE} for indices")
      break()
    endif()
  endforeach()

  if(NOT SUNDIALS_CINDEX_TYPE)
    print_error("No integer type of size 4 was found.\n\
                 Tried ${POSSIBLE_INT32}\n\
                 Try setting the advanced option SUNDIALS_INDEX_TYPE.")
  endif()

  # set Fortran integer size too
  set(SUNDIALS_FINDEX_TYPE "4")
else()
  print_error("Invalid index size.")
endif()
