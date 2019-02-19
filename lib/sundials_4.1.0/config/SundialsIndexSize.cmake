# -----------------------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
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

INCLUDE(CheckTypeSize)

IF(SUNDIALS_INDEX_SIZE MATCHES "64")
  SET(SUNDIALS_CINDEX_TYPE "")
  
  # if the user specified an index type use it, otherwise try the standard options
  IF (SUNDIALS_INDEX_TYPE)
    SET(POSSIBLE_INT64 ${SUNDIALS_INDEX_TYPE})
  ELSE()
    SET(POSSIBLE_INT64 int64_t;__int64;long long;long)
  ENDIF()
  
  FOREACH(INT64_TYPE ${POSSIBLE_INT64})
    string(REPLACE " " "_" INT64_TYPE_NOSPACE ${INT64_TYPE})
    CHECK_TYPE_SIZE("${INT64_TYPE}" HAS_${INT64_TYPE_NOSPACE})
    IF(HAS_${INT64_TYPE_NOSPACE} EQUAL "8")
      SET(SUNDIALS_CINDEX_TYPE ${INT64_TYPE})
      MESSAGE(STATUS "Using ${INT64_TYPE} for indices")
      BREAK()
    ENDIF()
  ENDFOREACH()
  
  IF(NOT SUNDIALS_CINDEX_TYPE)
    PRINT_ERROR("No integer type of size 8 was found.\n\
                 Tried ${POSSIBLE_INT64}.\n\
                 Try setting the advanced option SUNDIALS_INDEX_TYPE.")
  ENDIF()

  # set Fortran integer size too
  SET(SUNDIALS_FINDEX_TYPE "8")
  # prepare substitution variable INDEX_TYPE for sundials_config.h
  SET(INDEX_TYPE "#define SUNDIALS_INT${SUNDIALS_INDEX_SIZE}_T 1")
ELSEIF(SUNDIALS_INDEX_SIZE MATCHES "32")
  SET(SUNDIALS_CINDEX_TYPE "")
  
  # if the user specified an index type use it, otherwise try the standard options
  IF (SUNDIALS_INDEX_TYPE)
    SET(POSSIBLE_INT32 ${SUNDIALS_INDEX_TYPE})
  ELSE()
    SET(POSSIBLE_INT32 int32_t;int;long)
  ENDIF()
  
  FOREACH(INT32_TYPE ${POSSIBLE_INT32})
    string(REPLACE " " "_" INT32_TYPE_NOSPACE ${INT32_TYPE})
    CHECK_TYPE_SIZE("${INT32_TYPE}" HAS_${INT32_TYPE_NOSPACE})
    IF(HAS_${INT32_TYPE_NOSPACE} EQUAL "4")
      SET(SUNDIALS_CINDEX_TYPE ${INT32_TYPE})
      MESSAGE(STATUS "Using ${INT32_TYPE} for indices")
      BREAK()
    ENDIF()
  ENDFOREACH()
  
  IF(NOT SUNDIALS_CINDEX_TYPE)
    PRINT_ERROR("No integer type of size 4 was found.\n\
                 Tried ${POSSIBLE_INT32}\n\
                 Try setting the advanced option SUNDIALS_INDEX_TYPE.")
  ENDIF()
  
  # set Fortran integer size too
  SET(SUNDIALS_FINDEX_TYPE "4")
  # prepare substitution variable INDEX_TYPE for sundials_config.h
  SET(INDEX_TYPE "#define SUNDIALS_INT${SUNDIALS_INDEX_SIZE}_T 1")
ELSE()
  PRINT_ERROR("Invalid index size.")
ENDIF()
