# ---------------------------------------------------------------
# Programmer:  David J. Gardner @ LLNL
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
# Print warning is the user sets a deprecated CMake variable and
# copy the value into the correct CMake variable
# ---------------------------------------------------------------

# macro to print warning for deprecated CMake variable
MACRO(PRINT_DEPRECATED old_variable new_variable)
  PRINT_WARNING("${old_variable} is deprecated and will be removed in the future."
                "Copying value to ${new_variable}.")
ENDMACRO()

IF(DEFINED EXAMPLES_ENABLE)
  PRINT_DEPRECATED(EXAMPLES_ENABLE EXAMPLES_ENABLE_C)
  FORCE_VARIABLE(EXAMPLES_ENABLE_C BOOL "Build SUNDIALS C examples" ${EXAMPLES_ENABLE})
  UNSET(EXAMPLES_ENABLE CACHE)
ENDIF()

IF(DEFINED CXX_ENABLE)
  PRINT_DEPRECATED(CXX_ENABLE EXAMPLES_ENABLE_CXX)
  FORCE_VARIABLE(EXAMPLES_ENABLE_CXX BOOL "Build ARKode C++ examples" ${CXX_ENABLE})
  UNSET(CXX_ENABLE CACHE)
ENDIF()

IF(DEFINED F90_ENABLE)
  PRINT_DEPRECATED(F90_ENABLE EXAMPLES_ENABLE_F90)
  FORCE_VARIABLE(EXAMPLES_ENABLE_F90 BOOL "Build ARKode Fortran90 examples" ${F90_ENABLE})
  UNSET(F90_ENABLE CACHE)
ENDIF()

IF(DEFINED FCMIX_ENABLE)
  PRINT_DEPRECATED(FCMIX_ENABLE F77_INTERFACE_ENABLE)
  FORCE_VARIABLE(F77_INTERFACE_ENABLE BOOL "Build Fortran 77 interfaces" ${FCMIX_ENABLE})
  UNSET(FCMIX_ENABLE CACHE)
ENDIF()

# SUNDIALS_INDEX_TYPE got new behavior
if(SUNDIALS_INDEX_TYPE)
  string(TOUPPER ${SUNDIALS_INDEX_TYPE} tmp)

  if(tmp STREQUAL "INT32_T")
    PRINT_WARNING("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for."
    "Setting SUNDIALS_INDEX_SIZE to 32 and clearing SUNDIALS_INDEX_TYPE.")
    FORCE_VARIABLE(SUNDIALS_INDEX_SIZE STRING "SUNDIALS index size" 32)
    FORCE_VARIABLE(SUNDIALS_INDEX_TYPE STRING "SUNDIALS index type" "")
  elseif(tmp STREQUAL "INT64_T")
    PRINT_WARNING("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for."
    "Setting SUNDIALS_INDEX_SIZE to 64 and clearing SUNDIALS_INDEX_TYPE.")
    FORCE_VARIABLE(SUNDIALS_INDEX_SIZE STRING "SUNDIALS index size" 64)
    FORCE_VARIABLE(SUNDIALS_INDEX_TYPE STRING "SUNDIALS index type" "")
  else()
    PRINT_WARNING("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for." "")
  endif()
endif()

if(DEFINED MPI_MPICC)
  print_deprecated(MPI_MPICC MPI_C_COMPILER)
  force_variable(MPI_C_COMPILER FILEPATH "MPI C compiler" ${MPI_MPICC})
  unset(MPI_MPICC CACHE)
endif()

if(DEFINED MPI_MPICXX)
  print_deprecated(MPI_MPICXX MPI_CXX_COMPILER)
  force_variable(MPI_CXX_COMPILER FILPATH "MPI C++ compiler" ${MPI_MPICXX})
  unset(MPI_MPICXX CACHE)
endif()

if((DEFINED MPI_MPIF77) OR (DEFINED MPI_MPIF90))
  if(DEFINED MPI_MPIF90)
    print_warning("MPI_MPIF77 and MPI_MPIF90 are deprecated and will be removed in the future." "Copying MPI_MPIF90 value to MPI_Fortran_COMPILER")
    force_variable(MPI_Fortran_COMPILER FILEPATH "MPI Fortran compiler" ${MPI_MPIF90})
  else()
    print_warning("MPI_MPIF77 and MPI_MPIF90 are deprecated and will be removed in the future." "Copying MPI_MPIF77 value to MPI_Fortran_COMPILER")
    force_variable(MPI_Fortran_COMPILER FILEPATH "MPI Fortran compiler" ${MPI_MPIF77})
  endif()
  unset(MPI_MPIF77 CACHE)
  unset(MPI_MPIF90 CACHE)
endif()

if(DEFINED MPI_RUN_COMMAND)
  print_deprecated(MPI_RUN_COMMAND MPIEXEC_EXECUTABLE)
  force_variable(MPIEXEC_EXECUTABLE FILEPATH "MPI run command" ${MPI_RUN_COMMAND})
  unset(MPI_RUN_COMMAND CACHE)
endif()
