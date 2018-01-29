# ---------------------------------------------------------------
# Programmer:  David J. Gardner @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
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
  UNSET(EXAMPLES_ENABLE)
ENDIF()

IF(DEFINED CXX_ENABLE)
  PRINT_DEPRECATED(CXX_ENABLE EXAMPLES_ENABLE_CXX)
  FORCE_VARIABLE(EXAMPLES_ENABLE_CXX BOOL "Build ARKode C++ examples" ${CXX_ENABLE})
  UNSET(CXX_ENABLE)
ENDIF()

IF(DEFINED F90_ENABLE)
  PRINT_DEPRECATED(F90_ENABLE EXAMPLES_ENABLE_F90)
  FORCE_VARIABLE(EXAMPLES_ENABLE_F90 BOOL "Build ARKode Fortran90 examples" ${F90_ENABLE})
  UNSET(F90_ENABLE)
ENDIF()
