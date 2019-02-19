# ---------------------------------------------------------------
# Programmer: David J. Gardner @ LLNL
#             Radu Serban @ LLNL
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
# CMake macros used throughout the SUNDIALS build system
# ---------------------------------------------------------------

# Print variable and value
IF (NOT COMMAND PRINT_VARIABLE)
  FUNCTION(PRINT_VARIABLE var)
    MESSAGE("${var} = '${${var}}'")
  ENDFUNCTION()
ENDIF()

# Macros to hide/show cached variables.
# These two macros can be used to "hide" or "show" in the
# list of cached variables various variables and/or options
# that depend on other options.
# Note that once a variable is modified, it will preserve its
# value (hidding it merely makes it internal)

# hide variable (set as internal), retain its value, and
# leave the documentation string empty
MACRO(HIDE_VARIABLE var)
  IF(DEFINED ${var})
    SET(${var} "${${var}}" CACHE INTERNAL "")
  ENDIF(DEFINED ${var})
ENDMACRO(HIDE_VARIABLE)

# show variable (set as cache)
MACRO(SHOW_VARIABLE var type doc default)
  IF(DEFINED ${var})
    # ignores <default> and forces <var> to its current value
    SET(${var} "${${var}}" CACHE "${type}" "${doc}" FORCE)
  ELSE(DEFINED ${var})
    # set to <var> to the <default> value
    SET(${var} "${default}" CACHE "${type}" "${doc}")
  ENDIF(DEFINED ${var})
ENDMACRO(SHOW_VARIABLE)

# show variable (set as cache) and overwrite (force) its value
MACRO(FORCE_VARIABLE var type doc val)
  SET(${var} "${val}" CACHE "${type}" "${doc}" FORCE)
ENDMACRO(FORCE_VARIABLE)

# Macros to append a common suffix or prefix to the elements of a list

MACRO(ADD_SUFFIX rootlist suffix)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${root}${suffix})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_SUFFIX)

MACRO(ADD_PREFIX prefix rootlist)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${prefix}${root})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_PREFIX)

# Macro to print warning that some features will be disabled
# due to some failure.

MACRO(PRINT_WARNING message action)
  SET(MSG
  "------------------------------------------------------------------------\n"
  "WARNING: ${message}\n"
  "${action}\n"
  "------------------------------------------------------------------------")
  MESSAGE(WARNING ${MSG})
ENDMACRO()

# Macro to print fatal error messages. Fatal errors will NOT allow
# a makefile to be generated

MACRO(PRINT_ERROR message)
  SET(MSG
  "************************************************************************\n"
  "ERROR: ${message}\n"
  "************************************************************************")
  MESSAGE(FATAL_ERROR ${MSG})
ENDMACRO()

# Returns an unquoted string. Note that CMake will readily turn such
# strings back into lists, due to the duality of lists and
# semicolon-separated strings. So be careful how you use it.

MACRO(LIST2STRING alist astring)
  FOREACH(elem ${${alist}})
   SET(${astring} "${${astring}} ${elem}")
  ENDFOREACH(elem)
ENDMACRO(LIST2STRING)

# Returns a string of unique example names from a list of example tuples

MACRO(EXAMPLES2STRING example_list example_string)
  SET(tmp_list "")
  FOREACH(example_tuple ${${example_list}})
    LIST(GET example_tuple 0 example)
    LIST(APPEND tmp_list ${example})
  ENDFOREACH()
  LIST(REMOVE_DUPLICATES tmp_list)
  LIST2STRING(tmp_list ${example_string})
ENDMACRO(EXAMPLES2STRING)

