# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
#                Radu Serban @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
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
if (NOT COMMAND PRINT_VARIABLE)
  function(PRINT_VARIABLE var)
    message("${var} = '${${var}}'")
  endfunction()
endif()

# Macros to hide/show cached variables.
# These two macros can be used to "hide" or "show" in the
# list of cached variables various variables and/or options
# that depend on other options.
# Note that once a variable is modified, it will preserve its
# value (hidding it merely makes it internal)

# hide variable (set as internal), retain its value, and
# leave the documentation string empty
macro(HIDE_VARIABLE var)
  if(DEFINED ${var})
    set(${var} "${${var}}" CACHE INTERNAL "")
  endif(DEFINED ${var})
endmacro(HIDE_VARIABLE)

# show variable (set as cache)
macro(SHOW_VARIABLE var type doc default)
  if(DEFINED ${var})
    # ignores <default> and forces <var> to its current value
    set(${var} "${${var}}" CACHE "${type}" "${doc}" FORCE)
  else(DEFINED ${var})
    # set to <var> to the <default> value
    set(${var} "${default}" CACHE "${type}" "${doc}")
  endif(DEFINED ${var})
endmacro(SHOW_VARIABLE)

# show variable (set as cache) and overwrite (force) its value
macro(FORCE_VARIABLE var type doc val)
  set(${var} "${val}" CACHE "${type}" "${doc}" FORCE)
endmacro(FORCE_VARIABLE)

# Macros to append a common suffix or prefix to the elements of a list

macro(ADD_SUFFIX rootlist suffix)
  set(outlist )
  foreach(root ${${rootlist}})
    list(APPEND outlist ${root}${suffix})
  endforeach(root)
  set(${rootlist} ${outlist})
endmacro(ADD_SUFFIX)

macro(ADD_PREFIX prefix rootlist)
  set(outlist )
  foreach(root ${${rootlist}})
    list(APPEND outlist ${prefix}${root})
  endforeach(root)
  set(${rootlist} ${outlist})
endmacro(ADD_PREFIX)

# Macro to print warning that some features will be disabled
# due to some failure.

macro(PRINT_WARNING message action)
  set(MSG
  "------------------------------------------------------------------------\n"
  "WARNING: ${message}\n"
  "${action}\n"
  "------------------------------------------------------------------------")
  message(WARNING ${MSG})
endmacro()

# Macro to print fatal error messages. Fatal errors will NOT allow
# a makefile to be generated

macro(PRINT_ERROR message)
  set(MSG
  "************************************************************************\n"
  "ERROR: ${message}\n"
  "************************************************************************")
  message(FATAL_ERROR ${MSG})
endmacro()

# Returns an unquoted string. Note that CMake will readily turn such
# strings back into lists, due to the duality of lists and
# semicolon-separated strings. So be careful how you use it.

macro(LIST2STRING alist astring)
  foreach(elem ${${alist}})
   set(${astring} "${${astring}} ${elem}")
  endforeach(elem)
endmacro(LIST2STRING)

# Returns a string of unique example names from a list of example tuples

macro(EXAMPLES2STRING example_list example_string)
  set(tmp_list "")
  foreach(example_tuple ${${example_list}})
    list(GET example_tuple 0 example)
    list(APPEND tmp_list ${example})
  endforeach()
  list(REMOVE_DUPLICATES tmp_list)
  list2string(tmp_list ${example_string})
endmacro(EXAMPLES2STRING)

# Macros from other files
include(SundialsAddF2003InterfaceLibrary)
include(SundialsAddTest)
include(SundialsAddTestInstall)
include(SundialsCudaArchGeq)
include(SundialsOption)
