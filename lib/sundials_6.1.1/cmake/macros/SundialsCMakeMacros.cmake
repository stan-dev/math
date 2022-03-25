# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
#                Radu Serban @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
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

# Macro to force a cache variable to take on the value

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

# Macro to print warnings.

macro(print_warning message action)
  set(options )
  set(oneValueArgs MODE)
  set(multiValueArgs )

  # parse inputs and create variables print_warning_<keyword>
  cmake_parse_arguments(print_warning "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(print_warning_MODE)
    set(_mode ${print_warning_MODE})
  else()
    set(_mode WARNING)
  endif()

  set(MSG
  "------------------------------------------------------------------------\n"
  "WARNING: ${message}\n"
  "${action}\n"
  "------------------------------------------------------------------------")

  message(${_mode} ${MSG})
endmacro()

# Macro to print error messages.

macro(print_error message)
  set(options )
  set(oneValueArgs MODE)
  set(multiValueArgs )

  # parse inputs and create variables print_warning_<keyword>
  cmake_parse_arguments(print_error "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(print_error_MODE)
    set(_mode ${print_error_MODE})
  else()
    set(_mode FATAL_ERROR)
  endif()

  set(MSG
  "************************************************************************\n"
  "ERROR: ${message}\n"
  "************************************************************************")

  message(${_mode} ${MSG})
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

# Sets the SUNDIALS_GIT_VERSION variable

function(sundials_git_version)
  find_package(Git QUIET)

  set(_tmp "")

  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/.git AND ${GIT_FOUND})
    execute_process(COMMAND git describe --abbrev=12 --dirty --always --tags
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE _tmp)
    string(STRIP "${_tmp}" _tmp)
  endif()

  set(SUNDIALS_GIT_VERSION "${_tmp}" CACHE INTERNAL "")
  unset(_tmp)
endfunction()

# Macros from other files
include(SundialsAddExecutable)
include(SundialsAddLibrary)
include(SundialsAddTest)
include(SundialsAddTestInstall)
include(SundialsInstallExamples)
include(SundialsOption)
