# ---------------------------------------------------------------
# Programmer:  Cody Balos @ LLNL
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
# Based on the FindPackageMultipass module by Jed Brown.
# ---------------------------------------------------------------
# Copyright Jed Brown
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ---------------------------------------------------------------
# PackageMultipass - this module defines two macros
#
# FIND_PACKAGE_MULTIPASS (Name CURRENT
#  STATES VAR0 VAR1 ...
#  DEPENDENTS DEP0 DEP1 ...)
#
#  This function creates a cache entry <UPPERCASED-Name>_CURRENT which
#  the user can set to "NO" to trigger a reconfiguration of the package.
#  The first time this function is called, the values of
#  <UPPERCASED-Name>_VAR0, ... are saved.  If <UPPERCASED-Name>_CURRENT
#  is false or if any STATE has changed since the last time
#  FIND_PACKAGE_MULTIPASS() was called, then CURRENT will be set to "NO",
#  otherwise CURRENT will be "YES".  IF not CURRENT, then
#  <UPPERCASED-Name>_DEP0, ... will be FORCED to NOTFOUND.
#  Example:
#    find_path (FOO_DIR include/foo.h)
#    FIND_PACKAGE_MULTIPASS (Foo foo_current
#      STATES DIR
#      DEPENDENTS INCLUDES LIBRARIES)
#    if (NOT foo_current)
#      # Make temporary files, run programs, etc, to determine FOO_INCLUDES and FOO_LIBRARIES
#    endif (NOT foo_current)
#
# MULTIPASS_SOURCE_RUNS (Name INCLUDES LIBRARIES SOURCE RUNS LANGUAGE)
#  Always runs the given test, use this when you need to re-run tests
#  because parent variables have made old cache entries stale. The LANGUAGE
#  variable is either C or CXX indicating which compiler the test should
#  use.
# ---------------------------------------------------------------

macro (FIND_PACKAGE_MULTIPASS _name _current)

  # convert the package name to all caps
  string (TOUPPER ${_name} _NAME)

  # copy all input args and remove _name and _current
  set (_args ${ARGV})
  list (REMOVE_AT _args 0 1)

  # initialize package status to current
  set (_states_current "YES")

  # check if the keyword STATES was input
  list (GET _args 0 _cmd)
  if (_cmd STREQUAL "STATES")

    # remove STATE from the args list and get the first state variable
    list (REMOVE_AT _args 0)
    list (GET _args 0 _state)

    # loop over the state variables until none are left or DEPENDENTS is reached
    while (_state AND NOT _state STREQUAL "DEPENDENTS")

      # get the name of the variable with the stored state variable value
      set (_stored_var PACKAGE_MULTIPASS_${_NAME}_${_state})

      # if the stored value is different from the current value, signal that
      # reconfiguration is necessary
      if (NOT "${${_stored_var}}" STREQUAL "${${_NAME}_${_state}}")
        set (_states_current "NO")
      endif ()

      # update the stored value
      set (${_stored_var} "${${_NAME}_${_state}}" CACHE INTERNAL "Stored state for ${_name}." FORCE)

      # remove the current state variable from the args list and get the next one
      list (REMOVE_AT _args 0)
      list (GET _args 0 _state)

    endwhile ()

  endif ()

  # get the name of the package status variable
  set (_stored ${_NAME}_CURRENT)

  # check if reconfiguration was requested
  if (NOT ${_stored})
    set (${_stored} "YES" CACHE BOOL "Is the configuration for ${_name} current?  Set to \"NO\" to reconfigure." FORCE)
    set (_states_current "NO")
  endif ()

  # set the output package status
  set (${_current} ${_states_current})

  # check the if dependent variables need to be cleared so that the calling
  # module can reset their values
  if (NOT ${_current} AND PACKAGE_MULTIPASS_${_name}_CALLED)

    message (STATUS "Clearing ${_name} dependent variables")

    # check if the keyword DEPENDENTS was input
    list (GET _args 0 _cmd)
    if (_cmd STREQUAL "DEPENDENTS")

      # remove the DEPENDENTS from the list
      list (REMOVE_AT _args 0)

      # clear the value of each dependent variable in the list
      foreach (dep ${_args})
        set (${_NAME}_${dep} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
      endforeach ()

    endif ()

    # reset the package FOUND status
    set (${_NAME}_FOUND "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)

  endif ()

  # signal that multipass has previously been called for this package
  set (PACKAGE_MULTIPASS_${_name}_CALLED YES CACHE INTERNAL "Private" FORCE)

endmacro (FIND_PACKAGE_MULTIPASS)


macro (MULTIPASS_SOURCE_RUNS includes libraries source runs language)
  # This is a ridiculous hack. CHECK_${language}_SOURCE_*
  # thinks that if the *name* of the return variable doesn't change,
  # then the test does not need to be re-run.  We keep an internal count
  # which we increment to guarantee that every test name is unique. If we've
  # gotten here, then the configuration has changed enough that the
  # test *needs* to be rerun.
  if (NOT MULTIPASS_TEST_COUNT)
    set (MULTIPASS_TEST_COUNT 00)
  endif (NOT MULTIPASS_TEST_COUNT)
  math (EXPR _tmp "${MULTIPASS_TEST_COUNT} + 1") # Why can't I add to a cache variable?
  set (MULTIPASS_TEST_COUNT ${_tmp} CACHE INTERNAL "Unique test ID")
  set (testname MULTIPASS_TEST_${MULTIPASS_TEST_COUNT}_${runs})
  set (testdir ${PROJECT_BINARY_DIR}/PETSC_test)
  if (NOT EXISTS ${testdir})
    file(MAKE_DIRECTORY ${testdir})
  endif ()
  set (CMAKE_REQUIRED_INCLUDES ${includes})
  set (CMAKE_REQUIRED_LIBRARIES ${libraries})
  # if MPI is available, use it for the test
  if (MPI_${language}_COMPILER)
    set (REQUIRED_COMPILER ${MPI_${language}_COMPILER})
  else ()
    set (REQUIRED_COMPILER ${CMAKE_${language}_COMPILER})
  endif ()
  if(${language} STREQUAL "C")
    set (extension c)
  else ()
    set (extension cxx)
  endif ()
  # Create simple test code
  file(WRITE ${testdir}/src.${extension} "${source}")
  # Create a CMakeLists.txt file for the test code
  file(WRITE ${testdir}/CMakeLists.txt
    "cmake_minimum_required(VERSION ${CMAKE_VERSION})\n"
    "project(ctest ${language})\n"
    "set(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "set(CMAKE_${language}_COMPILER \"${REQUIRED_COMPILER}\")\n"
    "set(CMAKE_${language}_STANDARD \"${CMAKE_${language}_STANDARD}\")\n"
    "set(CMAKE_${language}_FLAGS \"${CMAKE_${language}_FLAGS}\")\n"
    "include_directories(${CMAKE_REQUIRED_INCLUDES})\n"
    "add_executable(ctest src.${extension})\n"
    "target_link_libraries(ctest ${CMAKE_REQUIRED_LIBRARIES})\n")
  # Attempt to compile the test code
  try_compile(${testname} ${testdir} ${testdir} ctest
    OUTPUT_VARIABLE _output)
  # Write output compiling the test code
  file(WRITE ${testdir}/src.out "${_output}")
  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${testdir}/CMakeFiles)
  # Process test result
  if (${testname})
    set (${runs} TRUE)
  else ()
    set (${runs} FALSE)
  endif ()
  unset (_output)
endmacro (MULTIPASS_SOURCE_RUNS)
