# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# Based on the FindPETSC module by Jed Brown.
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# Try to find PETSC. This has three usage modes.
#
# The first usage mode is to find PETSC by introspection.
# This case is triggered when PETSC_DIR is not set by the user.
# Setting the variables below change the behavior of the search in this mode:
#  PETSC_DIR     - directory in which PETSC resides
#  PETSC_ARCH    - build architecture
#  PETSC_CURRENT - (advanced) redo the find stage and executable tests
#  PETSC_WORKS   - (advanced) set to ON to ignore the output of the
#                  executable tests (not recommended)
#
# The second usage mode is to find PETSC based on the user-provided
# PETSC_DIR, and optionally PETSC_ARCH, variables. This case is triggered
# when just PETSC_DIR, and optionally PETSC_ARCH, are set by the user.
# Setting the variables below change the behavior of the search in this mode:
#  PETSC_DIR     - directory in which PETSC resides
#  PETSC_ARCH    - build architecture
#  PETSC_CURRENT - (advanced) redo the find stage and executable tests
#  PETSC_WORKS   - (advanced) set to ON to ignore the output of the
#                  executable tests (not recommended)
#
# The third usage mode is to 'find' PETSC based on the user-provided list
# of include directories and libraries. This mode will only use the includes
# and libraries provided in the PETSC_INCLUDES and PETSC_LIBRARIES variable.
# This case is triggered when PETSC_INCLUDES, and PETSC_LIBRARIES are set.
# Setting the variables below change the behavior of the search in this mode:
#  PETSC_LIBRARIES - (advanced) link these to use PETSC
#  PETSC_INCLUDES  - (advanced) the PETSC include directories
#  PETSC_CURRENT   - (advanced) redo the executable tests
#  PETSC_WORKS     - (advanced) set to ON to ignore the output of the
#                    executable tests (not recommended)
#
# Note that setting PETSC_LIBRARIES and PETSC_INCLUDES takes precedence over
# setting PETSC_DIR.
#
# Once done this will define the targets:
#
#  SUNDIALS::PETSC_ALL         - a CMake target for all of PETSc
#  SUNDIALS::PETSC_SYS         - a CMake target for the main PETSc library
#  SUNDIALS::PETSC_VEC         - a CMake target for the PETSc vector library
#  SUNDIALS::PETSC_MAT         - a CMake target for the PETSc matrix library
#  SUNDIALS::PETSC_DM          - a CMake target for the PETSc DM library
#  SUNDIALS::PETSC_KSP         - a CMake target for the PETSc KSP library
#  SUNDIALS::PETSC_SNES        - a CMake target for the PETSc SNES library
#  SUNDIALS::PETSC_TS          - a CMake target for the PETSc TS library
#
# It will also define the following, potentially useful, variables:
#
#  PETSC_COMPILER     - (advanced) Compiler used by PETSC, helpful to find a compatible MPI
#  PETSC_DEFINITIONS  - (advanced) Compiler switches for using PETSC
#  PETSC_MPIEXEC      - (advanced) Executable for running MPI programs
#  PETSC_INDEX_SIZE   - (internal) the size of indices in PETSC
#  PETSC_PRECISION    - (internal) the real type precision in PETSC
#  PETSC_VERSION      - (internal) Version string (MAJOR.MINOR.SUBMINOR)
#
# Usage:
#  find_package(PETSC COMPONENTS CXX)  - required if build --with-clanguage=C++ --with-c-support=0
#  find_package(PETSC COMPONENTS C)    - standard behavior of checking build using a C compiler
#  find_package(PETSC)                 - same as above
#
# Redistribution and use is allowed according to the terms of the BSD license.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# helper macros and functions
# ------------------------------------------------------------------------------

function (PETSC_GET_VERSION)
  if (EXISTS "${PETSC_INCLUDE_DIR}/petscversion.h")
    file (STRINGS "${PETSC_INCLUDE_DIR}/petscversion.h" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
    foreach (line ${vstrings})
      string (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
      list (GET fields 1 var)
      list (GET fields 2 val)
      set (${var} ${val} PARENT_SCOPE)
      set (${var} ${val})         # Also in local scope so we have access below
    endforeach ()
    if (PETSC_VERSION_RELEASE)
      if ($(PETSC_VERSION_PATCH) GREATER 0)
        set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}p${PETSC_VERSION_PATCH}" CACHE INTERNAL "PETSC version" FORCE)
      else ()
        set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}" CACHE INTERNAL "PETSC version" FORCE)
      endif ()
    else ()
      # make dev version compare higher than any patch level of a released version
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.99" CACHE INTERNAL "PETSC version" FORCE)
    endif ()
  else ()
    message (SEND_ERROR "${PETSC_INCLUDE_DIR}/petscversion.h does not exist")
  endif ()
endfunction ()

macro (PETSC_GET_VARIABLE name var)
  if (NOT DEFINED MAKE_EXECUTABLE)
    # need to find the make executable the first time this macro is used
    find_program (MAKE_EXECUTABLE NAMES make gmake)
    if (MAKE_EXECUTABLE MATCHES "NOTFOUND")
      message(SEND_ERROR "MAKE_EXECUTABLE could not be found (looked for `make` and `gmake`)")
    endif ()
  endif ()
  set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
  execute_process (COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f ${petsc_config_makefile} show VARIABLE=${name}
    OUTPUT_VARIABLE ${var}
    RESULT_VARIABLE petsc_return)
endmacro (PETSC_GET_VARIABLE)

macro (PETSC_TEST_RUNS includes libraries runs)
  if (PETSC_VERSION VERSION_GREATER 3.1)
    set (_PETSC_TSDestroy "TSDestroy(&ts)")
  else ()
    set (_PETSC_TSDestroy "TSDestroy(ts)")
  endif ()

  set (_PETSC_TEST_SOURCE "
static const char help[] = \"PETSC test program.\";
#include <petscts.h>
int main(int argc,char *argv[]) {
  PetscErrorCode ierr;
  TS ts;

  ierr = PetscInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = ${_PETSC_TSDestroy};CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
")

  multipass_source_runs ("${includes}" "${libraries}" "${_PETSC_TEST_SOURCE}" ${runs} "${PETSC_LANGUAGE_BINDINGS}")

  if (${${runs}})
    set (PETSC_EXECUTABLE_RUNS "YES" CACHE INTERNAL
      "The system can successfully run a PETSC executable" FORCE)
  else()
    set (PETSC_EXECUTABLE_RUNS "NO" CACHE INTERNAL
      "The system can NOT successfully run a PETSC executable" FORCE)
  endif ()
endmacro (PETSC_TEST_RUNS)

macro (PETSC_FIND_LIBRARY suffix name)
  # Clear any stale value, if we got here, we need to find it again
  set (PETSC_LIBRARY_${suffix} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)

  if (WIN32)
    set (libname lib${name}) # windows expects "libfoo", linux expects "foo"
  else (WIN32)
    set (libname ${name})
  endif (WIN32)

  find_library (PETSC_LIBRARY_${suffix} NAMES ${libname} HINTS ${petsc_lib_dir} NO_DEFAULT_PATH)
  set (PETSC_LIBRARIES_${suffix} "${PETSC_LIBRARY_${suffix}}" CACHE INTERNAL "PETSC ${suffix} libraries" FORCE)
  mark_as_advanced(PETSC_LIBRARY_${suffix})
endmacro (PETSC_FIND_LIBRARY suffix name)

macro (PETSC_FIND_LIBRARY_IN_LIST suffix names liblist)
  # Clear any stale value, if we got here, we need to find it again
  set (PETSC_LIBRARY_${suffix} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)

  foreach (name ${names})
    if (WIN32)
      set (libname lib${name}) # windows expects "libfoo", linux expects "foo"
    else (WIN32)
      set (libname ${name})
    endif (WIN32)
    foreach (lib ${${liblist}})
      if ("${lib}" MATCHES "${libname}[.].*")
        set (PETSC_LIBRARY_${suffix} ${lib} CACHE INTERNAL "" FORCE)
        list (REMOVE_ITEM ${liblist} ${lib})
        break ()
      endif ()
    endforeach ()
  endforeach ()
  set (PETSC_LIBRARIES_${suffix} "${PETSC_LIBRARY_${suffix}}" CACHE INTERNAL "PETSC ${suffix} libraries" FORCE)
  mark_as_advanced(PETSC_LIBRARY_${suffix})

endmacro (PETSC_FIND_LIBRARY_IN_LIST suffix names liblist)

macro (PETSC_JOIN libs deps)
  list (APPEND PETSC_LIBRARIES_${libs} ${PETSC_LIBRARIES_${deps}})
  # since list APPEND creates a new local variable in the current scope we need
  # to set the cache variable value to propagate the changes upwards
  set (PETSC_LIBRARIES_${libs} ${PETSC_LIBRARIES_${libs}} CACHE INTERNAL "PETSC ${libs} libraries" FORCE)
endmacro (PETSC_JOIN libs deps)

# ------------------------------------------------------------------------------
# FindPETSC
# ------------------------------------------------------------------------------

set (PETSC_VALID_COMPONENTS C CXX)

if (NOT PETSC_FIND_COMPONENTS)

  get_property (_enabled_langs GLOBAL PROPERTY ENABLED_LANGUAGES)
  list(FIND _enabled_langs "C" _c_index)
  if (${_c_index} GREATER -1)
    set (PETSC_LANGUAGE_BINDINGS "C")
  else ()
    set (PETSC_LANGUAGE_BINDINGS "CXX")
  endif ()

else()

  # Right now, this is designed for compatability with the --with-clanguage option, so
  # only allow one item in the components list.
  list(LENGTH ${PETSC_FIND_COMPONENTS} components_length)
  if(${components_length} GREATER 1)
    message(FATAL_ERROR "Only one component for PETSC is allowed to be specified")
  endif()
  # This is a stub for allowing multiple components should that time ever come. Perhaps
  # to also test Fortran bindings?
  foreach(component ${PETSC_FIND_COMPONENTS})
    list(FIND PETSC_VALID_COMPONENTS ${component} component_location)
    if(${component_location} EQUAL -1)
      message(FATAL_ERROR "\"${component}\" is not a valid PETSC component.")
    else()
      list(APPEND PETSC_LANGUAGE_BINDINGS ${component})
    endif()
  endforeach()

endif()

# Set which state variables to check to determine if the PETSC configuration is
# current and clear the other state variables
if (PETSC_INCLUDES OR PETSC_LIBRARIES)

  if (PETSC_INCLUDES AND PETSC_LIBRARIES)

    set (PETSC_STATES "LIBRARIES;INCLUDES" CACHE INTERNAL "" FORCE)
    set (PETSC_DIR  "" CACHE PATH "Path to the root of a PETSc installation" FORCE)
    set (PETSC_ARCH "" CACHE STRING "PETSc architecture" FORCE)

  else ()

    string (CONCAT msg
      "Both PETSC_INCLUDES and PETSC_LIBRARIES must be provided:\n"
      "  PETSC_INCLUDES=${PETSC_INCLUDES}\n"
      "  PETSC_LIBRARIES=${PETSC_LIBRARIES}")
    message (FATAL_ERROR ${msg})

  endif ()

else ()

  set (PETSC_STATES "DIR;ARCH" CACHE INTERNAL "" FORCE)
  set (PETSC_INCLUDES  "" CACHE STRING "Semi-colon separated list of PETSc include directories" FORCE)
  set (PETSC_LIBRARIES "" CACHE STRING "Semi-colon separated list of PETSc link libraries" FORCE)

endif ()

# Keep track of FindPETSC state so that we do not do the complete
# set of tests and variable lookups every time cmake is run.
include (FindPackageMultipass)
set (petsc_slaves LIBRARIES_SYS LIBRARIES_VEC LIBRARIES_MAT LIBRARIES_DM LIBRARIES_KSP LIBRARIES_SNES LIBRARIES_TS)
set (petsc_deps LIBRARY_DIR INCLUDE_DIR LIBRARIES_ INCLUDES_ COMPILER MPIEXEC EXECUTABLE_RUNS ${petsc_slaves})
find_package_multipass (PETSC petsc_config_current STATES ${PETSC_STATES} DEPENDENTS ${petsc_deps})

# This runs anytime the current configuration is not current.
# This happens either when a user sets PETSC_CURRENT=FALSE,
# or when one of the dependents given to find_package_multipass changes.
if (NOT petsc_config_current)

  if (PETSC_INCLUDES AND PETSC_LIBRARIES)

    message (STATUS "Finding PETSC using PETSC_INCLUDES and PETSC_LIBRARIES")

    # extract path from PETSC_INCLUDES
    foreach (_include_dir ${PETSC_INCLUDES})
      if (EXISTS "${_include_dir}/petsc.h")
        set (PETSC_INCLUDE_DIR "${_include_dir}" CACHE INTERNAL "Internal PETSc include directory" FORCE)
        break ()
      endif ()
    endforeach ()

    # check if the include directory was found
    if (NOT PETSC_INCLUDE_DIR)
      string (CONCAT msg
        "Could not determine PETSc include directory from PETSC_INCLUDES:\n"
        "  PETSC_INCLUDES=${PETSC_INCLUDES}\n")
      message (FATAL_ERROR ${msg})
    endif()

    # extract path from PETSC_LIBRARIES
    foreach (_library_path ${PETSC_LIBRARIES})
      get_filename_component (_library_name "${_library_path}" NAME)
      if (_library_name MATCHES "petsc")
        get_filename_component (_library_dir "${_library_path}" DIRECTORY)
        set (PETSC_LIBRARY_DIR "${_library_dir}" CACHE INTERNAL "Internal PETSc library directory" FORCE)
        break ()
      endif ()
    endforeach ()

    # check if the library directory was found
    if (NOT PETSC_LIBRARY_DIR)
      string (CONCAT msg
        "Could not DETERMINE PETSc library directory from PETSC_LIBRARIES:\n"
        "  PETSC_LIBRARIES=${PETSC_LIBRARIES}")
      message (FATAL_ERROR ${msg})
    endif()

    # set internal PETSC_DIR and PETSC_ARCH variables
    set (PETSC_DIR_  "${PETSC_LIBRARY_DIR}/.." CACHE INTERNAL "Internal PETSC_DIR"  FORCE)
    set (PETSC_ARCH_ ""                        CACHE INTERNAL "Internal PETSC_ARCH" FORCE)

  else()

    message (STATUS "Finding PETSC using PETSC_DIR")

    # find PETSC_DIR
    if (NOT PETSC_DIR)

      message (STATUS "Looking for PETSc in common install locations")

      # Debian uses versioned paths e.g /usr/lib/petscdir/3.5/
      file (GLOB DEB_PATHS "/usr/lib/petscdir/*")

      find_path (PETSC_DIR include/petsc.h
        HINTS ENV PETSC_DIR
        PATHS
        /usr/lib/petsc
        # Debian paths
        ${DEB_PATHS}
        # Arch Linux path
        /opt/petsc/linux-c-opt
        # MacPorts path
        /opt/local/lib/petsc
        $ENV{HOME}/petsc
        DOC "PETSC Directory")

      # check if PETSC_DIR was set/found
      if (NOT PETSC_DIR)

        string (CONCAT msg
          "Could not locate PETSc install directory please set:\n"
          "  - PETSC_DIR and (optionally) PETSC_ARCH\n"
          "or used the advanced options\n"
          "  - PETSC_INCLUDES and PETSC_LIBRARIES.")
        message (FATAL_ERROR ${msg})

      endif ()

    endif()

    # find PETSC_ARCH
    if (NOT PETSC_ARCH)

      set (_petsc_arches
        $ENV{PETSC_ARCH}                   # If set, use environment variable first
        linux-gnu-c-debug linux-gnu-c-opt  # Debian defaults
        x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
      set (PETSCCONF "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
      foreach (arch ${_petsc_arches})
        find_path (PETSCCONF petscconf.h
          HINTS ${PETSC_DIR}
          PATH_SUFFIXES ${arch}/include bmake/${arch}
          NO_DEFAULT_PATH)
        if (PETSCCONF)
          set (PETSC_ARCH "${arch}" CACHE STRING "PETSC build architecture" FORCE)
          break ()
        endif ()
      endforeach ()
      set (PETSCCONF "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)

    endif ()

    if (PETSC_ARCH)
      set (PETSC_INCLUDE_DIR "${PETSC_DIR}/${PETSC_ARCH}/include" CACHE INTERNAL "Internal PETSc include directory" FORCE)
      set (PETSC_LIBRARY_DIR "${PETSC_DIR}/${PETSC_ARCH}/lib"     CACHE INTERNAL "Internal PETSc library directory" FORCE)
    else ()
      set (PETSC_INCLUDE_DIR "${PETSC_DIR}/include" CACHE INTERNAL "Internal PETSc include directory" FORCE)
      set (PETSC_LIBRARY_DIR "${PETSC_DIR}/lib"     CACHE INTERNAL "Internal PETSc library directory" FORCE)
    endif ()

    # set internal PETSC_DIR and PETSC_ARCH variables
    set (PETSC_DIR_  "${PETSC_DIR}"  CACHE INTERNAL "Internal PETS_DIR"  FORCE)
    set (PETSC_ARCH_ "${PETSC_ARCH}" CACHE INTERNAL "Internal PETS_ARCH" FORCE)

  endif ()

  # Resolve the conf/rules and conf/variables files.
  # The location of these files has changed with different PETSc versions,
  # so look in a few different locations for them.
  if (EXISTS "${PETSC_LIBRARY_DIR}/petsc/conf/petscvariables") # > 3.5
    set (petsc_conf_rules "${PETSC_LIBRARY_DIR}/petsc/conf/rules")
    set (petsc_conf_variables "${PETSC_LIBRARY_DIR}/petsc/conf/variables")
  elseif (EXISTS "${PETSC_INCLUDE_DIR}/petscconf.h")   # > 2.3.3
    set (petsc_conf_rules "${PETSC_DIR_}/conf/rules")
    set (petsc_conf_variables "${PETSC_DIR_}/conf/variables")
  elseif (EXISTS "${PETSC_DIR_}/bmake/${PETSC_ARCH_}/petscconf.h") # <= 2.3.3
    set (petsc_conf_rules "${PETSC_DIR_}/bmake/common/rules")
    set (petsc_conf_variables "${PETSC_DIR_}/bmake/common/variables")
  elseif (PETSC_LIBRARIES AND PETSC_INCLUDES)
    message (FATAL_ERROR "PETSC_LIBRARIES=${PETSC_LIBRARIES} and PETSC_INCLUDES=${PETSC_INCLUDES} do not specify a valid PETSC installation")
  else ()
    message (FATAL_ERROR "PETSC_DIR=${PETSC_DIR} and PETSC_ARCH=${PETSC_ARCH} do not specify a valid PETSC installation")
  endif ()

  # ----------------------------------------------------------------------------
  # Probe the PETSc installation for information about how it was configured.
  # ----------------------------------------------------------------------------

  # Get the PETSc version
  petsc_get_version()

  # Put variables into environment since they are needed to get
  # configuration (petscvariables) in the PETSC makefile
  set (ENV{PETSC_DIR} "${PETSC_DIR_}")
  set (ENV{PETSC_ARCH} "${PETSC_ARCH_}")

  # A temporary makefile to probe the PETSC configuration
  set (petsc_config_makefile "${PROJECT_BINARY_DIR}/Makefile.petsc")
  file (WRITE "${petsc_config_makefile}"
"## This file was autogenerated by FindPETSC.cmake
# PETSC_DIR  = ${PETSC_DIR_}
# PETSC_ARCH = ${PETSC_ARCH_}
include ${petsc_conf_rules}
include ${petsc_conf_variables}
show :
\t-@echo -n \${\${VARIABLE}}
")

  # Extract information about the PETSC configuration
  petsc_get_variable (PETSC_LIB_DIR            petsc_lib_dir)
  petsc_get_variable (PETSC_EXTERNAL_LIB_BASIC petsc_libs_external)
  petsc_get_variable (PETSC_CCPPFLAGS          petsc_cpp_line)
  petsc_get_variable (PETSC_INCLUDE            petsc_include)
  petsc_get_variable (PCC                      petsc_cc)
  petsc_get_variable (PCC_FLAGS                petsc_cc_flags)
  petsc_get_variable (MPIEXEC                  petsc_mpiexec)
  petsc_get_variable (PETSC_INDEX_SIZE         petsc_index_size)
  petsc_get_variable (PETSC_PRECISION          petsc_precision)

  # We are done with the temporary Makefile, calling PETSC_GET_VARIABLE after this point is invalid!
  file (REMOVE ${petsc_config_makefile})

  # ----------------------------------------------------------------------------
  # Determine what libraries and includes are needed.
  # ----------------------------------------------------------------------------

  if (PETSC_INCLUDES AND PETSC_LIBRARIES)

    # If the user manually set PETSC_INCUDES and PETSC_LIBRARIES, we work off of
    # what they provided.

    # Make a copy of the user-provided library list to modify as libraries are
    # found and extracted
    set (PETSC_LIBRARIES_REMAINING ${PETSC_LIBRARIES})

    # Look for petscvec first, if it doesn't exist, we must be using single-library
    petsc_find_library_in_list (VEC petscvec PETSC_LIBRARIES_REMAINING)

    if (PETSC_LIBRARY_VEC)

      # libpetscsys is called libpetsc prior to 3.1 (when single-library was introduced)
      petsc_find_library_in_list (SYS  "petscsys;petsc" PETSC_LIBRARIES_REMAINING)
      petsc_find_library_in_list (MAT  petscmat         PETSC_LIBRARIES_REMAINING)
      petsc_find_library_in_list (DM   petscdm          PETSC_LIBRARIES_REMAINING)
      petsc_find_library_in_list (KSP  petscksp         PETSC_LIBRARIES_REMAINING)
      petsc_find_library_in_list (SNES petscsnes        PETSC_LIBRARIES_REMAINING)
      petsc_find_library_in_list (TS   petscts          PETSC_LIBRARIES_REMAINING)
      petsc_join (SYS  REMAINING)
      petsc_join (VEC  SYS)
      petsc_join (MAT  VEC)
      petsc_join (DM   MAT)
      petsc_join (KSP  DM)
      petsc_join (SNES KSP)
      petsc_join (TS   SNES)

      set (PETSC_LIBRARY_ALL   ${PETSC_LIBRARY_TS}   CACHE INTERNAL "All PETSC libraries" FORCE)
      set (PETSC_LIBRARIES_ALL ${PETSC_LIBRARIES_TS} CACHE INTERNAL "All PETSC libraries" FORCE)

      message (STATUS "Recognized PETSC install with separate libraries for each package")

    else ()

      # There is no libpetscvec
      set (PETSC_LIBRARY_VEC "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)

      petsc_find_library_in_list (SINGLE petsc PETSC_LIBRARIES_REMAINING)
      # Debian 9/Ubuntu 16.04 uses _real and _complex extensions when using libraries in /usr/lib/petsc.
      if (NOT PETSC_LIBRARY_SINGLE)
        petsc_find_library_in_list (SINGLE petsc_real PETSC_LIBRARIES_REMAINING)
      endif()
      if (NOT PETSC_LIBRARY_SINGLE)
        petsc_find_library_in_list (SINGLE petsc_complex PETSC_LIBRARIES_REMAINING)
      endif()

      foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
        set (PETSC_LIBRARIES_${pkg} "${PETSC_LIBRARY_SINGLE}" CACHE INTERNAL "PETSC ${pkg} libraries" FORCE)
      endforeach ()

      message (STATUS "Recognized PETSC install with single library for all packages")

    endif ()

    # At this point PETSC_LIBRARIES_REMAINING should only contain external
    # libraries needed by PETSc. These may (e.g., static build) or may not
    # (e.g., shared build) be needed to compile but are added to the package
    # libraries regardless.
    foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
      list (APPEND PETSC_LIBRARIES_${pkg} ${PETSC_LIBRARIES_REMAINING})
      # since list APPEND creates a new local variable in the current scope we need
      # to set the cache variable value to propagate the changes upwards
      set (PETSC_LIBRARIES_${pkg} ${PETSC_LIBRARIES_${pkg}} CACHE INTERNAL "PETSC ${pkg} libraries" FORCE)
    endforeach ()

    # Try to run a simple executable
    petsc_test_runs ("${PETSC_INCLUDES}" "${PETSC_LIBRARIES_TS}" petsc_works_userprovided)
    if (petsc_works_userprovided)
      message (STATUS "PETSC works with the includes and libraries given.")
    else ()
      message (STATUS "PETSC could not be used, maybe the install is broken.")
    endif ()

    # set include and library variables needed to create targets below
    set (petsc_includes_needed ${PETSC_INCLUDES})
    set (petsc_libraries_needed ${PETSC_LIBRARIES})

  else ()

    include (ResolveCompilerPaths)
    # Extract include paths and libraries from compile command line
    resolve_includes (petsc_includes_all "${petsc_cpp_line}")

    # On windows we need to make sure we're linking against the right
    # runtime library
    if (WIN32)
      if (petsc_cc_flags MATCHES "-MT")

        set (using_md False)
        foreach(flag_var
            CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
            CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
          if(${flag_var} MATCHES "/MD")
            set (using_md True)
          endif(${flag_var} MATCHES "/MD")
        endforeach(flag_var)
        if(${using_md} MATCHES "True")
          string(CONCAT msg "PETSC was built with /MT, but /MD is currently set.\n"
            "See http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F")
          message(WARNING ${msg})
        endif(${using_md} MATCHES "True")

      endif (petsc_cc_flags MATCHES "-MT")
    endif (WIN32)

    include (CorrectWindowsPaths)
    convert_cygwin_path(petsc_lib_dir)

    # Look for petscvec first, if it doesn't exist, we must be using single-library
    petsc_find_library (VEC petscvec)
    if (PETSC_LIBRARY_VEC)

      petsc_find_library (SYS  "petscsys;petsc") # libpetscsys is called libpetsc prior to 3.1 (when single-library was introduced)
      petsc_find_library (MAT  petscmat)
      petsc_find_library (DM   petscdm)
      petsc_find_library (KSP  petscksp)
      petsc_find_library (SNES petscsnes)
      petsc_find_library (TS   petscts)
      petsc_join (VEC  SYS)
      petsc_join (MAT  VEC)
      petsc_join (DM   MAT)
      petsc_join (KSP  DM)
      petsc_join (SNES KSP)
      petsc_join (TS   SNES)

      set (PETSC_LIBRARY_ALL   ${PETSC_LIBRARY_TS}   CACHE INTERNAL "All PETSC libraries" FORCE)
      set (PETSC_LIBRARIES_ALL ${PETSC_LIBRARIES_TS} CACHE INTERNAL "All PETSC libraries" FORCE)

      message (STATUS "Recognized PETSC install with separate libraries for each package")

    else ()

      # There is no libpetscvec
      set (PETSC_LIBRARY_VEC "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)

      petsc_find_library (SINGLE petsc)
      # Debian 9/Ubuntu 16.04 uses _real and _complex extensions when using libraries in /usr/lib/petsc.
      if (NOT PETSC_LIBRARY_SINGLE)
        petsc_find_library (SINGLE petsc_real)
      endif()
      if (NOT PETSC_LIBRARY_SINGLE)
        petsc_find_library (SINGLE petsc_complex)
      endif()

      foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
        set (PETSC_LIBRARIES_${pkg} "${PETSC_LIBRARY_SINGLE}" CACHE INTERNAL "PETSC ${pkg} libraries" FORCE)
      endforeach ()

      message (STATUS "Recognized PETSC install with single library for all packages")

    endif ()

    # determine the include and library variables needed to create targets below

    find_path (PETSC_INCLUDE_CONF petscconf.h HINTS "${PETSC_INCLUDE_DIR}" "${PETSC_DIR_}/bmake/${PETSC_ARCH_}" NO_DEFAULT_PATH)
    mark_as_advanced (PETSC_INCLUDE_CONF)

    set (petsc_includes_minimal ${PETSC_INCLUDE_CONF} ${PETSC_INCLUDE_DIR})

    petsc_test_runs ("${petsc_includes_minimal}" "${PETSC_LIBRARIES_TS}" petsc_works_minimal)
    if (petsc_works_minimal)

      message (STATUS "Minimal PETSC includes and libraries work.  This probably means we are building with shared libs.")
      set (petsc_includes_needed "${petsc_includes_minimal}")

    else (petsc_works_minimal) # Minimal includes fail, see if just adding full includes fixes it

      petsc_test_runs ("${petsc_includes_all}" "${PETSC_LIBRARIES_TS}" petsc_works_allincludes)
      if (petsc_works_allincludes) # It does, we just need all the includes

        string (CONCAT msg "PETSC requires extra include paths, but links correctly with only interface libraries.\n"
          "This is an unexpected configuration (but it seems to work fine).")
        message (STATUS ${msg})
        set (petsc_includes_needed ${petsc_includes_all})

      else (petsc_works_allincludes) # We are going to need to link the external libs explicitly

        resolve_libraries (petsc_libraries_external "${petsc_libs_external}")
        foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
          list (APPEND PETSC_LIBRARIES_${pkg} ${petsc_libraries_external})
          # since list APPEND creates a new local variable in the current scope we need
          # to set the cache variable value to propagate the changes upwards
          set (PETSC_LIBRARIES_${pkg} ${PETSC_LIBRARIES_${pkg}} CACHE INTERNAL "PETSC ${pkg} libraries" FORCE)
        endforeach (pkg)

        petsc_test_runs ("${petsc_includes_minimal}" "${PETSC_LIBRARIES_TS}" petsc_works_alllibraries)
        if (petsc_works_alllibraries)

          string (CONCAT msg "PETSC only need minimal includes, but requires explicit linking to all dependencies.\n"
            "This is expected when PETSC is built with static libraries.")
          message(STATUS ${msg})
          set (petsc_includes_needed ${petsc_includes_minimal})

        else (petsc_works_alllibraries)

          # It looks like we really need everything, should have listened to Matt
          set (petsc_includes_needed ${petsc_includes_all})
          petsc_test_runs ("${petsc_includes_all}" "${PETSC_LIBRARIES_TS}" petsc_works_all)
          if (petsc_works_all) # We fail anyways
            string (CONCAT msg "PETSC requires extra include paths and explicit linking to all dependencies.\n"
              "This probably means you have static libraries and something unexpected in PETSC headers.")
            message (STATUS ${msg})
          else (petsc_works_all) # We fail anyways
            message (STATUS "PETSC could not be used, maybe the install is broken.")
          endif (petsc_works_all)

        endif (petsc_works_alllibraries)

      endif (petsc_works_allincludes)

    endif (petsc_works_minimal)

    set (petsc_libraries_needed ${PETSC_LIBRARIES_ALL})

  endif ()

  # ----------------------------------------------------------------------------
  # Now we set all of the variables needed to build targets.
  # ----------------------------------------------------------------------------

  # If PETSC_WORKS is set override the executable test results. This variable
  # can be manually set to ON to force CMake to accept a given PETSC
  # configuration, but this will almost always result in a broken build.
  if (PETSC_WORKS)
    message (STATUS "Overwriting PETSc test results with PETSC_WORKS = ${PETSC_WORKS}")
    set (PETSC_EXECUTABLE_RUNS ${PETSC_WORKS} CACHE INTERNAL "Overwritten by PETSC_WORKS" FORCE)
  endif ()

  # We do an out-of-source build so __FILE__ will be an absolute path, hence __INSDIR__ is superfluous
  if (${PETSC_VERSION} VERSION_LESS 3.1)
    set (PETSC_DEFINITIONS "-D__SDIR__=\"\"" CACHE STRING "PETSC definitions" FORCE)
  else ()
    set (PETSC_DEFINITIONS "-D__INSDIR__=\"\"" CACHE STRING "PETSC definitions" FORCE)
  endif ()

  # Sometimes this can be used to assist FindMPI.cmake
  set (PETSC_COMPILER ${petsc_cc}      CACHE FILEPATH "PETSC compiler"                            FORCE)
  set (PETSC_MPIEXEC  ${petsc_mpiexec} CACHE FILEPATH "Executable for running PETSC MPI programs" FORCE)

  # Internal variables needed for configuring targets
  set (PETSC_INDEX_SIZE  ${petsc_index_size}       CACHE INTERNAL "PETSC index size"               FORCE)
  set (PETSC_PRECISION   ${petsc_precision}        CACHE INTERNAL "PETSC real type precision"      FORCE)
  set (PETSC_INCLUDES_   ${petsc_includes_needed}  CACHE INTERNAL "PETSC include paths to be used" FORCE)
  set (PETSC_LIBRARIES_  ${petsc_libraries_needed} CACHE INTERNAL "PETSC libraries to be used"     FORCE)

  # Note that we have forced values for all these choices.  If you
  # change these, you are telling the system to trust you that they
  # work.  It is likely that you will end up with a broken build.
  mark_as_advanced (PETSC_CURRENT PETSC_COMPILER PETSC_DEFINITIONS PETSC_MPIEXEC PETSC_EXECUTABLE_RUNS)

endif ()

find_package_handle_standard_args (PETSC
  REQUIRED_VARS PETSC_EXECUTABLE_RUNS
  VERSION_VAR PETSC_VERSION
  FAIL_MESSAGE "PETSC could not be found.")

# Create targets
if (PETSC_FOUND)
  if (PETSC_LIBRARY_SINGLE)
    foreach (suffix SYS VEC MAT DM KSP SNES TS ALL)
      if (NOT TARGET SUNDIALS::PETSC_${suffix})
        add_library (SUNDIALS::PETSC_${suffix} UNKNOWN IMPORTED)
        # add properties one-by-one for easier debugging
        set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDES_}")
        set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
          INTERFACE_LINK_LIBRARIES "${PETSC_LIBRARIES_}")
        set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
          INTERFACE_COMPILE_OPTIONS ${PETSC_DEFINITIONS})
        set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
          IMPORTED_LOCATION ${PETSC_LIBRARY_SINGLE})
      endif ()
    endforeach ()
  else ()
    foreach (suffix SYS VEC MAT DM KSP SNES TS ALL)
      if (PETSC_LIBRARY_${suffix} AND (NOT TARGET SUNDIALS::PETSC_${suffix}))
        add_library (SUNDIALS::PETSC_${suffix} UNKNOWN IMPORTED)
        # add properties one-by-one for easier debugging
        set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDES_}")
        set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
          INTERFACE_LINK_LIBRARIES "${PETSC_LIBRARIES_${suffix}}")
        set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
          INTERFACE_COMPILE_OPTIONS ${PETSC_DEFINITIONS})
        set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
          IMPORTED_LOCATION ${PETSC_LIBRARY_${suffix}})
      endif ()
    endforeach ()
  endif ()
endif (PETSC_FOUND)
