# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMake macro for installing examples.
# ---------------------------------------------------------------

# The macro:
#
#   SUNDIALS_INSTALL_EXAMPLES(<MODULE> <EXAMPLES_VAR>
#     DESTINATION path
#     CMAKE_TEMPLATE name
#     [MAKE_TEMPLATE name [SOLVER_LIBRARY target]]
#     [TEST_INSTALL target]
#     [SUNDIALS_TARGETS targets]
#     [EXTRA_FILES files]
#     [EXTRA_INCLUDES includes]
#   )
#
# adds an install target for examples in EXAMPLES_VAR that go with MODULE (e.g. arkode, nvecserial).
#
# The DESTINATION option is the path *within* EXAMPLES_INSTALL_PATH that the files should be installed.
#
# The CMAKE_TEMPLATE option is the name of the examples/templates CMake template to use (e.g. cmakelists_CXX_ex.in)
#
# The MAKE_TEMPLATE option is the name of the examples/templates Make template to use
#
# The SOLVER_LIBRARY option is used when a MAKE_TEMPLATE is provided. It should be the library name for SUNDIALS solver (e.g. arkode, cvode, ...)
#
# The TEST_INSTALL option adds a test_install target with the given target name for the MODULE.
#
# The SUNDIALS_TARGETS option is a list of CMake targets in the SUNDIALS:: namespace that the examples need to be linked to.
#
# The OTHER_TARGETS option is a list of CMake targets that the examples need to be linked to.
#
# The EXAMPLES_DEPENDENCIES option is a list of additional source files that the examples are dependent on.
#
# The EXTRA_FILES option is a list of files to install that are not example source code.
#
# The EXTRA_INCLUDES option is a list of additional includes to set with INCLUDE_DIRECTORIES.
#

macro(sundials_install_examples MODULE EXAMPLES_VAR)

  set(options )
  set(oneValueArgs SOLVER_LIBRARY DESTINATION CMAKE_TEMPLATE MAKE_TEMPLATE TEST_INSTALL)
  set(multiValueArgs SUNDIALS_TARGETS OTHER_TARGETS EXAMPLES_DEPENDENCIES EXTRA_FILES EXTRA_INCLUDES)

  # Parse keyword arguments/options
  cmake_parse_arguments(sundials_install_examples
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Install the extra files
  foreach(file ${sundials_install_examples_EXTRA_FILES})
    install(FILES ${file} DESTINATION ${EXAMPLES_INSTALL_PATH}/${sundials_install_examples_DESTINATION})
  endforeach()

  # Install the examples
  foreach(example_tuple ${${EXAMPLES_VAR}})
    list(GET example_tuple 0 example) # filename always has to be the first item in the example tuple
    get_filename_component(example_noext ${example} NAME_WE)
    file(GLOB example_out ${example_noext}*.out)
    install(FILES ${example} ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/${sundials_install_examples_DESTINATION})
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  string(TOUPPER "${MODULE}" SOLVER)
  set(SOLVER_LIB "${sundials_install_examples_SOLVER_LIBRARY}")
  set(EXAMPLES_DEPENDENCIES "${sundials_install_examples_EXAMPLES_DEPENDENCIES}")
  set(EXTRA_INCLUDES "${sundials_install_examples_EXTRA_INCLUDES}")
  examples2string(${EXAMPLES_VAR} EXAMPLES)

  # components for find_package
  list2string(sundials_install_examples_SUNDIALS_TARGETS
    EXAMPLES_CMAKE_COMPONENTS)

  set(target_list "")
  foreach(target ${sundials_install_examples_SUNDIALS_TARGETS})
    list(APPEND target_list SUNDIALS::${target})
  endforeach()
  foreach(target ${sundials_install_examples_OTHER_TARGETS})
    list(APPEND target_list ${target})
  endforeach()
  list2string(target_list EXAMPLES_CMAKE_TARGETS)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/${sundials_install_examples_CMAKE_TEMPLATE}
    ${PROJECT_BINARY_DIR}/examples/${sundials_install_examples_DESTINATION}/CMakeLists.txt
    @ONLY
    )
  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/${sundials_install_examples_DESTINATION}/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/${sundials_install_examples_DESTINATION}
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX AND (DEFINED sundials_install_examples_MAKE_TEMPLATE))
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/${sundials_install_examples_MAKE_TEMPLATE}
      ${PROJECT_BINARY_DIR}/examples/${sundials_install_examples_DESTINATION}/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/${sundials_install_examples_DESTINATION}/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/${sundials_install_examples_DESTINATION}
      RENAME Makefile
      )
  endif()

  # Add test_install target
  if(DEFINED sundials_install_examples_TEST_INSTALL)
    sundials_add_test_install(${MODULE} ${sundials_install_examples_TEST_INSTALL})
  endif()

endmacro()
