# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the F2003 CVODES serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(FCVODES_examples
  "cvs_analytic_fp_f2003\;\;develop"
  "cvsAdvDiff_FSA_non_f2003\;-sensi stg t\;develop"
  "cvsAdvDiff_FSA_non_f2003\;-sensi sim t\;develop"
  )

# Specify libraries to link against
set(CVODES_LIB sundials_fcvodes_mod)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODES_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${FCVODES_examples})
  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

   # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.f90 ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/F2003_serial)
  endif()
endforeach(example_tuple ${FCVODES_examples})

# create Makefile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/F2003_serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODES")
  set(SOLVER_LIB "sundials_cvodes")
  set(SOLVER_FLIB "sundials_fcvodes_mod")

  examples2string(FCVODES_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvodes/F2003_serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvodes/F2003_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/F2003_serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/F2003_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvodes/F2003_serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/F2003_serial
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(cvodes F2003_serial)

endif()
