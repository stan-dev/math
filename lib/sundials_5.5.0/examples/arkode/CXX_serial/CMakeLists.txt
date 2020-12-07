# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                David J. Gardner @ LLNL
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
# CMakeLists.txt file for ARKODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples
  "ark_analytic_sys\;\;develop"
  "ark_heat2D\;\;develop"
  "ark_kpr_Mt\;0 5\;develop"
  "ark_kpr_Mt\;1 4\;develop"
  "ark_kpr_Mt\;2 8 0 -10\;develop"
  "ark_kpr_Mt\;0 4 1\;develop"
  "ark_kpr_Mt\;0 -4\;develop"
  "ark_kpr_Mt\;1 -5\;develop"
  "ark_kpr_Mt\;2 -5 0 -10\;develop"
  "ark_kpr_Mt\;1 -3 0 -10 0\;develop"
  "ark_kpr_Mt\;0 3 0 -10 0\;develop"
  "ark_kpr_Mt\;2 4 0 -10 0\;develop"
  )

# Auxiliary files to install
set(ARKODE_extras
  plot_heat2D.py
  plot_sol.py
  )

# Specify libraries to link against
set(ARKODE_LIB sundials_arkode)
set(NVECS_LIB sundials_nvecserial)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})


# Add the build and install targets for each example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # add example executable
  if(NOT TARGET ${example})
    add_executable(${example} ${example}.cpp)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    target_compile_features(${example} PRIVATE cxx_std_11)

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
    install(FILES ${example}.cpp ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_serial)
  endif()

endforeach(example_tuple ${ARKODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_serial)

  # Install the extra files
  foreach(extrafile ${ARKODE_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_serial)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "ARKODE")
  set(SOLVER_LIB "sundials_arkode")
  set(EXAMPLES_CMAKE_TARGETS "arkode")

  examples2string(ARKODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/CXX_serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/CXX_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/CXX_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/CXX_serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_serial
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(arkode CXX_serial)

endif()
