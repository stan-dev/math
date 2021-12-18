# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# CMakeLists.txt file for the ARKode F2003 module parallel examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(FARKODE_examples
  "ark_brusselator1D_task_local_nls_f2003\;--monitor\;1\;4\;develop\;2"
  "ark_brusselator1D_task_local_nls_f2003\;--monitor --global-nls\;1\;4\;develop\;2"
  "ark_brusselator1D_task_local_nls_f2003\;--monitor --explicit --tf 3\;1\;4\;develop\;2")

if(MPI_Fortran_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link against
set(ARKODE_LIBS
  sundials_arkode
  sundials_farkode_mod
  sundials_nvecmpiplusx
  sundials_fnvecmpiplusx_mod
  sundials_nvecmpimanyvector
  sundials_fnvecmpimanyvector_mod)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${ARKODE_LIBS} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${FARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)
  list(GET example_tuple 5 example_precision)

  if(NOT TARGET ${example})
    # Install fortran modules to a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

    # example source files
    add_executable(${example} ${example}.f90)

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})

    set_target_properties(${example} PROPERTIES FOLDER "Examples")
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
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type}
    FLOAT_PRECISION ${example_precision})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.f90 ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_parallel)
  endif()

endforeach(example_tuple ${FARKODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_parallel)

  # Install the extra files
  foreach(extrafile ${ARKODE_extras})
    install(FILES ${extrafile}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_parallel)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "ARKODE")

  # Makefile: convert semi-colon separated target list to space separated string
  list2string(ARKODE_LIBS EXAMPLE_LIBS)

  # CMakeLists: replace sundials_ prefix and convert to space separted string
  list(TRANSFORM ARKODE_LIBS REPLACE "sundials_" "SUNDIALS::"
    OUTPUT_VARIABLE EXAMPLES_CMAKE_TARGETS_tmp)
  list2string(EXAMPLES_CMAKE_TARGETS_tmp EXAMPLES_CMAKE_TARGETS)

  examples2string(FARKODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/F2003_parallel/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/F2003_parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_parallel
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/F2003_parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/F2003_parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F2003_parallel
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(arkode F2003_parallel)

endif()
