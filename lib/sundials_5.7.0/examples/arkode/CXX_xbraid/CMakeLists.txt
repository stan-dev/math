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
# CMakeLists.txt file for ARKODE C++ XBraid examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples
  "ark_heat2D_xbraid\;--x_print_level 0\;1\;8\;develop"
  "ark_heat2D_p_xbraid\;--np 2 2 2 --x_print_level 0\;1\;8\;develop"
  )

if(ENABLE_HYPRE AND HYPRE_FOUND)
  set(ARKODE_examples_hypre
    "ark_heat2D_hypre_pfmg_xbraid\;--np 2 2 2 --x_print_level 0\;1\;8\;develop"
    )
endif()

# Auxiliary files to install
set(ARKODE_extras
  plot_heat2D.py
  )

# Specify libraries to link against
set(ARKODE_LIB sundials_arkode)
set(ARKODE_XBRAID_LIB sundials_arkode_xbraid)
set(NVECP_LIB sundials_nvecparallel)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${ARKODE_LIB} ${ARKODE_XBRAID_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})

    # example source files
    add_executable(${example} ${example}.cpp)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} MPI::MPI_CXX)

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
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cpp ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_xbraid)
  endif()

endforeach(example_tuple ${ARKODE_examples})



if(ENABLE_HYPRE AND HYPRE_FOUND)

  # Add the build and install targets for each example
  foreach(example_tuple ${ARKODE_examples_hypre})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_args)
    list(GET example_tuple 2 number_of_nodes)
    list(GET example_tuple 3 number_of_tasks)
    list(GET example_tuple 4 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})

      # example source files
      add_executable(${example} ${example}.cpp)

      # folder to organize targets in an IDE
      set_target_properties(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      target_link_libraries(${example} ${SUNDIALS_LIBS})
      target_link_libraries(${example} SUNDIALS::HYPRE)
      target_link_libraries(${example} SUNDIALS::XBRAID)

      if(NOT MPI_CXX_COMPILER)
        target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
      endif()

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
      EXAMPLE_TYPE ${example_type})

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

    # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.cpp ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_xbraid)
    endif()

  endforeach(example_tuple ${ARKODE_examples})

endif()


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README.md DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_xbraid)

  # Install the extra files
  foreach(extrafile ${ARKODE_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_xbraid)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "ARKODE")
  set(SOLVER_LIB "sundials_arkode")

  examples2string(ARKODE_examples EXAMPLES)

  examples2string(ARKODE_examples_hypre EXAMPLES_HYPRE)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_xbraid_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/CXX_xbraid/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/CXX_xbraid/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_xbraid
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)

    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_xbraid_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/CXX_xbraid/Makefile_ex
      @ONLY
      )

    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/CXX_xbraid/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_xbraid
      RENAME Makefile
      )

  endif(UNIX)

  # add test_install target
  sundials_add_test_install(arkode CXX_xbraid)

endif(EXAMPLES_INSTALL)
