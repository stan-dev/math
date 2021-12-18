# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for IDA Trilinos examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
if(Trilinos_INTERFACE_MPI_CXX_FOUND)
  set(IDA_examples
      "idaHeat2D_kry_tpetra\;1\;1\;develop"
      "idaHeat2D_kry_p_tpetra\;1\;4\;develop"
     )
else()
  set(IDA_examples
      "idaHeat2D_kry_tpetra\;1\;1\;develop"
     )
endif()

# Set Trilinos compilers/flags
set(CMAKE_CXX_COMPILER ${Trilinos_INTERFACE_CXX_COMPILER})
set(CMAKE_C_COMPILER   ${Trilinos_INTERFACE_C_COMPILER})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_INTERFACE_CXX_COMPILER_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Trilinos_INTERFACE_C_COMPILER_FLAGS}")

set(MPIEXEC_EXECUTABLE ${Trilinos_INTERFACE_MPIEXEC})

# Specify libraries to link against
set(IDA_LIB sundials_ida)
set(NVECS_LIB sundials_nvectrilinos)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${IDA_LIB} ${NVECS_LIB})

# Add the build and install targets for each IDA example
foreach(example_tuple ${IDA_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # add example source files
  add_executable(${example} ${example}.cpp)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  sundials_add_test(${example} ${example}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS} ${EXTRA_LINK_LIBS})

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cpp ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/trilinos)
  endif()

endforeach(example_tuple ${IDA_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/trilinos)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "IDA")
  set(SOLVER_LIB "sundials_ida")

  examples2string(IDA_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_trilinos_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/ida/trilinos/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/ida/trilinos/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/trilinos
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_trilinos_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/trilinos/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/trilinos/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/trilinos
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(ida trilinos)

endif()
