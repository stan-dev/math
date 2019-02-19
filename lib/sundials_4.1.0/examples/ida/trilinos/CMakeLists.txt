# ---------------------------------------------------------------
# Programmer:  Slaven Peles and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
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
IF(Trilinos_INTERFACE_MPI_CXX_FOUND)
  SET(IDA_examples
      "idaHeat2D_kry_tpetra\;1\;1\;develop"
      "idaHeat2D_kry_p_tpetra\;1\;4\;develop"
     )
ELSE()
  SET(IDA_examples
      "idaHeat2D_kry_tpetra\;1\;1\;develop"
     )
ENDIF()

# Set Trilinos compilers/flags
SET(CMAKE_CXX_COMPILER ${Trilinos_INTERFACE_CXX_COMPILER})
SET(CMAKE_C_COMPILER   ${Trilinos_INTERFACE_C_COMPILER})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_INTERFACE_CXX_COMPILER_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Trilinos_INTERFACE_C_COMPILER_FLAGS}")

SET(MPIEXEC_EXECUTABLE ${Trilinos_INTERFACE_MPIEXEC})

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_static)
  SET(NVECS_LIB sundials_nvectrilinos_static)
ELSE()
  SET(IDA_LIB sundials_ida_shared)
  SET(NVECS_LIB sundials_nvectrilinos_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${IDA_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each IDA example
FOREACH(example_tuple ${IDA_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # add example source files
  ADD_EXECUTABLE(${example} ${example}.cpp)

  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.cpp ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/trilinos)
  ENDIF()

ENDFOREACH(example_tuple ${IDA_examples})


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/trilinos)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "IDA")
  SET(SOLVER_LIB "sundials_ida")

  EXAMPLES2STRING(IDA_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_trilinos_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/ida/trilinos/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/ida/trilinos/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/trilinos
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_trilinos_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/trilinos/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/trilinos/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/trilinos
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
