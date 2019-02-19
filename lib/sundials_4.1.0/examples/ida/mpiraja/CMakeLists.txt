# ---------------------------------------------------------------
# Programmer:  Slaven Peles, and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for IDA raja examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
SET(IDA_examples
  "idaHeat2D_kry_p_mpiraja\;1\;4\;develop"
  )

IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
ENDIF()

# Set C++ compiler flags to include RAJA flags
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${RAJA_COMPILE_FLAGS})

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ${RAJA_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_static)
  SET(NVECS_LIB sundials_nvecmpicudaraja_static)
ELSE()
  SET(IDA_LIB sundials_ida_shared)
  SET(NVECS_LIB sundials_nvecmpicudaraja_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${IDA_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS} -L${RAJA_LIB_PATH} -lRAJA)


# Add the build and install targets for each IDA example
FOREACH(example_tuple ${IDA_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  set_source_files_properties(${example}.cu PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

  # add example source files
  CUDA_ADD_EXECUTABLE(${example} ${example}.cu)

  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} ${MPI_CXX_LIBRARIES})
  # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.cu ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/mpiraja)
  ENDIF()

ENDFOREACH(example_tuple ${IDA_examples})


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/mpiraja)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "IDA")
  SET(SOLVER_LIB "sundials_ida")
  SET(NVECTOR_LIB "sundials_nvecmpicudaraja")

  EXAMPLES2STRING(IDA_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_raja_CUDA_ex.in
    ${PROJECT_BINARY_DIR}/examples/ida/mpiraja/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/ida/mpiraja/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/mpiraja
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_raja_CUDA_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/mpiraja/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/mpiraja/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/mpiraja
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
