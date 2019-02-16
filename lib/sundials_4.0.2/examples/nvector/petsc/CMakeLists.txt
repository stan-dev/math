# ---------------------------------------------------------------
# Programmer:  Steven Smith @ LLNL
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
# CMakeLists.txt file for PETSc nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS PETSc nvector
SET(nvector_petsc_examples
  "test_nvector_petsc\;1000 0\;\;\;"    # run sequentially
  "test_nvector_petsc\;1000 0\;1\;4\;"  # run 4 procs
  )

# Dependencies for nvector examples
SET(nvector_examples_dependencies
  test_nvector
  sundials_nvector
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(. ..)
INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR})

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(NVECS_LIB sundials_nvecpetsc_static)
ELSE()
  SET(NVECS_LIB sundials_nvecpetsc_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${NVECS_LIB} ${EXTRA_LINK_LIBS})

LIST(APPEND SUNDIALS_LIBS ${PETSC_LIBRARIES})

IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()


# Add the build and install targets for each nvector example
FOREACH(example_tuple ${nvector_petsc_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_args)
  LIST(GET example_tuple 2 number_of_nodes)
  LIST(GET example_tuple 3 number_of_tasks)
  LIST(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  IF(NOT TARGET ${example})
    # example source files
    ADD_EXECUTABLE(${example} ${example}.c
      ../test_nvector.c ../../../src/sundials/sundials_nvector.c)

    # folder to organize targets in an IDE
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

    IF(NOT MPI_C_COMPILER)
      TARGET_LINK_LIBRARIES(${example} ${MPI_LIBRARIES})
    ENDIF()
  ENDIF()

  # check if example args are provided and set the test name
  IF("${example_args}" STREQUAL "")
    IF("${number_of_tasks}" STREQUAL "")
      SET(test_name ${example})
    ELSE()
      STRING(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    ENDIF()
  ELSE()
    IF("${number_of_tasks}" STREQUAL "")
      STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    ELSE()
      STRING(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks}_${example_args})
    ENDIF()
  ENDIF()  

  # add example to regression tests
  SUNDIALS_ADD_TEST(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c
      ../test_nvector.c
      ../test_nvector.h
      ../../../src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/petsc)
  ENDIF()

ENDFOREACH(example_tuple ${nvector_petsc_examples})


IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/petsc)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER_LIB "sundials_nvecpetsc")

  EXAMPLES2STRING(nvector_petsc_examples EXAMPLES)
  EXAMPLES2STRING(nvector_examples_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_petsc_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/nvector/petsc/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/nvector/petsc/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/petsc
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_petsc_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/nvector/petsc/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/nvector/petsc/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/petsc
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
