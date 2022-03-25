# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for sunlinsol LAPACK band examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using the LAPACK band linear solver
set(sunlinsol_lapackband_examples
  "test_sunlinsol_lapackband\;10 2 3 0\;"
  "test_sunlinsol_lapackband\;300 7 4 0\;"
  "test_sunlinsol_lapackband\;1000 8 8 0\;"
  "test_sunlinsol_lapackband\;5000 3 100 0\;"
  )

# Dependencies for nvector examples
set(sunlinsol_lapackband_dependencies
  test_sunlinsol
  sundials_nvector
  sundials_matrix
  sundials_linearsolver
  )

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(NVEC_LIB sundials_nvecserial)
set(SUNMAT_LIB sundials_sunmatrixband)
set(SUNLINSOL_LIB sundials_sunlinsollapackband)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVEC_LIB} ${SUNMAT_LIB} ${SUNLINSOL_LIB} ${EXTRA_LINK_LIBS})

list(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})

# Add the build and install targets for each sunlinsol example
foreach(example_tuple ${sunlinsol_lapackband_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.c ../test_sunlinsol.c)

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
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c
      ../test_sunlinsol.h
      ../test_sunlinsol.c
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_linearsolver.c
      ../../../src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/lapackband)
  endif()

endforeach(example_tuple ${sunlinsol_lapackband_examples})


if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/lapackband)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER_LIB "sundials_sunlinsollapackband")
  set(LIBS "${LIBS} -lsundials_sunmatrixband -lsundials_sunmatrixdense")

  # Set the link directory for the band and dense sunmatrix libraries
  # The generated CMakeLists.txt does not use find_library() locate them
  set(EXTRA_LIBS_DIR "${libdir}")

  examples2string(sunlinsol_lapackband_examples EXAMPLES_BL)
  examples2string(sunlinsol_lapackband_dependencies EXAMPLES_DEPENDENCIES_BL)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunlinsol/lapackband/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/lapackband/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/lapackband
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/lapackband/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/lapackband/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/lapackband
      RENAME Makefile
      )
  endif()

endif()
