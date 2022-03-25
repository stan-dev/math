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
# CMakeLists.txt file for sunlinsol SPBCGS examples
# ---------------------------------------------------------------

# Set tolerance for linear solver test based on Sundials precision
if(SUNDIALS_PRECISION MATCHES "SINGLE")
  set(TOL "1e-5")
elseif(SUNDIALS_PRECISION MATCHES "DOUBLE")
  set(TOL "1e-13")
else()
  set(TOL "1e-16")
endif()

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS SPBCGS linear solver
set(sunlinsol_spbcgs_examples
  "test_sunlinsol_spbcgs_serial\;100 1 100 ${TOL} 0\;"
  "test_sunlinsol_spbcgs_serial\;100 2 100 ${TOL} 0\;"
  )

# Dependencies for nvector examples
set(sunlinsol_spbcgs_dependencies
  test_sunlinsol
  sundials_nvector
  sundials_linearsolver
  )

# If building F2003 tests
if (BUILD_FORTRAN_MODULE_INTERFACE)
  set(sunlinsol_spbcgs_fortran_examples
    "test_fsunlinsol_spbcgs_mod_serial\;\;")
endif()

# Add source directory to include directories
include_directories(. ../..)

# Specify libraries to link against
set(NVEC_LIB sundials_nvecserial)
set(SUNLINSOL_LIB sundials_sunlinsolspbcgs)
if (BUILD_FORTRAN_MODULE_INTERFACE)
  list(APPEND NVEC_LIB sundials_fnvecserial_mod)
  list(APPEND SUNLINSOL_LIB sundials_fsunlinsolspbcgs_mod)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVEC_LIB} ${SUNLINSOL_LIB} ${EXTRA_LINK_LIBS})


# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_spbcgs_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.c ../../test_sunlinsol.c)

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
      ../../test_sunlinsol.h
      ../../test_sunlinsol.c
      ${SUNDIALS_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
      ${SUNDIALS_SOURCE_DIR}/src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/spbcgs/serial)
  endif()

endforeach(example_tuple ${sunlinsol_spbcgs_examples})

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_spbcgs_fortran_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.f90
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/test_utilities.f90
      ${SUNDIALS_SOURCE_DIR}/examples/sunlinsol/test_sunlinsol.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # set fortran module directory to avoid name collisions
    set_target_properties(${example} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})

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
  endif()

endforeach(example_tuple ${sunlinsol_spbcgs_fortran_examples})

if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/spbcgs/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER_LIB "sundials_sunlinsolspbcgs")

  examples2string(sunlinsol_spbcgs_examples EXAMPLES)
  examples2string(sunlinsol_spbcgs_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunlinsol/spbcgs/serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/spbcgs/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/spbcgs/serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/spbcgs/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/spbcgs/serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/spbcgs/serial
      RENAME Makefile
      )
  endif()

endif()
