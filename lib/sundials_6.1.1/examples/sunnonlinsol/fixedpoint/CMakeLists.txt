# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# CMakeLists.txt file for sunnonlinsol fixedpoint examples
# ------------------------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Example programs
set(examples
  "test_sunnonlinsol_fixedpoint\;\;"
  "test_sunnonlinsol_fixedpoint\;2\;"
  "test_sunnonlinsol_fixedpoint\;2 0.5\;"
)

# if building F2003 tests
if (BUILD_FORTRAN_MODULE_INTERFACE)
  set(fortran_examples
    "test_fsunnonlinsol_fixedpoint_mod\;\;"
  )
endif()

# Add source directory to include directories
include_directories(.)

# Specify libraries to link against
set(SUNDIALS_LIBS sundials_nvecserial)
list(APPEND SUNDIALS_LIBS sundials_sunnonlinsolfixedpoint)
if(BUILD_FORTRAN_MODULE_INTERFACE)
  list(APPEND SUNDIALS_LIBS sundials_fnvecserial_mod)
  list(APPEND SUNDIALS_LIBS sundials_fsunnonlinsolfixedpoint_mod)
endif()

# Set-up linker flags and link libraries
list(APPEND SUNDIALS_LIBS ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.c)

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

  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunnonlinsol/fixedpoint)
  endif()

endforeach(example_tuple ${examples})

# Add the build and install targets for each example
foreach(example_tuple ${fortran_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.f90
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/test_utilities.f90)

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

endforeach(example_tuple ${fortran_examples})


if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunnonlinsol/fixedpoint)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER_LIB "sundials_sunnonlinsolfixedpoint")
  set(LIBS "${LIBS} -lsundials_sunmatrixdense -lsundials_sunlinsoldense")

  # Set the link directory for the dense sunmatrix and linear solver library
  # The generated CMakeLists.txt does not use find_library() locate it
  set(EXTRA_LIBS_DIR "${libdir}")

  examples2string(examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunnonlinsol/fixedpoint/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/sunnonlinsol/fixedpoint/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunnonlinsol/fixedpoint
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunnonlinsol/fixedpoint/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/sunnonlinsol/fixedpoint/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunnonlinsol/fixedpoint
      RENAME Makefile
      )
  endif()

endif()
