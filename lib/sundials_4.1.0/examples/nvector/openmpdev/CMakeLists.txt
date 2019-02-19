# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Shelby Lockhart @ LLNL
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
# CMakeLists.txt file for OpenMP DEV nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS OpenMP DEV nvector
set(nvector_openmpdev_examples
  "test_nvector_openmpdev\;1000 0\;"
  "test_nvector_openmpdev\;10000 0\;"
  )

# Dependencies for nvector examples
set(nvector_examples_dependencies
  test_nvector
  sundials_nvector
  )

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
if(LINK_LIBRARY_TYPE MATCHES "static")
  set(NVECS_LIB sundials_nvecopenmpdev_static)
else()
  set(NVECS_LIB sundials_nvecopenmpdev_shared)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXTRA_LINK_LIBS})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# Use C flags for linker as well.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")


# Add the build and install targets for each example
foreach(example_tuple ${nvector_openmpdev_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.c
      ../test_nvector.c ../../../src/sundials/sundials_nvector.c)

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
  SUNDIALS_ADD_TEST(${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c
      ../test_nvector.c
      ../test_nvector.h
      ../../../src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/openmpdev)
  endif()

endforeach()


if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/openmpdev)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER_LIB "sundials_nvecopenmpdev")

  EXAMPLES2STRING(nvector_openmpdev_examples EXAMPLES)
  EXAMPLES2STRING(nvector_examples_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_openmpdev_ex.in
    ${PROJECT_BINARY_DIR}/examples/nvector/openmpdev/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/nvector/openmpdev/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/openmpdev
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_openmpdev_ex.in
      ${PROJECT_BINARY_DIR}/examples/nvector/openmpdev/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/nvector/openmpdev/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/openmpdev
      RENAME Makefile
      )
  endif()

endif()
