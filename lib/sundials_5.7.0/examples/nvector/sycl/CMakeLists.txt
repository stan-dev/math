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
# CMakeLists.txt file for SYCL nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS sycl nvector
set(nvector_sycl_examples
  "test_nvector_sycl\;7 32 0\;\;\;"
  "test_nvector_sycl\;500 128 0\;\;\;"
  "test_nvector_sycl\;1000 0 0\;\;\;"
  )

# Dependencies for nvector examples
set(nvector_examples_dependencies
  test_nvector.c
  )

# Add source directory to include directories
include_directories(. .. ${PROJECT_SOURCE_DIR}/examples/utilities)

# Specify libraries to link against
set(NVECS_LIB sundials_nvecsycl)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_sycl_examples})

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
    add_executable(${example} ${example}.cpp ../test_nvector.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cpp
      ../test_nvector.c
      ../test_nvector.h
      ${SUNDIALS_SOURCE_DIR}/src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/sycl)
  endif()

endforeach(example_tuple ${nvector_sycl_examples})


if(EXAMPLES_INSTALL)

  # Install extra files
  install(FILES custom_memory_helper.h DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/sycl)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "nvector_sycl")
  set(EXAMPLES_CMAKE_TARGETS "nvecsycl")
  set(EXAMPLES_LIBS "-lsundials_nvecsycl")

  examples2string(nvector_sycl_examples EXAMPLES)
  examples2string(nvector_examples_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/nvector/sycl/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/nvector/sycl/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/sycl
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/nvector/sycl/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/nvector/sycl/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/sycl
      RENAME Makefile
      )
  endif()

endif()
