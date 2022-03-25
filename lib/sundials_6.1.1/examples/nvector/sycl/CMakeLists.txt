# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# CMakeLists.txt file for SYCL nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS sycl nvector
set(nvector_sycl_examples
  "test_nvector_sycl.cpp\;7 32 0\;"
  "test_nvector_sycl.cpp\;500 128 0\;"
  "test_nvector_sycl.cpp\;1000 0 0\;"
  )

# Add source directory to include directories
include_directories(. .. ${PROJECT_SOURCE_DIR}/examples/utilities)

# Add the build and install targets for each example
foreach(example_tuple ${nvector_sycl_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example_target})

    # example source files
    add_executable(${example_target} ${example})

    # folder to organize targets in an IDE
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example_target}
      PRIVATE
      test_nvector_obj
      sundials_nvecsycl
      ${EXTRA_LINK_LIBS})

  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    string(REGEX REPLACE " " "_" test_name ${example_target})
  else()
    string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example_target}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach()


if(EXAMPLES_INSTALL)

  sundials_install_examples(nvector_sycl nvector_sycl_examples
    EXAMPLES_DEPENDENCIES
      test_nvector.c
    EXTRA_FILES
      ${SUNDIALS_SOURCE_DIR}/examples/nvector/test_nvector.c
      ${SUNDIALS_SOURCE_DIR}/examples/nvector/test_nvector.h
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/custom_memory_helper_sycl.h
    CMAKE_TEMPLATE
      cmakelists_CXX_ex.in
    SUNDIALS_TARGETS
      nvecsycl
    DESTINATION
      nvector/sycl
    )

endif()
