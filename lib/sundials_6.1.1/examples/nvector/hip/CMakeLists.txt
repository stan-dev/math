# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for HIP nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

if(HIP_PLATFORM STREQUAL "nvcc")
  set(_warp_size 32)
else()
  set(_warp_size 64)
endif()

# Examples using SUNDIALS hip nvector
set(nvector_hip_examples
  "test_nvector_hip.cpp\;3 ${_warp_size} 0\;\;\;"
  "test_nvector_hip.cpp\;500 128 0\;\;\;"
  "test_nvector_hip.cpp\;1000 0 0\;\;\;"
  )

# Add source directory to include directories
include_directories(. .. ${PROJECT_SOURCE_DIR}/examples/utilities)

# Specify libraries to link against
set(NVECS_LIB sundials_nvechip)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_hip_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example_target})
    # example source files
    add_executable(${example_target} ${example})

    # link vector test utilties
    target_link_libraries(${example_target} PRIVATE test_nvector_obj)

    # folder to organize targets in an IDE
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example_target} PRIVATE ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example_target})
    else()
      string(REGEX REPLACE " " "_" test_name ${example_target}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name ${example_target}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example_target}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${nvector_hip_examples})

if(EXAMPLES_INSTALL)

  sundials_install_examples(nvector nvector_hip_examples
    SUNDIALS_TARGETS
      nvechip
    EXAMPLES_DEPENDENCIES
      test_nvector.c
    EXTRA_FILES
      ../test_nvector.c
      ../test_nvector.h
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/custom_memory_helper_gpu.h
    CMAKE_TEMPLATE
      cmakelists_HIP_ex.in
    DESTINATION
      nvector/hip)

endif(EXAMPLES_INSTALL)
