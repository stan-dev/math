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
# CMakeLists.txt file for CVODE HIP examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(all_examples
  "cvAdvDiff_kry_hip.cpp\;\;develop"
  "cvAdvDiff_diag_hip.cpp\;0 0\;develop"
  "cvAdvDiff_diag_hip.cpp\;0 1\;develop"
  "cvAdvDiff_diag_hip.cpp\;1 1\;develop"
  )

# Specify libraries to link against
set(SUNDIALS_LIBS
  sundials_cvode
  sundials_nvechip)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  set(SUNDIALS_LIBS ${SUNDIALS_LIBS} sundials_cvode_fused_hip)
endif()

# Add source directory to include directories
include_directories(.)

# Add the build and install targets for each CVODE example
foreach(example_tuple ${all_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  if (NOT TARGET ${example_target})
    # example source files
    add_executable(${example_target} ${example})

    # folder for IDEs
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example_target} PRIVATE hip::device ${SUNDIALS_LIBS} ${EXTRA_LINK_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example_target})
  else()
    string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example_target}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})
endforeach(example_tuple ${all_examples})


# Install the targets
if(EXAMPLES_INSTALL)

  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(_fused_kernels_if_enabled cvode_fused_hip)
  else()
    set(_fused_kernels_if_enabled )
  endif()

  sundials_install_examples(cvode all_examples
    CMAKE_TEMPLATE
      cmakelists_HIP_ex.in
    SUNDIALS_TARGETS
      cvode
      nvechip
      ${_fused_kernels_if_enabled}
    DESTINATION
      cvode/hip)

endif()
