# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# ARKode unit_tests level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# List of test tuples of the form "name\;args"
set(ARKODE_unit_tests
  "ark_test_arkstepsetforcing\;1 0"
  "ark_test_arkstepsetforcing\;1 1"
  "ark_test_arkstepsetforcing\;1 2"
  "ark_test_arkstepsetforcing\;1 3"
  "ark_test_arkstepsetforcing\;1 4"
  "ark_test_arkstepsetforcing\;1 5"
  "ark_test_arkstepsetforcing\;1 3 2.0 10.0"
  "ark_test_arkstepsetforcing\;1 3 2.0 10.0 2.0 8.0"
  "ark_test_arkstepsetforcing\;1 3 2.0 10.0 1.0 5.0"
  "ark_interp_check\;-100"
  "ark_interp_check\;-10000"
  "ark_interp_check\;-1000000"
  )

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
if(LINK_LIBRARY_TYPE MATCHES "static")
  set(ARKODE_LIB sundials_arkode_static)
  set(NVECS_LIB sundials_nvecserial_static)
else()
  set(ARKODE_LIB sundials_arkode_shared)
  set(NVECS_LIB sundials_nvecserial_shared)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each ARKODE test
foreach(test_tuple ${ARKODE_unit_tests})

  # parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 test_args)

  # check if this test has already been added, only need to add
  # test source files once for testing with different inputs
  if(NOT TARGET ${test})

    # test source files
    add_executable(${test} ${test}.c)

    set_target_properties(${test} PROPERTIES FOLDER "unit_tests")

    # include location of public and private header files
    target_include_directories(${test} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_include_directories(${test} PRIVATE ${CMAKE_SOURCE_DIR}/src)

    # libraries to link against
    target_link_libraries(${test} ${SUNDIALS_LIBS})

  endif()

  # check if test args are provided and set the test name
  if("${test_args}" STREQUAL "")
    set(test_name ${test})
  else()
    string(REPLACE " " "_" test_name "${test}_${test_args}")
    string(REPLACE " " ";" test_args "${test_args}")
  endif()

  # add test to regression tests
  add_test(NAME ${test_name} COMMAND ${test} ${test_args})

endforeach()
