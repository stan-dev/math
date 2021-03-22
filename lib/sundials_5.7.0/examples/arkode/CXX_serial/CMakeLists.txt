# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                David J. Gardner @ LLNL
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
# CMakeLists.txt file for ARKODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples
  "ark_analytic_sys.cpp\;\;develop"
  "ark_heat2D.cpp\;\;develop"
  "ark_kpr_Mt.cpp\;0 5\;develop"
  "ark_kpr_Mt.cpp\;1 4\;develop"
  "ark_kpr_Mt.cpp\;2 8 0 -10\;develop"
  "ark_kpr_Mt.cpp\;0 4 1\;develop"
  "ark_kpr_Mt.cpp\;0 -4\;develop"
  "ark_kpr_Mt.cpp\;1 -5\;develop"
  "ark_kpr_Mt.cpp\;2 -5 0 -10\;develop"
  "ark_kpr_Mt.cpp\;1 -3 0 -10 0\;develop"
  "ark_kpr_Mt.cpp\;0 3 0 -10 0\;develop"
  "ark_kpr_Mt.cpp\;2 4 0 -10 0\;develop"
  )

# Auxiliary files to install
set(ARKODE_extras
  plot_heat2D.py
  plot_sol.py
  )

# Specify libraries to link against
set(ARKODE_LIB sundials_arkode)
set(NVECS_LIB sundials_nvecserial)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # add example executable
  if(NOT TARGET ${example_target})
    add_executable(${example_target} ${example})

    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    target_compile_features(${example_target} PRIVATE cxx_std_11)

    # libraries to link against
    target_link_libraries(${example_target} ${SUNDIALS_LIBS})
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

endforeach(example_tuple ${ARKODE_examples})

if(EXAMPLES_INSTALL)

  sundials_install_examples(arkode ARKODE_examples
    CMAKE_TEMPLATE
      cmakelists_CXX_ex.in
    MAKE_TEMPLATE
      makefile_serial_CXX_ex.in
    SOLVER_LIBRARY
      sundials_arkode
    SUNDIALS_TARGETS
      arkode
    DESTINATION
      arkode/CXX_serial
    EXTRA_FILES
      README
      ${ARKODE_extras}
    TEST_INSTALL
      CXX_serial
  )

endif()
