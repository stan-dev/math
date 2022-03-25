# ------------------------------------------------------------------------------
# Programmer(s): Shelby Lockhart @ UIUC/LLNL
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
# CMakeLists.txt file for Kinsol MPI examples
# ------------------------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the type is
# develop for examples excluded from 'make test' in releases.

# Examples to build
set(KINSOL_examples
  "kin_heat2D_nonlin_p.cpp\;--np 2 2\;1\;4\;develop"
  "kin_em_p.cpp\;\;1\;2\;develop")

# Header files to install
set(KINSOL_headers
  "kin_heat2D_nonlin_p.hpp"
  "kin_em_p.hpp")

# Add the build target for each example
foreach(example_tuple ${KINSOL_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  if(NOT TARGET ${example_target})

    # add example source files
    add_executable(${example_target} ${example})

    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example_target} PRIVATE
      sundials_kinsol
      sundials_nvecparallel
      sundials_nvecmpiplusx
      sundials_nvecserial
      MPI::MPI_CXX
      ${EXTRA_LINK_LIBS})

  endif()

  # set the test name
  if(example_args)
    string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
  else()
    set(test_name ${example_target})
  endif()

  # add regression test
  sundials_add_test(${test_name} ${example_target}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

endforeach()


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(kinsol KINSOL_examples
    CMAKE_TEMPLATE
      cmakelists_CXX_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_CXX_ex.in
    SOLVER_LIBRARY
      sundials_kinsol
    SUNDIALS_TARGETS
      kinsol
      nvecparallel
      nvecmpiplusx
      nvecserial
    DESTINATION
      kinsol/CXX_parallel
    EXTRA_FILES
      ${KINSOL_headers}
      README
    TEST_INSTALL
      CXX_parallel
  )

endif()
