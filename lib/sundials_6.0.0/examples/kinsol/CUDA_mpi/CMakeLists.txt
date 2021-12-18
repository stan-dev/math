# ------------------------------------------------------------------------------
# Programmer(s): Shelby Lockhart @ UIUC/LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# CMakeLists.txt file for KINSOL MPI+CUDA examples
# ------------------------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the type is develop
# for examples excluded from 'make test' in releases.

# Examples to build
set(KINSOL_examples
  "kin_em_mpicuda.cu\;1\;2\;develop")

# Header files to install
set(KINSOL_headers
  "kin_em_mpicuda.hpp")

# Add the build target for each example
foreach(example_tuple ${KINSOL_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # add example source files
  add_executable(${example_target} ${example})

  set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

  # include directories
  target_include_directories(${example_target} PRIVATE ${MPI_CXX_INCLUDE_DIRS})

  # libraries to link against
  target_link_libraries(${example_target} PRIVATE
    sundials_kinsol
    sundials_nvecmpiplusx
    sundials_nveccuda
    ${MPI_CXX_LIBRARIES}
    ${EXTRA_LINK_LIBS})

  # add example to regression tests
  sundials_add_test(${example_target} ${example_target}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example_target}.out
    EXAMPLE_TYPE ${example_type})

endforeach()


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(kinsol KINSOL_examples
    CMAKE_TEMPLATE
      cmakelists_CUDA_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_CUDA_ex.in
    SOLVER_LIBRARY
      sundials_kinsol
    SUNDIALS_TARGETS
      kinsol
      nvecmpiplusx
      nveccuda
    DESTINATION
      kinsol/CUDA_mpi
    EXTRA_FILES
      ${KINSOL_headers}
      README
    TEST_INSTALL
      CUDA_mpi
  )

endif()
