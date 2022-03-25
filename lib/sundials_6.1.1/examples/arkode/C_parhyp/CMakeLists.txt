# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for ARKODE parhyp C examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples
  "ark_diurnal_kry_ph.c\;1\;4\;develop"
  )

# Add the build and install targets for each ARKODE example
foreach(example_tuple ${ARKODE_examples})

  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # create target
  add_executable(${example_target} ${example})

  set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

  # libraries to link against
  target_link_libraries(${example_target}
    sundials_arkode
    sundials_nvecparhyp
    ${EXTRA_LINK_LIBS})

  # add to regression tests
  sundials_add_test(${example_target} ${example_target}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example_target}.out
    EXAMPLE_TYPE ${example_type})

endforeach(example_tuple ${ARKODE_examples})


if(EXAMPLES_INSTALL)

  sundials_install_examples(arkode ARKODE_examples
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parhyp_C_ex.in
    SOLVER_LIBRARY
      sundials_arkode
    SUNDIALS_TARGETS
      arkode
      nvecparhyp
    DESTINATION
      arkode/C_parhyp
    EXTRA_FILES
      README
    TEST_INSTALL
      C_parhyp
  )

endif()
