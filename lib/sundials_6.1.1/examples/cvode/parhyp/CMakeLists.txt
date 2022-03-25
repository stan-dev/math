# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
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
# CMakeLists.txt file for CVODE parhyp examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cvAdvDiff_non_ph.c\;2\;2\;develop"
  )

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  set(_fused_stubs sundials_cvode_fused_stubs)
endif()

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # example source files
  add_executable(${example_target} ${example})

  set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

  # libraries to link against
  target_link_libraries(${example_target}
    PRIVATE
    sundials_cvode
    sundials_nvecparhyp
    ${_fused_stubs}
    ${EXTRA_LINK_LIBS})

  # add example to regression tests
  sundials_add_test(${example_target} ${example_target}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example_target}.out
    EXAMPLE_TYPE ${example_type})

endforeach(example_tuple ${CVODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(_fused_stubs_target cvode_fused_stubs)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  sundials_install_examples(cvode CVODE_examples
    CMAKE_TEMPLATE
      cmakelists_C_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parhyp_C_ex.in
    SOLVER_LIBRARY
      sundials_cvode
    SUNDIALS_TARGETS
      cvode
      nvecparhyp
      ${_fused_stubs_target}
    DESTINATION
      cvode/parhyp
    EXTRA_FILES
      README
    TEST_INSTALL
      parhyp
  )

endif()
