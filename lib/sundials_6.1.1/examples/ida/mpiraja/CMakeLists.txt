# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for IDA raja examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(IDA_examples
  "idaHeat2D_kry_p_mpiraja.cpp\;1\;4\;develop"
  )

if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_lang CUDA)
else()
  set(_lang CXX)
endif()

if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  set(OTHER_LIBS OpenMP::OpenMP_CXX)
endif()

# Add the build and install targets for each IDA example
foreach(example_tuple ${IDA_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # add example source files
  set_source_files_properties(${example} PROPERTIES LANGUAGE ${_lang})

  add_executable(${example_target} ${example})

  set_target_properties(${example_traget} PROPERTIES FOLDER "Examples")

  target_include_directories(${example_target} PUBLIC ${MPI_CXX_INCLUDE_DIRS})

  # libraries to link against
  target_link_libraries(${example_target}
    PRIVATE
    sundials_ida
    sundials_nvecmpiplusx
    sundials_nvecraja
    ${MPI_CXX_LIBRARIES}
    ${EXTRA_LINK_LIBS}
    ${OTHER_LIBS})

  # add example to regression tests
  sundials_add_test(${example_target} ${example_target}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example_target}.out
    EXAMPLE_TYPE ${example_type})

endforeach(example_tuple ${IDA_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
    set(EXAMPLES_FIND_PACKAGE "find_package(OpenMP REQUIRED)\n")
    set(_openmp_target OpenMP::OpenMP_CXX)
  endif()

  if(RAJA_NEEDS_THREADS)
    set(EXAMPLES_FIND_PACKAGE "${EXAMPLES_FIND_PACKAGE}find_package(Threads REQUIRED)\n")
  endif()

  sundials_install_examples(ida IDA_examples
    CMAKE_TEMPLATE
      cmakelists_${_lang}_MPI_ex.in
    MAKE_TEMPLATE
      makefile_parallel_RAJA_ex.in
    SOLVER_LIBRARY
      sundials_ida
    SUNDIALS_TARGETS
      ida
      nvecmpiplusx
      nvecraja
    OTHER_TARGETS
      ${_openmp_target}
    DESTINATION
      ida/mpiraja
    EXTRA_FILES
      README
    TEST_INSTALL
      mpiraja
  )

endif()
