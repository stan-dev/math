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
# CMakeLists.txt file for MPIPlusX, X = CUDA NVECTOR examples.
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS MPI+cuda nvector
set(nvector_cuda_examples
  "test_nvector_mpicuda\;1000 0\;\;\;"    # run sequentially
  "test_nvector_mpicuda\;1000 0\;1\;4\;"  # run parallel on 4
  )

# Dependencies for nvector examples
set(nvector_examples_dependencies
  test_nvector
  test_mpinvector
  )

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(NVECS_LIB sundials_nvecmpiplusx sundials_nveccuda)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_cuda_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.cu)

    # link vector test utilties
    target_link_libraries(${example} PRIVATE test_nvector_obj test_nvectormpi_obj)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # Workaround bug in CMake < 3.17.3 when using MPI::MPI_CXX
    target_include_directories(${example} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(${example} PRIVATE ${MPI_CXX_LIBRARIES})

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cu
      ../test_nvector.c
      ../test_mpinvector.c
      ../test_nvector.h
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/mpicuda)
  endif()

endforeach(example_tuple ${nvector_cuda_examples})


if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/mpicuda)

  examples2string(nvector_cuda_examples EXAMPLES)
  examples2string(nvector_examples_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_CUDA_ex.in
    ${PROJECT_BINARY_DIR}/examples/nvector/mpicuda/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/nvector/mpicuda/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/mpicuda
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_CUDA_ex.in
      ${PROJECT_BINARY_DIR}/examples/nvector/mpicuda/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/nvector/mpicuda/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/mpicuda
      RENAME Makefile
      )
  endif()

endif()
