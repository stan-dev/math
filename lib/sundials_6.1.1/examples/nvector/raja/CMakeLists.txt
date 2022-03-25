# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody Balos @ LLNL
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
# CMakeLists.txt file for raja nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS raja nvector
set(nvector_raja_examples
  "test_nvector_raja\;1000 0\;\;\;"
  )

# Dependencies for nvector examples
set(nvector_examples_dependencies
  test_nvector
  sundials_math
)

# Add source directory to include directories
include_directories(. .. ${SUNDIALS_SOURCE_DIR}/examples/utilities)

if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_lang CUDA)
else()
  set(_lang CXX)
endif()

if(SUNDIALS_RAJA_BACKENDS MATCHES "SYCL")
  set(_custom_memhelper custom_memory_helper_sycl.h)
else()
  set(_custom_memhelper custom_memory_helper_gpu.h)
endif()

if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  set(OTHER_LIBS OpenMP::OpenMP_CXX)
endif()

# Add the build and install targets for each example
foreach(example_tuple ${nvector_raja_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})

    set_source_files_properties(${example}.cpp PROPERTIES LANGUAGE ${_lang})

    add_executable(${example} ${example}.cpp)

    # link vector test utilties
    target_link_libraries(${example} PRIVATE test_nvector_obj)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    target_link_libraries(${example}
      PRIVATE
      sundials_nvecraja
      ${EXTRA_LINK_LIBS}
      ${OTHER_LIBS})

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
    install(FILES ${example}.cpp
      ../test_nvector.c
      ../test_nvector.h
      ${SUNDIALS_SOURCE_DIR}/examples/utilities/${_custom_memhelper}
      ${SUNDIALS_SOURCE_DIR}/src/sundials/sundials_math.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/raja)
  endif()

endforeach(example_tuple ${nvector_raja_examples})


if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/raja)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
    set(NVECTOR_LIB "sundials_nveccudaraja")
  elseif(SUNDIALS_RAJA_BACKENDS MATCHES "HIP")
    set(NVECTOR_LIB "sundials_nvechipraja")
  else()
    set(NVECTOR_LIB "sundials_nvecsyclraja")
  endif()

  examples2string(nvector_raja_examples EXAMPLES)
  examples2string(nvector_examples_dependencies EXAMPLES_DEPENDENCIES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_RAJA_ex.in
    ${PROJECT_BINARY_DIR}/examples/nvector/raja/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/nvector/raja/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/raja
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_RAJA_ex.in
      ${PROJECT_BINARY_DIR}/examples/nvector/raja/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/nvector/raja/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/nvector/raja
      RENAME Makefile
      )
  endif()

endif()
