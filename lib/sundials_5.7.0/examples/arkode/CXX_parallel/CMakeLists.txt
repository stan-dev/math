# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU,
#                Cody J. Balos @ LLNL
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
# -----------------------------------------------------------------
# CMakeLists.txt file for ARKODE C++ parallel examples
# -----------------------------------------------------------------

# Example lists are tuples "name\;compile defs\;args\;nodes\;tasks\;type\;float precision"
# where the type is develop for examples excluded from 'make test' in releases.

# Examples to run.
set(serial_examples
  "ark_heat2D_p\;\;--np 2 2\;1\;4\;develop\;default")
set(serial_raja_examples )
set(cuda_raja_examples )
set(hip_raja_examples )
set(raja_raja_examples )
set(openmpdev_raja_examples )

# Header files to install
set(ARKODE_headers )

# Auxiliary files to install
set(ARKODE_extras
  plot_heat2D_p.py
)

# More examples that only build with double precision.
if(ENABLE_RAJA AND (SUNDIALS_PRECISION MATCHES "DOUBLE") AND
    (SUNDIALS_INDEX_SIZE MATCHES "32"))
  list(APPEND serial_raja_examples
    "ark_brusselator1D_task_local_nls\;\;--monitor\;1\;4\;develop\;2"
    "ark_brusselator1D_task_local_nls\;\;--monitor --global-nls\;1\;4\;develop\;2"
    "ark_brusselator1D_task_local_nls\;\;--monitor --explicit --tf 1\;1\;4\;develop\;2")
  if(BUILD_NVECTOR_CUDA AND (RAJA_BACKENDS MATCHES "CUDA"))
    list(APPEND cuda_raja_examples
      "ark_brusselator1D_task_local_nls\;USE_CUDA_NVEC\;--monitor\;1\;4\;develop\;2"
      "ark_brusselator1D_task_local_nls\;USE_CUDAUVM_NVEC\;--monitor\;1\;4\;exclude\;2"
      "ark_brusselator1D_task_local_nls\;USE_CUDA_NVEC\;--monitor --global-nls\;1\;4\;develop\;2")
  endif()
  if(BUILD_NVECTOR_HIP AND (RAJA_BACKENDS MATCHES "HIP"))
    list(APPEND hip_raja_examples
      "ark_brusselator1D_task_local_nls\;USE_HIP_NVEC\;--monitor\;1\;4\;exclude\;2")
  endif()
  if(BUILD_NVECTOR_RAJA AND (RAJA_BACKENDS MATCHES "CUDA"))
    list(APPEND raja_raja_examples
      "ark_brusselator1D_task_local_nls\;USE_RAJA_NVEC\;--monitor\;1\;4\;exclude\;2")
  endif()
  if(BUILD_NVECTOR_OPENMPDEV AND (RAJA_BACKENDS MATCHES "TARGET_OPENMP"))
    list(APPEND openmpdev_raja_examples
      "ark_brusselator1D_task_local_nls\;USE_OMPDEV_NVEC\;--monitor\;1\;4\;exclude\;2")
  endif()
  list(APPEND ARKODE_headers ark_brusselator1D.h)
  list(APPEND ARKODE_extras plot_brusselator1D.py)
endif()

if(NOT hip_raja_examples)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  if(CMAKE_CUDA_HOST_COMPILER)
    set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
  endif()
endif()

# macro to loop through each set of examples and build/install them
macro(build_examples examples_to_build lang)

  # Add the build and install targets for each example
  foreach(example_tuple ${${examples_to_build}})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_defines)
    list(GET example_tuple 2 example_args)
    list(GET example_tuple 3 number_of_nodes)
    list(GET example_tuple 4 number_of_tasks)
    list(GET example_tuple 5 example_type)
    list(GET example_tuple 6 example_precision)

    if("${example_defines}" STREQUAL "")
      set(example_target ${example})
    else()
      string(REGEX REPLACE " " "." example_target ${example}.${example_defines})
    endif()

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example_target})
      set_source_files_properties(${example}.cpp PROPERTIES LANGUAGE ${lang})

      # example source files
      add_executable(${example_target} ${example}.cpp)

      # folder to organize targets in an IDE
      set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

      # ensure the linker language is reset to CXX
      set_target_properties(${example_target} PROPERTIES LINKER_LANGUAGE CXX)

      # libraries to link against
      if(NOT (${CMAKE_${lang}_COMPILER} STREQUAL ${MPI_CXX_COMPILER}))
        target_include_directories(${example_target} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
        target_link_libraries(${example_target} PRIVATE ${OTHER_LIBS} ${SUNDIALS_LIBS} ${MPI_CXX_LIBRARIES})
      else()
        target_link_libraries(${example_target} PRIVATE ${OTHER_LIBS} ${SUNDIALS_LIBS})
      endif()

      # compile definitions
      if(NOT "${example_defines}" STREQUAL "")
        target_compile_definitions(${example_target} PUBLIC ${example_defines})
      endif()
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
      MPI_NPROCS ${number_of_tasks}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type}
      FLOAT_PRECISION ${example_precision})

    if(EXAMPLES_INSTALL)
      # Find all .out files for this example
      file(GLOB example_out ${example_target}*.out)

      # Install example source and out files
      install(FILES ${example}.cpp ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_parallel)
    endif()
  endforeach()
endmacro()

# Specify libraries to link against for each configuration
# and build the examples.
set(SUNDIALS_LIBS sundials_arkode sundials_nvecparallel)
set(OTHER_LIBS )
build_examples(serial_examples CXX)

set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx)
# if RAJA has HIP enabled, we have to link to HIP even if we dont use it
set(OTHER_LIBS RAJA ${EXTRA_LINK_LIBS})
if(RAJA_BACKENDS MATCHES "HIP")
  set(OTHER_LIBS hip::device ${OTHER_LIBS})
endif()
build_examples(serial_raja_examples CXX)

if(BUILD_NVECTOR_CUDA)
  set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx sundials_nveccuda)
  set(OTHER_LIBS RAJA ${EXTRA_LINK_LIBS})
  build_examples(cuda_raja_examples CUDA)
endif()

if(BUILD_NVECTOR_HIP)
  set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx sundials_nvechip)
  set(OTHER_LIBS RAJA ${EXTRA_LINK_LIBS})
  build_examples(hip_raja_examples CXX)
endif()

if(BUILD_NVECTOR_RAJA)
  set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx sundials_nvecraja)
  set(OTHER_LIBS RAJA ${EXTRA_LINK_LIBS})
  build_examples(raja_raja_examples CUDA)
endif()

if(BUILD_NVECTOR_OPENMPDEV)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx sundials_nvecopenmpdev)
  set(OTHER_LIBS OpenMP::OpenMP_CXX RAJA ${EXTRA_LINK_LIBS})
  build_examples(openmpdev_raja_examples CXX)
endif()

if(EXAMPLES_INSTALL)
  # Install the header files
  foreach(extrafile ${ARKODE_headers})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_parallel)
  endforeach()

  # Install the extra files
  foreach(extrafile ${ARKODE_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_parallel)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "ARKODE")
  set(EXAMPLES_CMAKE_TARGETS "SUNDIALS::arkode;SUNDIALS::nvecparallel;SUNDIALS::nvecmpiplusx")
  set(examples_to_install "${serial_examples}")
  if(serial_raja_examples)
    list(APPEND EXAMPLES_CMAKE_TARGETS "RAJA")
    list(APPEND examples_to_install "${serial_raja_examples}")
  endif()
  examples2string(examples_to_install EXAMPLES)

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_MPI_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/CXX_parallel/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/CXX_parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_parallel
    )

  # add test_install target
  sundials_add_test_install(arkode CXX_parallel)
endif()
