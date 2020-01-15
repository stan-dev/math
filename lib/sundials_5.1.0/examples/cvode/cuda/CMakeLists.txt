# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# Based on CMakeLists.txt in cvode/serial by Radu Serban @ LLNL
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
# CMakeLists.txt file for CVODE CUDA examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cvAdvDiff_kry_cuda\;develop"
  "cvAdvDiff_kry_cuda_managed\;develop"
  )

# Examples using cuSolverSP linear solvers
set(CVODE_examples_cusolver
  "cvRoberts_block_cusolversp_batchqr\;develop"
)

if(SUNDIALS_INDEX_SIZE MATCHES "32")
  sundials_cuda_arch_geq(cuda_arch_ok "sm_60")
  if(cuda_arch_ok)
    set(all_examples "${CVODE_examples_cusolver}")
    # Specify libraries to link against
    if(LINK_LIBRARY_TYPE MATCHES "static")
      set(SUNDIALS_LIBS sundials_sunlinsolcusolversp_static)
    else()
      set(SUNDIALS_LIBS sundials_sunlinsolcusolversp_shared)
    endif()
  else()
    message(STATUS "cvRoberts_block_cusolversp_batchqr example excluded due to insufficient CUDA_ARCH")
  endif()
endif()

list(APPEND all_examples "${CVODE_examples}")

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
if(LINK_LIBRARY_TYPE MATCHES "static")
  set(SUNDIALS_LIBS ${SUNDIALS_LIBS}
                    sundials_cvode_static
                    sundials_nveccuda_static
                    ${EXTRA_LINK_LIBS})
else()
  set(SUNDIALS_LIBS ${SUNDIALS_LIBS}
                    sundials_cvode_shared
                    sundials_nveccuda_shared
                    ${EXTRA_LINK_LIBS})
endif()

# Add source directory to include directories
include_directories(.)

# Add the build and install targets for each CVODE example
foreach(example_tuple ${all_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_type)

  set_source_files_properties(${example}.cu PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

  # example source files
  add_executable(${example} ${example}.cu)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  sundials_add_test(${example} ${example}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cu ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/cuda)
  endif()

endforeach(example_tuple ${CVODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/cuda)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  set(NVECTOR_LIB "sundials_nveccuda")

  examples2string(CVODE_examples EXAMPLES)

  if(SUNDIALS_INDEX_SIZE MATCHES "32" AND cuda_arch_ok)
    set(SUNLS_LIB "sundials_sunlinsolcusolversp")
    examples2string(CVODE_examples_cusolver EXAMPLES_CUSOLVER)
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_CUDA_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/cuda/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/cuda/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/cuda
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_CUDA_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/cuda/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/cuda/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/cuda
      RENAME Makefile
      )
  endif(UNIX)

  # add test_install target
  sundials_add_test_install(cvode cuda)

endif(EXAMPLES_INSTALL)
