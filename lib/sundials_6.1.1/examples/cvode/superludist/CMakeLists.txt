# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for CVODE SuperLU_DIST examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;type" where the
# type is develop for examples excluded from 'make test' in releases

if(SUPERLUDIST_VERSION VERSION_GREATER_EQUAL "6.3.0")
  list(APPEND CVODE_examples "cvAdvDiff_sludist\;2 1\;2\;develop")
else()
  print_warning("The cvAdvDiff_sludist example requires SuperLU_DIST >= 6.3.0, but the version found is less."
                "SUPERLUDIST_VERSION: ${SUPERLUDIST_VERSION}")
endif()

if(MPI_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  # disable C++ extensions (for known wrappers)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -DLAM_BUILDING")
else()
  include_directories(${MPI_INCLUDE_PATH})
endif()

if(SUPERLUDIST_OpenMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Specify libraries to link against
set(CVODE_LIB sundials_cvode)
set(NVECP_LIB sundials_nvecparallel)
set(SUNMAT_LIB sundials_sunmatrixslunrloc)
set(SUNLS_LIB sundials_sunlinsolsuperludist)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_LIB sundials_cvode_fused_stubs)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECP_LIB} ${SUNMAT_LIB} ${SUNLS_LIB})

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 example_type)

  # example source files
  add_executable(${example} ${example}.cpp)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  sundials_add_test(${example} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_nodes}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}_${SUNDIALS_INDEX_SIZE}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS} ${EXTRA_LINK_LIBS})

  if(NOT MPI_C_COMPILER)
    target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  endif()

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cpp ${example}_${SUNDIALS_INDEX_SIZE}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/superludist)
  endif()

endforeach(example_tuple ${CVODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/superludist)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  string(CONCAT SUNDIALS_LIBS "-lsundials_cvode "
    "-lsundials_nvecparallel "
    "-lsundials_sunmatrixslunrloc "
    "-lsundials_sunlinsolsuperludist")
  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  examples2string(CVODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_superludist_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/superludist/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/superludist/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/superludist
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_superludist_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/superludist/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/superludist/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/superludist
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(cvode superludist)

endif()
