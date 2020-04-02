# ---------------------------------------------------------------
# Programmer(s): Shelby Lockhart and Cody Balos @ LLNL
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
# -----------------------------------------------------------------
# CMakeLists.txt file for CVODE OpenMPDEV examples
# -----------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cvAdvDiff_kry_ompdev\;4\;develop"
  )

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
if(LINK_LIBRARY_TYPE MATCHES "static")
  set(CVODE_LIB sundials_cvode_static)
  set(NVECOMP_LIB sundials_nvecopenmpdev_static)
else()
  set(CVODE_LIB sundials_cvode_shared)
  set(NVECOMP_LIB sundials_nvecopenmpdev_shared)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECOMP_LIB} ${EXTRA_LINK_LIBS})

# update the compilation flags to include OpenMP support
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # example source files
  add_executable(${example} ${example}.c)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # add to regression tests
  sundials_add_test(${example} ${example}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS})

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/C_openmpdev)
  endif()

endforeach(example_tuple ${CVODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/C_openmpdev)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")

  examples2string(CVODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_openmpdev_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/C_openmpdev/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/C_openmpdev/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/C_openmpdev
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_openmpdev_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/C_openmpdev/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/C_openmpdev/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/C_openmpdev
      RENAME Makefile
      )
  endif(UNIX)

  # add test_install target
  sundials_add_test_install(cvode C_openmpdev)

endif(EXAMPLES_INSTALL)
