# ---------------------------------------------------------------
# Programmer: Shelby Lockhart and Cody Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for ARKODE OpenMPDEV examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
SET(ARKODE_examples
  "ark_analytic_nonlin_ompdev\;4\;develop"
  "ark_heat1D_ompdev\;4\;develop"
  "ark_heat1D_adapt_ompdev\;4\;develop"
  )

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(ARKODE_LIB sundials_arkode_static)
  SET(NVECOMP_LIB sundials_nvecopenmpdev_static)
ELSE()
  SET(ARKODE_LIB sundials_arkode_shared)
  SET(NVECOMP_LIB sundials_nvecopenmpdev_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${ARKODE_LIB} ${NVECOMP_LIB} ${EXTRA_LINK_LIBS})

# update the compilation flags to include OpenMP support
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAFS}")
# Use C flags for linker as well.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

# Add the build and install targets for each ARKODE example
FOREACH(example_tuple ${ARKODE_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_args)
  LIST(GET example_tuple 2 example_type)

  # example source files
  ADD_EXECUTABLE(${example} ${example}.c)

  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

  # add to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

  # install example source and out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_openmpdev)
  ENDIF()

ENDFOREACH(example_tuple ${ARKODE_examples})


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_openmpdev)

  # Install the extra files
  FOREACH(extrafile ${ARKODE_extras_OMP})
    INSTALL(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_openmpdev)
  ENDFOREACH()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "ARKODE")
  SET(SOLVER_LIB "sundials_arkode")

  EXAMPLES2STRING(ARKODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_openmpdev_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/C_openmpdev/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/C_openmpdev/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_openmpdev
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_openmpdev_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/C_openmpdev/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/C_openmpdev/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_openmpdev
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
