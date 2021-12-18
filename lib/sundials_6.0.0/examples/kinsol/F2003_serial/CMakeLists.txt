# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# ---------------------------------------------------------------
# CMakeLists.txt file for the F2003 KINSOL serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(FKINSOL_examples
  "kinRoboKin_dns_f2003\;develop"
  "kinLaplace_bnd_f2003\;develop"
  "kinLaplace_picard_kry_f2003\;develop"
  )

# Specify libraries to link against
set(KINSOL_LIB sundials_fkinsol_mod)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${KINSOL_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${FKINSOL_examples})
  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_type)

  # Install fortran modules to a unique directory to avoid naming collisions
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

  # example source files
  add_executable(${example} ${example}.f90)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  sundials_add_test(${example} ${example}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS})

   # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.f90 ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/F2003_serial)
  endif()
endforeach(example_tuple ${FKINSOL_examples})

# create Makefile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/F2003_serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "KINSOL")
  set(SOLVER_LIB "sundials_kinsol")
  set(SOLVER_FLIB "sundials_fkinsol_mod")

  examples2string(FKINSOL_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/kinsol/F2003_serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/kinsol/F2003_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/F2003_serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/kinsol/F2003_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/kinsol/F2003_serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/F2003_serial
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(kinsol F2003_serial)

endif()
