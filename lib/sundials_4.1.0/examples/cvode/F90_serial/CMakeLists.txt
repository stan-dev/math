# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the FCVODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(FCVODE_examples
  "cv_analytic_fp\;develop"
  "cv_analytic_sys_dns\;develop"
  "cv_analytic_sys_dns_jac\;develop"
  "cv_brusselator_dns\;develop"
  )

# Examples using LAPACK linear solvers
set(FCVODE_examples_BL
  )

# Add sparse solvers examples for 64-bit indextype configurations only,
# not compatible with 32-bit indextype
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  # Examples using KLU linear solver
  set(FCVODE_examples_KLU
    "cv_analytic_sys_klu\;develop"
    )
endif()

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
if(LINK_LIBRARY_TYPE MATCHES "static")
  set(CVODE_LIB sundials_fcvode_mod_static)
else()
  set(CVODE_LIB sundials_fcvode_mod_shared)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODE_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${FCVODE_examples})
  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_type)

  # Install fortran modules to a unique directory to avoid naming collisions
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

  # example source files
  add_executable(${example} ${example}.f90)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  SUNDIALS_ADD_TEST(${example} ${example}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS})

   # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.f90 ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F90_serial)
  endif()
endforeach(example_tuple ${FCVODE_examples})

# Add the build and install targets for each KLU example (if needed)
if(KLU_FOUND)

  # Sundials KLU linear solver module
  if(LINK_LIBRARY_TYPE MATCHES "static")
    set(SUNLINSOLKLU_LIBS
      sundials_fsunlinsolklu_mod_static)
  else()
    set(SUNLINSOLKLU_LIBS
      sundials_fsunlinsolklu_mod_shared)
  endif()

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS ${KLU_LIBRARIES})

  foreach(example_tuple ${FCVODE_examples_KLU})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # Install fortran modules to a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

    # example source files
    add_executable(${example} ${example}.f90)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    SUNDIALS_ADD_TEST(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})

     # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.f90 ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F90_serial)
    endif()

  endforeach(example_tuple ${FCVODE_examples_KLU})

endif(KLU_FOUND)

# create Makefile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F90_serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  set(SOLVER_FLIB "sundials_fcvode_mod")

  EXAMPLES2STRING(FCVODE_examples EXAMPLES)

  if(KLU_FOUND)
    EXAMPLES2STRING(FCVODE_examples_KLU EXAMPLES_KLU)
    set(SUNLS_LIB "sundials_sunlinsolklu")
    set(SUNLS_FLIB "sundials_fsunlinsolklu_mod")
  else()
    set(EXAMPLES_KLU "")
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/F90_serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/F90_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F90_serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/F90_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/F90_serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F90_serial
      RENAME Makefile
      )
  endif(UNIX)

endif(EXAMPLES_INSTALL)
