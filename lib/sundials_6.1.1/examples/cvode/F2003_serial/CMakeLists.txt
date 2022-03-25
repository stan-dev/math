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
# CMakeLists.txt file for the FCVODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(FCVODE_examples
  "cv_analytic_fp_f2003\;develop"
  "cv_analytic_sys_dns_f2003\;develop"
  "cv_analytic_sys_dns_jac_f2003\;develop"
  "cv_brusselator_dns_f2003\;develop"
  )

# Examples using LAPACK linear solvers
set(FCVODE_examples_BL
  )

# Add sparse solvers examples for 64-bit indextype configurations only,
# not compatible with 32-bit indextype
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  # Examples using KLU linear solver
  set(FCVODE_examples_KLU
    "cv_analytic_sys_klu_f2003\;develop"
    )
endif()

# Specify libraries to link against
set(CVODE_LIB sundials_cvode sundials_fcvode_mod)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_LIB sundials_cvode_fused_stubs)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODE_LIB})

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
  sundials_add_test(${example} ${example}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS} ${EXTRA_LINK_LIBS})

   # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.f90 ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F2003_serial)
  endif()
endforeach(example_tuple ${FCVODE_examples})

# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  # Sundials KLU linear solver module
  set(SUNLINSOLKLU_LIBS sundials_fsunlinsolklu_mod)

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
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})

     # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.f90 ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F2003_serial)
    endif()

  endforeach(example_tuple ${FCVODE_examples_KLU})

endif()

# create Makefile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F2003_serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  set(SOLVER_FLIB "sundials_fcvode_mod")
  set(NVEC_LIB "sundials_nvecserial")
  set(NVEC_FLIB "sundials_fnvecserial_mod")
  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  examples2string(FCVODE_examples EXAMPLES)

  if(BUILD_SUNLINSOL_KLU)
    examples2string(FCVODE_examples_KLU EXAMPLES_KLU)
    set(SUNLS_LIB "sundials_sunlinsolklu")
    set(SUNLS_FLIB "sundials_fsunlinsolklu_mod")
  else()
    set(EXAMPLES_KLU "")
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/F2003_serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/F2003_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F2003_serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/F2003_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/F2003_serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/F2003_serial
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(cvode F2003_serial)

endif()
