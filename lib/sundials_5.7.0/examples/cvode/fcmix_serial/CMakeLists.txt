# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
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
# CMakeLists.txt file for the FCVODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(FCVODE_examples
  "fcvAdvDiff_bnd\;develop"
  "fcvDiurnal_kry_bp\;develop"
  "fcvDiurnal_kry\;develop"
  "fcvRoberts_dns\;develop"
  "fcvRoberts_dns_constraints\;develop"
  )

# Examples using LAPACK linear solvers
set(FCVODE_examples_BL
  "fcvRoberts_dnsL\;develop"
  )

# Add sparse solvers examples for 64-bit indextype configurations only,
# not compatible with 32-bit indextype
if(SUNDIALS_INDEX_SIZE MATCHES "64")

  # Examples using KLU linear solver
  set(FCVODE_examples_KLU
    "fcvRoberts_klu\;develop"
    )

  # Examples using SuperLU_MT linear solver
  set(FCVODE_examples_SUPERLUMT
    "fcvRoberts_sps\;develop"
    )

endif()

# Specify libraries to link against
set(FNVECS_LIB sundials_fnvecserial)

# Only static FCMIX libraries are available
set(FCVODE_LIB sundials_fcvode${_STATIC_LIB_SUFFIX})

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  list(APPEND FCVODE_LIB sundials_cvode_fused_stubs)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${FCVODE_LIB} ${FNVECS_LIB})

# Add the build and install targets for each example
foreach(example_tuple ${FCVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_type)

  # example source files
  add_executable(${example} ${example}.f)

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
    install(FILES ${example}.f ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)
  endif()

endforeach(example_tuple ${FCVODE_examples})


# Add the build and install targets for each LAPACK example (if needed)
if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)

  set(SUNLINSOLLAPACK_LIBS
    sundials_sunlinsollapackband
    sundials_sunlinsollapackdense
    sundials_fsunlinsollapackband
    sundials_fsunlinsollapackdense)

  # LAPACK libraries
  list(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  foreach(example_tuple ${FCVODE_examples_BL})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # example source files
    add_executable(${example} ${example}.f)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLLAPACK_LIBS})

     # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.f ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)
    endif()

  endforeach(example_tuple ${FCVODE_examples_BL})

endif()


# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  set(SUNLINSOLKLU_LIBS
    sundials_sunlinsolklu
    sundials_fsunlinsolklu)

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS)

  foreach(example_tuple ${FCVODE_examples_KLU})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # example source files
    add_executable(${example} ${example}.f)

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
      install(FILES ${example}.f ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)
    endif()

  endforeach(example_tuple ${FCVODE_examples_KLU})

endif()


# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  set(SUNLINSOLSLUMT_LIBS
    sundials_sunlinsolsuperlumt
    sundials_fsunlinsolsuperlumt)

  foreach(example_tuple ${FCVODE_examples_SUPERLUMT})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # example source files
    add_executable(${example} ${example}.f)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLSLUMT_LIBS})

     # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.f ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)
    endif()

  endforeach(example_tuple ${FCVODE_examples_SUPERLUMT})

endif()

# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  set(SOLVER_FLIB "sundials_fcvode")
  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  examples2string(FCVODE_examples EXAMPLES)

  if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)
    examples2string(FCVODE_examples_BL EXAMPLES_BL)
  else()
    set(EXAMPLES_BL "")
  endif()

  if(BUILD_SUNLINSOL_KLU)
    examples2string(FCVODE_examples_KLU EXAMPLES_KLU)
  else()
    set(EXAMPLES_KLU "")
  endif()

  if(BUILD_SUNLINSOL_SUPERLUMT)
    examples2string(FCVODE_examples_SUPERLUMT EXAMPLES_SLUMT)
  else()
    set(EXAMPLES_SLUMT "")
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F77_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F77_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_serial
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(cvode fcmix_serial)

endif()
