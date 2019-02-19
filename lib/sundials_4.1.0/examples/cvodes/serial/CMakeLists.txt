# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
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
# CMakeLists.txt file for CVODES serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
SET(CVODES_examples
  "cvsRoberts_dns\;\;"
  "cvsRoberts_dns_constraints\;\;develop"
  "cvsRoberts_dns_uw\;\;develop"
  "cvsRoberts_FSA_dns\;-sensi stg1 t\;develop"
  "cvsRoberts_FSA_dns_Switch\;\;develop"
  "cvsRoberts_FSA_dns_constraints\;-sensi stg1 t\;develop"
  "cvsAdvDiff_ASAi_bnd\;\;develop"
  "cvsAdvDiff_bnd\;\;develop"
  "cvsAdvDiff_FSA_non\;-sensi stg t\;develop"
  "cvsAdvDiff_FSA_non\;-sensi sim t\;develop"
  "cvsDirectDemo_ls\;\;develop"
  "cvsDiurnal_FSA_kry\;-sensi stg t\;develop"
  "cvsDiurnal_FSA_kry\;-sensi sim t\;develop"
  "cvsDiurnal_kry_bp\;\;develop"
  "cvsDiurnal_kry\;\;develop"
  "cvsFoodWeb_ASAi_kry\;\;develop"
  "cvsFoodWeb_ASAp_kry\;\;develop"
  "cvsHessian_ASA_FSA\;\;develop"
  "cvsKrylovDemo_ls\;\;develop"
  "cvsKrylovDemo_prec\;\;develop"
  "cvsRoberts_ASAi_dns\;\;develop"
  "cvsRoberts_ASAi_dns_constraints\;\;develop"
  )

# Examples using LAPACK linear solvers
SET(CVODES_examples_BL
  "cvsAdvDiff_bndL\;\;develop"
  "cvsRoberts_dnsL\;\;develop"
  )

# Examples using KLU linear solver
SET(CVODES_examples_KLU
  "cvsRoberts_klu\;\;develop"
  "cvsRoberts_ASAi_klu\;\;develop"
  "cvsRoberts_FSA_klu\;-sensi stg1 t\;develop"
  )

# Examples using SuperLU_MT linear solver
SET(CVODES_examples_SUPERLUMT
  "cvsRoberts_sps\;\;develop"
  "cvsRoberts_ASAi_sps\;\;develop"
  "cvsRoberts_FSA_sps\;-sensi stg1 t\;develop"
  )

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE()
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF()

# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${CVODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})


# Add the build and install targets for each example
FOREACH(example_tuple ${CVODES_examples})

  # parse the example tuple
  LIST(GET example_tuple 0 example)
  LIST(GET example_tuple 1 example_args)
  LIST(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  IF(NOT TARGET ${example})
    # example source files
    ADD_EXECUTABLE(${example} ${example}.c)

    # folder to organize targets in an IDE
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  ENDIF()

  # check if example args are provided and set the test name
  IF("${example_args}" STREQUAL "")
    SET(test_name ${example})
  ELSE()
    STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  ENDIF()

  # add example to regression tests
  SUNDIALS_ADD_TEST(${test_name} ${example}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  FILE(GLOB example_out ${example}.out*)

  # install example source and .out files
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
  ENDIF()

ENDFOREACH(example_tuple ${CVODES_examples})


# Add the build and install targets for each LAPACK example (if needed)
IF(LAPACK_FOUND)

  # Sundials LAPACK linear solver modules
  IF(LINK_LIBRARY_TYPE MATCHES "static")
    SET(SUNLINSOLLAPACK_LIBS
      sundials_sunlinsollapackband_static
      sundials_sunlinsollapackdense_static)
  ELSE()
    SET(SUNLINSOLLAPACK_LIBS
      sundials_sunlinsollapackband_shared
      sundials_sunlinsollapackdense_shared)
  ENDIF()

  # LAPACK libraries
  LIST(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  # BLAS libraries
  IF(BLAS_FOUND)
    LIST(APPEND SUNLINSOLLAPACK_LIBS ${BLAS_LIBRARIES})
  ENDIF(BLAS_FOUND)

  FOREACH(example_tuple ${CVODES_examples_BL})

    # parse the example tuple
    LIST(GET example_tuple 0 example)
    LIST(GET example_tuple 1 example_args)
    LIST(GET example_tuple 2 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    IF(NOT TARGET ${example})
      # example source files
      ADD_EXECUTABLE(${example} ${example}.c)

      # folder to organize targets in an IDE
      SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} ${SUNLINSOLLAPACK_LIBS})
    ENDIF()

    # check if example args are provided and set the test name
    IF("${example_args}" STREQUAL "")
      SET(test_name ${example})
    ELSE()
      STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    ENDIF()

    # add example to regression tests
    SUNDIALS_ADD_TEST(${test_name} ${example}
      TEST_ARGS ${example_args}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type})

    # find all .out files for this example
    FILE(GLOB example_out ${example}.out*)

    # install example source and .out files
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF()

  ENDFOREACH(example_tuple ${CVODES_examples_BL})

ENDIF(LAPACK_FOUND)


# Add the build and install targets for each KLU example (if needed)
IF(KLU_FOUND)

  # Sundials KLU linear solver module
  IF(LINK_LIBRARY_TYPE MATCHES "static")
    SET(SUNLINSOLKLU_LIBS sundials_sunlinsolklu_static)
  ELSE()
    SET(SUNLINSOLKLU_LIBS sundials_sunlinsolklu_shared)
  ENDIF()

  # KLU libraries
  LIST(APPEND SUNLINSOLKLU_LIBS ${KLU_LIBRARIES})

  FOREACH(example_tuple ${CVODES_examples_KLU})

    # parse the example tuple
    LIST(GET example_tuple 0 example)
    LIST(GET example_tuple 1 example_args)
    LIST(GET example_tuple 2 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    IF(NOT TARGET ${example})
      # add example source files
      ADD_EXECUTABLE(${example} ${example}.c)

      # folder to organize targets in an IDE
      SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})
    ENDIF()

    # check if example args are provided and set the test name
    IF("${example_args}" STREQUAL "")
      SET(test_name ${example})
    ELSE()
      STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    ENDIF()

    # add example to regression tests
    SUNDIALS_ADD_TEST(${test_name} ${example}
      TEST_ARGS ${example_args}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type})

    # find all .out files for this example
    FILE(GLOB example_out ${example}.out*)

    # install example source and .out files
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF()

  ENDFOREACH(example_tuple ${CVODES_examples_KLU})

ENDIF(KLU_FOUND)


# Add the build and install targets for each SuperLU_MT example (if needed)
IF(SUPERLUMT_FOUND)

  # Sundials SuperLU_MT linear solver module
  IF(LINK_LIBRARY_TYPE MATCHES "static")
    SET(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt_static)
  ELSE()
    SET(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt_shared)
  ENDIF()

  # SuperLU_MT libraries
  LIST(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  # BLAS libraries
  IF(BLAS_FOUND)
    LIST(APPEND SUNLINSOLSLUMT_LIBS ${BLAS_LIBRARIES})
  ENDIF(BLAS_FOUND)

  FOREACH(example_tuple ${CVODES_examples_SUPERLUMT})

    # parse the example tuple
    LIST(GET example_tuple 0 example)
    LIST(GET example_tuple 1 example_args)
    LIST(GET example_tuple 2 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    IF(NOT TARGET ${example})
      # add example source files
      ADD_EXECUTABLE(${example} ${example}.c)

      # folder to organize targets in an IDE
      SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} ${SUNLINSOLSLUMT_LIBS})
    ENDIF()

    # check if example args are provided and set the test name
    IF("${example_args}" STREQUAL "")
      SET(test_name ${example})
    ELSE()
      STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    ENDIF()

    # add example to regression tests
    SUNDIALS_ADD_TEST(${test_name} ${example}
      TEST_ARGS ${example_args}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type})

    # find all .out files for this example
    FILE(GLOB example_out ${example}.out*)

    # install example source and .out files
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)
    ENDIF()

  ENDFOREACH(example_tuple ${CVODES_examples_SUPERLUMT})

ENDIF(SUPERLUMT_FOUND)


# create Makfile and CMakeLists.txt for examples
IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "CVODES")
  SET(SOLVER_LIB "sundials_cvodes")

  EXAMPLES2STRING(CVODES_examples EXAMPLES)

  IF(LAPACK_FOUND)
    EXAMPLES2STRING(CVODES_examples_BL EXAMPLES_BL)
  ELSE()
    SET(EXAMPLES_BL "")
  ENDIF()

  IF(KLU_FOUND)
    EXAMPLES2STRING(CVODES_examples_KLU EXAMPLES_KLU)
  ELSE()
    SET(EXAMPLES_KLU "")
  ENDIF()

  IF(SUPERLUMT_FOUND)
    EXAMPLES2STRING(CVODES_examples_SUPERLUMT EXAMPLES_SLUMT)
  ELSE()
    SET(EXAMPLES_SLUMTU "")
  ENDIF()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/cvodes/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/cvodes/serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvodes/serial
      RENAME Makefile
      )
  ENDIF(UNIX)

  # add test_install target
  SUNDIALS_ADD_TEST_INSTALL(cvodes cvsRoberts_dns
    EXAMPLE_DIR ${EXAMPLES_INSTALL_PATH}/cvodes/serial)

ENDIF(EXAMPLES_INSTALL)
