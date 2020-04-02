# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
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
# CMakeLists.txt file for IDAS serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(IDAS_examples
  "idasRoberts_dns\;\;"
  "idasRoberts_FSA_dns\;-sensi stg t\;develop"
  "idasRoberts_ASAi_dns\;\;develop"
  "idasAkzoNob_dns\;\;develop"
  "idasAkzoNob_ASAi_dns\;\;develop"
  "idasFoodWeb_bnd\;\;develop"
  "idasHeat2D_bnd\;\;develop"
  "idasHeat2D_kry\;\;develop"
  "idasHessian_ASA_FSA\;\;develop"
  "idasKrylovDemo_ls\;\;develop"
  "idasSlCrank_dns\;\;develop"
  "idasSlCrank_FSA_dns\;\;develop"
  )

# Examples using LAPACK linear solvers
set(IDAS_examples_BL
  )

# Example using KLU linear solver
set(IDAS_examples_KLU
  "idasRoberts_klu\;\;develop"
  "idasRoberts_ASAi_klu\;\;develop"
  "idasRoberts_FSA_klu\;-sensi stg t\;develop"
  )

# Examples using SuperLU_MT linear solver
set(IDAS_examples_SUPERLUMT
  "idasRoberts_sps\;\;develop"
  "idasRoberts_ASAi_sps\;\;develop"
  "idasRoberts_FSA_sps\;-sensi stg t\;develop"
  )

# Specify libraries to link against (through the target that was used to
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
if(LINK_LIBRARY_TYPE MATCHES "static")
  set(IDAS_LIB sundials_idas_static)
  set(NVECS_LIB sundials_nvecserial_static)
else()
  set(IDAS_LIB sundials_idas_shared)
  set(NVECS_LIB sundials_nvecserial_shared)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${IDAS_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})


# Add the build and install targets for each example
foreach(example_tuple ${IDAS_examples})

  # first item is example
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

  # install example source and .out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
  endif()

endforeach(example_tuple ${IDAS_examples})


# Add the build and install targets for each LAPACK example (if needed)
if(LAPACK_FOUND)

  # Sundials LAPACK linear solver modules
  if(LINK_LIBRARY_TYPE MATCHES "static")
    set(SUNLINSOLLAPACK_LIBS
      sundials_sunlinsollapackband_static
      sundials_sunlinsollapackdense_static)
  else()
    set(SUNLINSOLLAPACK_LIBS
      sundials_sunlinsollapackband_shared
      sundials_sunlinsollapackdense_shared)
  endif()

  # LAPACK libraries
  list(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  # BLAS libraries
  if(BLAS_FOUND)
    list(APPEND SUNLINSOLLAPACK_LIBS ${BLAS_LIBRARIES})
  endif(BLAS_FOUND)

  foreach(example_tuple ${IDAS_examples_BL})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_args)
    list(GET example_tuple 2 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})
      # example source files
      add_executable(${example} ${example}.c)

      # folder to organize targets in an IDE
      set_target_properties(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLLAPACK_LIBS})
    endif()

    # check if example args are provided and set the test name
    # IF("${example_args}" STREQUAL "")
    #   SET(test_name ${example})
    # ELSE()
    #   STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    # ENDIF()

    # add example to regression tests
    sundials_add_test(${example} ${example}
      TEST_ARGS ${example_args}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

    # install example source and .out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
    endif()

  endforeach(example_tuple ${IDAS_examples_BL})

endif(LAPACK_FOUND)


# Add the build and install targets for each KLU example (if needed)
if(KLU_FOUND)

  # Sundials KLU linear solver module
  if(LINK_LIBRARY_TYPE MATCHES "static")
    set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu_static)
  else()
    set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu_shared)
  endif()

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS ${KLU_LIBRARIES})

  foreach(example_tuple ${IDAS_examples_KLU})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_args)
    list(GET example_tuple 2 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})
      # example source files
      add_executable(${example} ${example}.c)

      # folder to organize targets in an IDE
      set_target_properties(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})
    endif()

    # check if example args are provided and set the test name
    if("${example_args}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    endif()

    # add example to regression tests
    sundials_add_test(${test_name} ${example}
      TEST_ARGS ${example_args}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type})

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

    # install example source and .out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
    endif()

  endforeach(example_tuple ${IDAS_examples_KLU})

endif(KLU_FOUND)


# Add the build and install targets for each SuperLU_MT example (if needed)
if(SUPERLUMT_FOUND)

  # Sundials SuperLU_MT linear solver module
  if(LINK_LIBRARY_TYPE MATCHES "static")
    set(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt_static)
  else()
    set(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt_shared)
  endif()

  # SuperLU_MT libraries
  list(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  # BLAS libraries
  if(BLAS_FOUND)
    list(APPEND SUNLINSOLSLUMT_LIBS ${BLAS_LIBRARIES})
  endif(BLAS_FOUND)

  foreach(example_tuple ${IDAS_examples_SUPERLUMT})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_args)
    list(GET example_tuple 2 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})
      # example source files
      add_executable(${example} ${example}.c)

      # folder to organize targets in an IDE
      set_target_properties(${example} PROPERTIES FOLDER "Examples")

      # libraries to link against
      target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLSLUMT_LIBS})
    endif()

    # check if example args are provided and set the test name
    if("${example_args}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    endif()

    # Do not include SuperLUMT examples in testing when the indextype is int64_t.
    # Answer files were generated with int32_t and minor differences in output
    # occur causing a false positive when testing. These tests can be re-enabled
    # when type specific answer files are added.
    if(SUNDIALS_INDEX_SIZE MATCHES "32")

      # add example to regression tests
      sundials_add_test(${test_name} ${example}
        TEST_ARGS ${example_args}
        ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        ANSWER_FILE ${test_name}.out
        EXAMPLE_TYPE ${example_type})

    endif(SUNDIALS_INDEX_SIZE MATCHES "32")

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

    # install example source and .out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)
    endif()

  endforeach(example_tuple ${IDAS_examples_SUPERLUMT})

endif(SUPERLUMT_FOUND)



if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "IDAS")
  set(SOLVER_LIB "sundials_idas")

  examples2string(IDAS_examples EXAMPLES)

  if(LAPACK_FOUND)
    examples2string(IDAS_examples_BL EXAMPLES_BL)
  else()
    set(EXAMPLES_BL "")
  endif()

  if(KLU_FOUND)
    examples2string(IDAS_examples_KLU EXAMPLES_KLU)
  else()
    set(EXAMPLES_KLU "")
  endif()

  if(SUPERLUMT_FOUND)
    examples2string(IDAS_examples_SUPERLUMT EXAMPLES_SLUMT)
    if (SUNDIALS_SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
      set(THREAD_LIBRARY_SLUMT ${CMAKE_THREAD_LIBS_INIT})
    else()
      set(THREAD_LIBRARY_SLUMT "")
    endif()
  else()
    set(EXAMPLES_SLUMTU "")
    set(THREAD_LIBRARY_SLUMT "")
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/idas/serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/idas/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/idas/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/idas/serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/idas/serial
      RENAME Makefile
      )
  endif(UNIX)

  # add test_install target
  sundials_add_test_install(idas serial EXECUTABLE idasRoberts_dns)

endif(EXAMPLES_INSTALL)
