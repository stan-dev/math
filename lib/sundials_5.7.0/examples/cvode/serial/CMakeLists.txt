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
# CMakeLists.txt file for CVODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cvAdvDiff_bnd\;\;develop"
  "cvDirectDemo_ls\;\;develop"
  "cvDiurnal_kry_bp\;\;develop"
  "cvDiurnal_kry\;\;develop"
  "cvDisc_dns\;\;develop"
  "cvKrylovDemo_ls\;\;develop"
  "cvKrylovDemo_ls\;1\;develop"
  "cvKrylovDemo_ls\;2\;develop"
  "cvKrylovDemo_prec\;\;develop"
  "cvParticle_dns\;\;develop"
  "cvPendulum_dns\;\;develop"
  "cvRoberts_dns\;\;"
  "cvRoberts_dns_uw\;\;develop"
  "cvRoberts_dns_negsol\;\;develop"
  "cvRoberts_dns_constraints\;\;develop"
  #cvAdvDiffReac_kry\;develop" # not released
  )

if(SUNDIALS_BUILD_WITH_MONITORING)
  list(APPEND CVODE_examples "cvKrylovDemo_ls\;0 1\;develop")
endif()

# Examples using LAPACK linear solvers
set(CVODE_examples_BL
  "cvAdvDiff_bndL\;\;develop"
  "cvRoberts_dnsL\;\;develop"
  )

# Examples using KLU linear solver
set(CVODE_examples_KLU
  "cvRoberts_klu\;\;develop"
  "cvRoberts_block_klu\;\;develop"
  )

# Examples using SuperLU_MT linear solver
set(CVODE_examples_SUPERLUMT
  "cvRoberts_sps\;\;develop"
  )

# Auxiliary files to install
set(CVODE_extras
  plot_cvParticle.py
  plot_cvPendulum.py
  )

# Specify libraries to link against
set(CVODE_LIB sundials_cvode)
set(NVECS_LIB sundials_nvecserial)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_LIB sundials_cvode_fused_stubs)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

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
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/serial)
  endif()

endforeach(example_tuple ${CVODE_examples})


# Add the build and install targets for each LAPACK example (if needed)
if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)

  # Sundials LAPACK linear solver modules
  set(SUNLINSOLLAPACK_LIBS
    sundials_sunlinsollapackband
    sundials_sunlinsollapackdense)

  # LAPACK libraries
  list(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  foreach(example_tuple ${CVODE_examples_BL})

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
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/serial)
    endif()

  endforeach(example_tuple ${CVODE_examples_BL})

endif()


# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  # Sundials KLU linear solver module
  set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu)

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS)

  foreach(example_tuple ${CVODE_examples_KLU})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_args)
    list(GET example_tuple 2 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})
      # add example source files
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
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/serial)
    endif()

  endforeach(example_tuple ${CVODE_examples_KLU})

endif()


# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  # Sundials SuperLU_MT linear solver module
  set(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt)

  # SuperLU_MT libraries
  list(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  foreach(example_tuple ${CVODE_examples_SUPERLUMT})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_args)
    list(GET example_tuple 2 example_type)

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example})
      # add example source files
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
    endif()

    # find all .out files for this example
    file(GLOB example_out ${example}*.out)

    # install example source and .out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/serial)
    endif()

  endforeach(example_tuple ${CVODE_examples_SUPERLUMT})

endif()


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/serial)

  # Install the extra files
  foreach(extrafile ${CVODE_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/serial)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  examples2string(CVODE_examples EXAMPLES)

  if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)
    examples2string(CVODE_examples_BL EXAMPLES_BL)
  else()
    set(EXAMPLES_BL "")
  endif()

  if(BUILD_SUNLINSOL_KLU)
    examples2string(CVODE_examples_KLU EXAMPLES_KLU)
  else()
    set(EXAMPLES_KLU "")
  endif()

  if(BUILD_SUNLINSOL_SUPERLUMT)
    examples2string(CVODE_examples_SUPERLUMT EXAMPLES_SLUMT)
    if (SUNDIALS_SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
      set(THREAD_LIBRARY_SLUMT ${CMAKE_THREAD_LIBS_INIT})
    else()
      set(THREAD_LIBRARY_SLUMT "")
    endif()
  else()
    set(EXAMPLES_SLUMT "")
    set(THREAD_LIBRARY_SLUMT "")
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/serial
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(cvode serial EXECUTABLE cvRoberts_dns)

endif()
