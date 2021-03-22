# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
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
# CMakeLists.txt file for the FCVODE parallel examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(FCVODE_examples
  "fcvDiag_non_p\;1\;2\;develop"
  "fcvDiag_kry_bbd_p\;1\;2\;develop"
  "fcvDiag_kry_p\;1\;2\;develop"
  )

if(MPI_Fortran_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link againstx
set(FNVECP_LIB sundials_fnvecparallel)

# Only static FCMIX libraries are available
set(FCVODE_LIB sundials_fcvode${_STATIC_LIB_SUFFIX})

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  list(APPEND FCVODE_LIB sundials_cvode_fused_stubs)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${FCVODE_LIB} ${FNVECP_LIB})

# Add the build and install targets for each example
foreach(example_tuple ${FCVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # example source files
  add_executable(${example} ${example}.f)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  sundials_add_test(${example} ${example}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS} ${EXTRA_LINK_LIBS})

  if(NOT MPI_Fortran_COMPILER)
    target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  endif()

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.f ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_parallel)
  endif()

endforeach(example_tuple ${FCVODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_parallel)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  set(SOLVER_FLIB "sundials_fcvode")
  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  examples2string(FCVODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_F77_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_parallel/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_parallel
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_F77_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/fcmix_parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/fcmix_parallel
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(cvode fcmix_parallel)

endif()
