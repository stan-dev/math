# ---------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the ARKODE parallel F90 examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
if(SUNDIALS_PRECISION MATCHES "double")
  # The fark_heat2D only supports double precision
  set(FARKODE_examples
    "fark_heat2D\;1\;2\;develop"
    )
  # Auxiliary files to install
  set(FARKODE_extras
    plot_heat2D.py
    )
endif()

if(MPI_Fortran_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link against
set(ARKODE_LIB sundials_arkode)
set(NVECP_LIB sundials_nvecparallel)
set(FNVECP_LIB sundials_fnvecparallel)

# Only static FCMIX libraries are available
set(FARKODE_LIB sundials_farkode${_STATIC_LIB_SUFFIX})

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${FARKODE_LIB} ${ARKODE_LIB} ${FNVECP_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each ARKODE example
foreach(example_tuple ${FARKODE_examples})

  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # example source files
  add_executable(${example} ${example}.f90)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")
  set_target_properties(${example} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # add example to regression tests
  sundials_add_test(${example} ${example}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS})

  if(NOT MPI_Fortran_COMPILER)
    target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
  endif()

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.f90 ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel)
  endif()

endforeach(example_tuple ${FARKODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel)

  # Install the extra files
  foreach(extrafile ${FARKODE_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "ARKODE")
  set(SOLVER_LIB "sundials_arkode")
  set(SOLVER_FLIB "sundials_farkode")

  examples2string(FARKODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_F90_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/F90_parallel/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/F90_parallel/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_F90_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/F90_parallel/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/F90_parallel/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/F90_parallel
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(arkode F90_parallel)

endif()
