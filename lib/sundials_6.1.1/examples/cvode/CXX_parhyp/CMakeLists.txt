# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for CVODE C++ parhyp examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the
# type is develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cv_heat2D_hypre_ls\;--np 2 2\;1\;4\;develop"
  "cv_heat2D_hypre_pfmg\;--np 2 2\;1\;4\;develop"
  )

# Auxiliary files to install
set(CVODE_extras
  plot_heat2D_p.py
  )

if(MPI_CXX_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

  # disable C++ extensions (for known wrappers)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -DLAM_BUILDING")
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Specify libraries to link against
set(CVODE_LIB sundials_cvode)
set(NVECP_LIB sundials_nvecparallel)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_LIB sundials_cvode_fused_stubs)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECP_LIB} ${EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    add_executable(${example} ${example}.cpp)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} SUNDIALS::HYPRE)

    if(NOT MPI_CXX_COMPILER)
      target_link_libraries(${example} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARIES})
    endif()
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
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.cpp ${example_out}
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/CXX_parhyp)
  endif()

endforeach(example_tuple ${CVODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/CXX_parhyp)

  # Install the extra files
  foreach(extrafile ${CVODE_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/CXX_parhyp)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  examples2string(CVODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parhyp_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/CXX_parhyp/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/cvode/CXX_parhyp/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/CXX_parhyp
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parhyp_CXX_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/CXX_parhyp/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/CXX_parhyp/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/cvode/CXX_parhyp
      RENAME Makefile
      )
  endif(UNIX)

  # add test_install target
  sundials_add_test_install(cvode CXX_parhyp)

endif(EXAMPLES_INSTALL)
