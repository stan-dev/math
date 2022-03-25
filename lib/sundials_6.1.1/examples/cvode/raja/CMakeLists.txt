# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for CVODE RAJA examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
  "cvAdvDiff_kry_raja.cpp\;develop"
  )

# Add source directory to include directories
include_directories(.)

if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
  set(_fused_stubs sundials_cvode_fused_stubs)
endif()

if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_lang CUDA)
else()
  set(_lang CXX)
endif()

if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  set(OTHER_LIBS OpenMP::OpenMP_CXX)
endif()

# Add the build and install targets for each CVODE example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  set_source_files_properties(${example} PROPERTIES LANGUAGE ${_lang})

  # add example source files
  add_executable(${example_target} ${example})

  set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

  # libraries to link against
  target_link_libraries(${example_target}
    PRIVATE
    sundials_cvode
    sundials_nvecraja
    ${_fused_stubs}
    ${EXTRA_LINK_LIBS}
    ${OTHER_LIBS})

  # add example to regression tests
  sundials_add_test(${example_target} ${example_target}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example_target}.out
    EXAMPLE_TYPE ${example_type})

endforeach(example_tuple ${CVODE_examples})


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(_fused_stubs_target cvode_fused_stubs)
    set(LIBS "-lsundials_cvode_fused_stubs ${LIBS}")
  endif()

  if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
    set(EXAMPLES_FIND_PACKAGE "find_package(OpenMP REQUIRED)\n")
    set(_openmp_target OpenMP::OpenMP_CXX)
  endif()

  if(RAJA_NEEDS_THREADS)
    set(EXAMPLES_FIND_PACKAGE "${EXAMPLES_FIND_PACKAGE}find_package(Threads REQUIRED)\n")
  endif()

  sundials_install_examples(cvode CVODE_examples
    CMAKE_TEMPLATE
      cmakelists_${_lang}_ex.in
    MAKE_TEMPLATE
      makefile_serial_RAJA_ex.in
    SOLVER_LIBRARY
      sundials_cvode
    SUNDIALS_TARGETS
      cvode
      nvecraja
      ${_fused_stubs_target}
    OTHER_TARGETS
      ${_openmp_target}
    DESTINATION
      cvode/raja
    EXTRA_FILES
      README
    TEST_INSTALL
      raja
  )

endif()
