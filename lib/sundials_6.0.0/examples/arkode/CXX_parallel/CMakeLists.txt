# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU,
#                Cody J. Balos @ LLNL
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
# -----------------------------------------------------------------
# CMakeLists.txt file for ARKODE C++ parallel examples
# -----------------------------------------------------------------

# macro to loop through a set of examples and build/install them
macro(build_examples examples_to_build lang)

  # Add the build and install targets for each example
  foreach(example_tuple ${${examples_to_build}})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_defines)
    list(GET example_tuple 2 example_args)
    list(GET example_tuple 3 number_of_nodes)
    list(GET example_tuple 4 number_of_tasks)
    list(GET example_tuple 5 example_type)
    list(GET example_tuple 6 example_precision)

    # extract the file name without extension
    get_filename_component(example_target ${example} NAME_WE)

    if(example_defines)
      set(example_target "${example_target}.${example_defines}")
    endif()

    # check if this example has already been added, only need to add
    # example source files once for testing with different inputs
    if(NOT TARGET ${example_target})

      set_source_files_properties(${example} PROPERTIES LANGUAGE ${lang})

      # example source files
      add_executable(${example_target} ${example})

      # folder to organize targets in an IDE
      set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

      # ensure the linker language is reset to CXX
      set_target_properties(${example_target} PROPERTIES LINKER_LANGUAGE CXX)

      # libraries to link against
      target_include_directories(${example_target}
        PRIVATE ${MPI_CXX_INCLUDE_DIRS})

      target_link_libraries(${example_target}
        PRIVATE ${OTHER_LIBS} ${SUNDIALS_LIBS} ${MPI_CXX_LIBRARIES})

      # compile definitions
      if(example_defines)
        target_compile_definitions(${example_target} PUBLIC ${example_defines})
      endif()

    endif()

    # check if example args are provided and set the test name
    if(example_args)
      string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
    else()
      set(test_name ${example_target})
    endif()

    # add example to regression tests
    sundials_add_test(${test_name} ${example_target}
      TEST_ARGS ${example_args}
      MPI_NPROCS ${number_of_tasks}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type}
      FLOAT_PRECISION ${example_precision})

    # Not all examples are added to the generated CMake file installed the
    # output and source files here
    if(EXAMPLES_INSTALL)

      # Find all .out files for this example
      file(GLOB example_out ${example_target}*.out)

      # Install example source and out files
      install(FILES ${example} ${example_out}
        DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/CXX_parallel)

    endif()

  endforeach()
endmacro()

# -----------------------------
# Examples to build and install
# -----------------------------

# Example lists are tuples "name\;compile defs\;args\;nodes\;tasks\;type\;float precision"
# where the type develop is for examples excluded from 'make test' in releases.

# List of headers to install (appended to below)
set(ARKODE_headers )

# List of additional files to install (appended to below)
set(ARKODE_extras)

# Additional libraries to link against (appended to below)
set(OTHER_LIBS ${EXTRA_LINK_LIBS})

# ------------
# MPI examples
# ------------

set(serial_examples
  "ark_heat2D_p.cpp\;\;--np 2 2\;1\;4\;develop\;default")
set(SUNDIALS_LIBS sundials_arkode sundials_nvecparallel)
build_examples(serial_examples CXX)

# Auxiliary files to install
list(APPEND ARKODE_extras plot_heat2D_p.py)

# ---------------------------
# ARKODE + CVODE MPI Examples
# ---------------------------

if(BUILD_CVODE)
  set(examples_cvode
    "ark_diffusion_reaction_p.cpp\;\;--np 2 2 --imex\;1\;4\;develop\;default"
    "ark_diffusion_reaction_p.cpp\;\;--np 2 2 --mri-arkstep\;1\;4\;develop\;default"
    "ark_diffusion_reaction_p.cpp\;\;--np 2 2 --mri-cvode-global\;1\;4\;develop\;default"
    "ark_diffusion_reaction_p.cpp\;\;--np 2 2 --mri-cvode-local\;1\;4\;develop\;default")
  set(SUNDIALS_LIBS sundials_arkode sundials_cvode sundials_nvecparallel sundials_nvecmpiplusx)
  build_examples(examples_cvode CXX)
endif()

# -------------
# RAJA Examples
# -------------

if(ENABLE_RAJA AND (SUNDIALS_PRECISION MATCHES "DOUBLE") AND
    (SUNDIALS_INDEX_SIZE MATCHES "32"))

  # Header files to install
  list(APPEND ARKODE_headers ark_brusselator1D.h)
  list(APPEND ARKODE_extras plot_brusselator1D.py)

  # All examples below link to RAJA
  set(OTHER_LIBS RAJA ${OTHER_LIBS})

  # If RAJA has HIP enabled, we have to link to HIP even if we dont use it
  if(RAJA_BACKENDS MATCHES "HIP")
    set(OTHER_LIBS hip::device ${OTHER_LIBS})
  endif()

  # If RAJA has OpenMP enabled, we have to link to OpenMP even if we dont use it
  if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
    set(OTHER_LIBS OpenMP::OpenMP_CXX ${OTHER_LIBS})
  endif()

  # -------------
  # RAJA + Serial
  # -------------

  set(serial_raja_examples
    "ark_brusselator1D_task_local_nls.cpp\;\;--monitor\;1\;4\;develop\;2"
    "ark_brusselator1D_task_local_nls.cpp\;\;--monitor --global-nls\;1\;4\;develop\;2"
    "ark_brusselator1D_task_local_nls.cpp\;\;--monitor --explicit --tf 1\;1\;4\;develop\;2")
  set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx)
  build_examples(serial_raja_examples CXX)

  # -----------
  # RAJA + CUDA
  # -----------

  if(RAJA_BACKENDS MATCHES "CUDA")

    if(BUILD_NVECTOR_CUDA)
      set(cuda_raja_examples
        "ark_brusselator1D_task_local_nls.cpp\;USE_CUDA_NVEC\;--monitor\;1\;4\;develop\;2"
        "ark_brusselator1D_task_local_nls.cpp\;USE_CUDAUVM_NVEC\;--monitor\;1\;4\;exclude\;2"
        "ark_brusselator1D_task_local_nls.cpp\;USE_CUDA_NVEC\;--monitor --global-nls\;1\;4\;develop\;2")
      set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx sundials_nveccuda)
      build_examples(cuda_raja_examples CUDA)
    endif()

    if(BUILD_NVECTOR_RAJA AND (SUNDIALS_RAJA_BACKENDS MATCHES "CUDA"))
      set(raja_raja_examples
        "ark_brusselator1D_task_local_nls.cpp\;USE_RAJA_NVEC\;--monitor\;1\;4\;exclude\;2")
      set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx sundials_nvecraja)
      build_examples(raja_raja_examples CUDA)
    endif()

    # ----------
    # RAJA + HIP
    # ----------

    if(BUILD_NVECTOR_HIP AND (RAJA_BACKENDS MATCHES "HIP"))
      set(hip_raja_examples
        "ark_brusselator1D_task_local_nls.cpp\;USE_HIP_NVEC\;--monitor\;1\;4\;exclude\;2")
      set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx sundials_nvechip)
      build_examples(hip_raja_examples CXX)
    endif()

    # ---------------------
    # RAJA + OpenMP Offload
    # ---------------------

    if(BUILD_NVECTOR_OPENMPDEV AND (RAJA_BACKENDS MATCHES "TARGET_OPENMP"))
      set(openmpdev_raja_examples
        "ark_brusselator1D_task_local_nls.cpp\;USE_OMPDEV_NVEC\;--monitor\;1\;4\;exclude\;2")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(SUNDIALS_LIBS sundials_arkode sundials_nvecmpiplusx sundials_nvecopenmpdev)
      build_examples(openmpdev_raja_examples CXX)
    endif()

  endif()
endif()


if(EXAMPLES_INSTALL)

  set(examples_to_install "${serial_examples}")
  set(_sundials_targets arkode nvecparallel)

  if(examples_cvode)
    list(APPEND examples_to_install "${examples_cvode}")
    list(APPEND _sundials_targets cvode nvecmpiplusx)
  endif()

  if(serial_raja_examples)
    list(APPEND examples_to_install "${serial_raja_examples}")
    list(APPEND _sundials_targets nvecmpiplusx)

    if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
      set(EXAMPLES_FIND_PACKAGE "find_package(OpenMP REQUIRED)\n")
    endif()

    if(RAJA_NEEDS_THREADS)
      set(EXAMPLES_FIND_PACKAGE "${EXAMPLES_FIND_PACKAGE}find_package(Threads REQUIRED)\n")
    endif()
  endif()

  sundials_install_examples(arkode examples_to_install
    CMAKE_TEMPLATE
      cmakelists_CXX_MPI_ex.in
    SUNDIALS_TARGETS
      ${_sundials_targets}
    OTHER_TARGETS
      ${OTHER_LIBS}
    DESTINATION
      arkode/CXX_parallel
    EXTRA_FILES
      ${ARKODE_headers}
      ${ARKODE_extras}
      README
    TEST_INSTALL
      CXX_parallel
  )

endif()
