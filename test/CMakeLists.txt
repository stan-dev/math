# Stan's vendored gtest_main initializes the MPI cluster and includes Stan
# headers even when MPI is disabled. Use the test flags instead of
# GoogleTest's internal warning policy for its own sources.
target_link_libraries(gtest_main PRIVATE stan_test_flags)
set_property(TARGET gtest_main PROPERTY COMPILE_FLAGS "")
add_library(stan_test_options INTERFACE)
target_link_libraries(stan_test_options INTERFACE stan_test_flags gtest_main)

# Set out to TRUE when path contains every STAN_TEST_NAME_FILTERS substring.
function(stan_matches_name_filters path out)
  foreach(filter IN LISTS STAN_TEST_NAME_FILTERS)
    string(FIND "${path}" "${filter}" match_index)
    if(match_index EQUAL -1)
      set(${out} FALSE PARENT_SCOPE)
      return()
    endif()
  endforeach()
  set(${out} TRUE PARENT_SCOPE)
endfunction()

# Every test executable is created here so each gets the same dependencies.
function(stan_add_test_binary target)
  # Xcode needs a real source alongside object-library references.
  add_executable(${target} "${PROJECT_SOURCE_DIR}/cmake/empty.cpp" ${ARGN})
  if(STAN_MPI)
    # A static archive can discard the unreferenced registration initializer.
    target_sources(${target} PRIVATE $<TARGET_OBJECTS:stan_mpi>)
  endif()
  stan_configure_test(${target} stan_test_options "")
  if(WIN32)
    add_custom_command(TARGET ${target} POST_BUILD
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
        "$<TARGET_FILE:TBB::tbb>" "$<TARGET_FILE_DIR:${target}>"
      VERBATIM)
  endif()
endfunction()

# Add the CTest test `name` for one source, labeled with ARGN and the
# source's repository path. Its objects (`name`_objects) are shared with the
# executables of stan_add_test_group.
function(stan_add_gtest name source repo_relative pch)
  add_library(${name}_objects OBJECT EXCLUDE_FROM_ALL "${source}")
  stan_configure_test(${name}_objects stan_test_options "${pch}")
  get_source_file_property(generator_target "${source}" STAN_GENERATOR_TARGET)
  if(generator_target)
    # Every consuming executable needs generation before object compilation.
    add_dependencies(${name}_objects ${generator_target})
  endif()
  stan_add_test_binary(${name} "$<TARGET_OBJECTS:${name}_objects>")
  if(STAN_MPI)
    set(mpi_ranks 1)
    if(repo_relative MATCHES "mpi_")
      set(mpi_ranks 2)
    endif()
    add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${mpi_ranks}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}> ${MPIEXEC_POSTFLAGS})
  else()
    add_test(NAME ${name} COMMAND ${name})
  endif()
  set_tests_properties(${name} PROPERTIES
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    LABELS "${ARGN};${repo_relative}")
  if(STAN_SELECTED_TEST_REGEX AND repo_relative MATCHES "${STAN_SELECTED_TEST_REGEX}")
    set_property(GLOBAL APPEND PROPERTY STAN_SELECTED_TESTS ${name})
  endif()
endfunction()

# Build the tests in ARGN with `prefix`_subtests, or as one executable with
# `prefix`_all_tests. The stan_build: label lets runTests.py find the group.
function(stan_add_test_group prefix)
  add_custom_target(${prefix}_subtests DEPENDS ${ARGN})
  set_property(TEST ${ARGN} APPEND PROPERTY LABELS "stan_build:${prefix}_subtests")
  set(objects "")
  foreach(test IN LISTS ARGN)
    list(APPEND objects "$<TARGET_OBJECTS:${test}_objects>")
  endforeach()
  stan_add_test_binary(${prefix}_all_tests ${objects})
  set_target_properties(${prefix}_all_tests PROPERTIES EXCLUDE_FROM_ALL TRUE)
endfunction()

# Add a test for each source in SOURCES, or for each *test.cpp below
# DIRECTORY, named after its path below REPO_PATH. Add a group for every
# directory level, starting with PREFIX. Return the selected tests in
# `PREFIX`_selected_tests. Sources use PCH, or `LAYER_PCH_PREFIX`_<layer>_pch
# when it exists for their first directory. With AD_PCH, read each source to
# choose AD_PCH or no PCH.
function(add_gtest_grouped_test)
  cmake_parse_arguments(PARSE_ARGV 0 arg ""
    "PREFIX;REPO_PATH;DIRECTORY;PCH;LAYER_PCH_PREFIX;AD_PCH" "SOURCES")
  set(sources ${arg_SOURCES})
  if(NOT sources)
    file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${arg_DIRECTORY}/*test.cpp")
  endif()
  set(directory_keys "")
  set(scanned_sources "")
  foreach(src IN LISTS sources)
    file(RELATIVE_PATH relative_source "${arg_DIRECTORY}" "${src}")
    set(repo_relative "${arg_REPO_PATH}/${relative_source}")
    if(NOT STAN_OPENCL AND repo_relative MATCHES "/opencl/")
      continue()
    endif()
    if(NOT STAN_MPI AND repo_relative MATCHES "/[^/]*mpi[^/]*test\\.cpp$")
      continue()
    endif()
    string(REPLACE "/" "_" name "${repo_relative}")
    string(REGEX REPLACE "\\.cpp$" "" name "${name}")
    set(pch "${arg_PCH}")
    if(arg_LAYER_PCH_PREFIX)
      string(REGEX MATCH "^[^/]+" layer "${relative_source}")
      if(TARGET ${arg_LAYER_PCH_PREFIX}_${layer}_pch)
        set(pch ${arg_LAYER_PCH_PREFIX}_${layer}_pch)
      endif()
    endif()
    if(arg_AD_PCH)
      # Precompile the AD helpers only for sources that already include them;
      # their global aliases must not leak into unrelated tests. A source
      # that sets STAN_/EIGEN_ macros before including Stan/Eigen gets no PCH.
      file(STRINGS "${src}" directives
        REGEX "^#(include <test/unit/math/test_ad\\.hpp>|[ \t]*(define|undef)[ \t]+(STAN_|EIGEN_))")
      list(APPEND scanned_sources "${src}")
      if(directives MATCHES "(^|;)#[ \t]*(define|undef)")
        set(pch "")
      elseif(directives AND TARGET ${arg_AD_PCH})
        set(pch ${arg_AD_PCH})
      endif()
    endif()
    get_filename_component(relative_dir "${relative_source}" DIRECTORY)
    string(REPLACE "/" ";" labels "${relative_dir}")
    stan_add_gtest(${name} "${src}" "${repo_relative}" "${pch}" ${arg_PREFIX} ${labels})

    stan_matches_name_filters("${repo_relative}" selected)
    if(NOT selected)
      continue()
    endif()
    set(key "${relative_dir}")
    while(TRUE)
      if(key STREQUAL "")
        set(key ROOT)
      endif()
      string(REPLACE "/" "_" key "${key}")
      if(NOT DEFINED subtests_${key})
        list(APPEND directory_keys "${key}")
      endif()
      list(APPEND subtests_${key} "${name}")
      if(key STREQUAL "ROOT")
        break()
      endif()
      # The next ancestor comes from the path, not the underscore target name.
      get_filename_component(relative_dir "${relative_dir}" DIRECTORY)
      set(key "${relative_dir}")
    endwhile()
  endforeach()
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${scanned_sources})
  foreach(key IN LISTS directory_keys)
    if(key STREQUAL "ROOT")
      stan_add_test_group(${arg_PREFIX} ${subtests_ROOT})
    else()
      stan_add_test_group(${arg_PREFIX}_${key} ${subtests_${key}})
    endif()
  endforeach()
  set(${arg_PREFIX}_selected_tests "${subtests_ROOT}" PARENT_SCOPE)
endfunction()

add_subdirectory(unit)
if(STAN_BUILD_PROB_TESTS)
  add_subdirectory(prob)
endif()
if(STAN_BUILD_EXPRESSION_TESTS)
  set(STAN_EXPRESSION_TEST_DIR "${PROJECT_BINARY_DIR}/expression-sources"
    CACHE PATH "Directory containing generated expression tests")
  file(GLOB expression_sources CONFIGURE_DEPENDS "${STAN_EXPRESSION_TEST_DIR}/*test.cpp")
  if(NOT expression_sources)
    message(FATAL_ERROR "Generate expression tests before configuring STAN_BUILD_EXPRESSION_TESTS")
  endif()
  # The expression helper has no include guard; precompile its heavy guarded
  # dependencies so each generated source can still include the helper once.
  stan_add_pch(expression_pch stan_test_options <stan/math/mix.hpp> <gtest/gtest.h>)
  add_gtest_grouped_test(PREFIX expressions REPO_PATH test/expressions
    DIRECTORY "${STAN_EXPRESSION_TEST_DIR}" PCH expression_pch SOURCES ${expression_sources})
endif()

# Makefile generators build separate top-level targets one at a time; one
# aggregate builds a selection of individual tests with every job.
get_property(selected_tests GLOBAL PROPERTY STAN_SELECTED_TESTS)
if(selected_tests)
  add_custom_target(stan_selected_tests DEPENDS ${selected_tests})
endif()

if(STAN_TEST_HEADERS)
  # A PCH would conceal missing includes. Force each header into an otherwise
  # empty translation unit, and do not inherit the tests' precompiled headers.
  file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/stan/*.hpp")
  add_custom_target(test-headers)
  foreach(header IN LISTS headers)
    file(RELATIVE_PATH relative_header "${PROJECT_SOURCE_DIR}" "${header}")
    if(NOT STAN_OPENCL AND relative_header MATCHES "/opencl/")
      continue()
    endif()
    if(NOT STAN_MPI AND relative_header MATCHES "/mpi_[^/]+\\.hpp$")
      continue()
    endif()
    string(REPLACE "/" "_" name "${relative_header}")
    string(REGEX REPLACE "\\.hpp$" "-test" name "${name}")
    add_library(${name} OBJECT EXCLUDE_FROM_ALL "${PROJECT_SOURCE_DIR}/cmake/empty.cpp")
    target_link_libraries(${name} PRIVATE stanmath)
    if(MSVC)
      target_compile_options(${name} PRIVATE "/FI${header}")
    else()
      target_compile_options(${name} PRIVATE -O0 -include "${header}" -Wunused-local-typedefs)
    endif()
    add_dependencies(test-headers ${name})
  endforeach()
endif()
