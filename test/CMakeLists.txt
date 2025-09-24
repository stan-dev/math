enable_testing()
# Unit Tests
## Get all properties that cmake supports
execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
## Convert command output into a CMake list
STRING(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
STRING(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)

function(print_target_properties tgt)
    if(NOT TARGET ${tgt})
      message("There is no target named '${tgt}'")
      return()
    endif()

    foreach (prop ${CMAKE_PROPERTY_LIST})
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
        get_target_property(propval ${tgt} ${prop})
        if (propval)
            message ("${tgt} ${prop} = ${propval}")
        endif()
    endforeach(prop)
endfunction(print_target_properties)

# Compile one test manually so that we can reuse precompile headers
add_executable(dae_test ${CMAKE_CURRENT_SOURCE_DIR}/unit/math/rev/functor/pph_dae_typed_test.cpp)
target_include_directories(dae_test PRIVATE ${CMAKE_SOURCE_DIR} ${CLCPP_SOURCE_DIR}/include)
target_precompile_headers(dae_test PRIVATE
    "$<$<COMPILE_LANGUAGE:CXX>:<gtest/gtest.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:C>:<stddef.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<test/unit/math/test_ad.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/fwd.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/mix.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/rev.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/prim.hpp$<ANGLE-R>>")

target_link_libraries(dae_test
  gtest_main benchmark::benchmark TBB::tbb
  Eigen3::Eigen sundials_kinsol_static sundials_cvodes_static
  sundials_idas_static sundials_nvecserial_static
  Boost::math Boost::numeric_odeint Boost::lexical_cast Boost::optional)
if (STAN_OPENCL)
  target_link_libraries(dae_test OpenCL::OpenCL)
  target_include_directories(dae_test PRIVATE ${CLCPP_SOURCE_DIR}/include)
endif()
if (STAN_MPI)
  target_link_libraries(dae_test Boost::mpi MPI::MPI_CXX)
endif()

define_property(TARGET PROPERTY SUB_TESTS BRIEF_DOCS "Subtests for a test target")

# Adds a test target for
# 1. Each directory under the specified test directory
# 2. A single target for each test in the test directory
# 3. A single target for all tests in the test directory
# The test target is named after the directory path, with / replaced by _
# @param test_directory The directory containing the tests
# @param folder_test_name The name of the test for the whole folder, will have _test appended to it
# @param target_pch The precompiled header target to reuse
function(add_gtest_grouped_test test_directory folder_test_name target_pch)
  # Collect all test sources (*test.cpp) once.
  file(GLOB_RECURSE ALL_TEST_SOURCES "${test_directory}/*test.cpp")

  if (NOT ALL_TEST_SOURCES)
    message(WARNING "No *test.cpp files found under: ${test_directory}")
    return()
  endif()

  # Common libs for all test targets.
  set(_COMMON_TEST_LIBS
      gtest_main
      benchmark::benchmark
      TBB::tbb
      Eigen3::Eigen
      sundials_kinsol_static
      sundials_cvodes_static
      sundials_idas_static
      sundials_nvecserial_static
      Boost::math Boost::numeric_odeint Boost::lexical_cast Boost::optional)

  # ---------------------------------------------
  # 1) Per-file subtests and map them to folders.
  #    We accumulate, for every directory (and its
  #    ancestors up to test_directory), the list of
  #    subtest executables that live under it.
  # ---------------------------------------------
  set(_DIR_KEYS "")  # stores the set of directory keys we’ll later materialize as *_subtests targets

  foreach(src IN LISTS ALL_TEST_SOURCES)
    # Unique name from repo-root-relative path.
    string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_FILE "${src}")
    string(REPLACE "/" "_" SUB_TEST_TARGET_NAME "${RELATIVE_FILE}")
    string(REPLACE ".cpp" "" SUB_TEST_TARGET_NAME "${SUB_TEST_TARGET_NAME}")

    if (NOT TARGET ${SUB_TEST_TARGET_NAME})
      add_executable(${SUB_TEST_TARGET_NAME} "${src}")
      target_compile_options(${SUB_TEST_TARGET_NAME} PRIVATE -O1 -Wno-misleading-indentation)
      target_include_directories(${SUB_TEST_TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
      if (TARGET ${target_pch})
        target_precompile_headers(${SUB_TEST_TARGET_NAME} REUSE_FROM ${target_pch})
      endif()
      target_link_libraries(${SUB_TEST_TARGET_NAME} PRIVATE ${_COMMON_TEST_LIBS})
      if (STAN_OPENCL)
        target_link_libraries(${SUB_TEST_TARGET_NAME} PRIVATE OpenCL::OpenCL)
        target_include_directories(${SUB_TEST_TARGET_NAME} PRIVATE ${CLCPP_SOURCE_DIR}/include)
      endif()
      if (STAN_MPI)
        target_link_libraries(${SUB_TEST_TARGET_NAME} PRIVATE Boost::mpi MPI::MPI_CXX)
      endif()

      add_test(NAME ${SUB_TEST_TARGET_NAME} COMMAND ${SUB_TEST_TARGET_NAME})
      set_property(TEST ${SUB_TEST_TARGET_NAME} PROPERTY LABELS "${folder_test_name}")
    endif()

    # Assign this subtest to every ancestor directory (recursive) up to test_directory.
    get_filename_component(_cur_dir "${src}" DIRECTORY)
    while (TRUE)
      # REL_DIR == "." or "" means root (test_directory itself).
      file(RELATIVE_PATH _rel_dir "${test_directory}" "${_cur_dir}")
      if (_rel_dir STREQUAL "." OR _rel_dir STREQUAL "")
        set(_key "ROOT")
      else()
        string(REPLACE "/" "_" _key "${_rel_dir}")
      endif()

      # Track keys and append this subtest into the list for this key.
      list(APPEND _DIR_KEYS "${_key}")
      set(_varname "SUBTESTS_${_key}")
      set(${_varname} "${${_varname}};${SUB_TEST_TARGET_NAME}")

      if (_cur_dir STREQUAL "${test_directory}")
        break()
      endif()
      get_filename_component(_cur_dir "${_cur_dir}" DIRECTORY)
    endwhile()
  endforeach()

  list(REMOVE_DUPLICATES _DIR_KEYS)

  # ---------------------------------------------
  # 2) Per-directory grouped executables (recursive)
  #    and per-directory *_subtests meta targets.
  # ---------------------------------------------
  # Build set of unique directories we saw while scanning sources (for _all_tests).
  set(UNIQUE_DIRS "")
  foreach(src IN LISTS ALL_TEST_SOURCES)
    get_filename_component(_d "${src}" DIRECTORY)
    while (TRUE)
      list(APPEND UNIQUE_DIRS "${_d}")
      if (_d STREQUAL "${test_directory}")
        break()
      endif()
      get_filename_component(_d "${_d}" DIRECTORY)
    endwhile()
  endforeach()
  list(REMOVE_DUPLICATES UNIQUE_DIRS)

  set(ALL_GROUP_TARGETS "")
  foreach(dir IN LISTS UNIQUE_DIRS)
    file(GLOB_RECURSE DIR_TESTS "${dir}/*test.cpp")
    if (NOT DIR_TESTS)
      continue()
    endif()

    file(RELATIVE_PATH REL_DIR "${test_directory}" "${dir}")
    if (REL_DIR STREQUAL "." OR REL_DIR STREQUAL "")
      set(GROUP_TARGET "${folder_test_name}_all_tests")
      set(GROUP_LABEL  "${folder_test_name}")
      set(_key "ROOT")
    else()
      string(REPLACE "/" "_" REL_NAME "${REL_DIR}")
      set(GROUP_TARGET "${folder_test_name}_${REL_NAME}_all_tests")
      set(GROUP_LABEL  "${folder_test_name};${REL_NAME}")
      set(_key "${REL_NAME}")
    endif()

    if (NOT TARGET ${GROUP_TARGET})
      add_executable(${GROUP_TARGET} ${DIR_TESTS})
      target_compile_options(${GROUP_TARGET} PRIVATE -O1 -Wno-misleading-indentation)
      target_include_directories(${GROUP_TARGET} PRIVATE ${CMAKE_SOURCE_DIR})
      if (TARGET ${target_pch})
        target_precompile_headers(${GROUP_TARGET} REUSE_FROM ${target_pch})
      endif()
      target_link_libraries(${GROUP_TARGET} PRIVATE ${_COMMON_TEST_LIBS})
      if (STAN_OPENCL)
        target_link_libraries(${GROUP_TARGET} PRIVATE OpenCL::OpenCL)
        target_include_directories(${GROUP_TARGET} PRIVATE ${CLCPP_SOURCE_DIR}/include)
      endif()
      if (STAN_MPI)
        target_link_libraries(${GROUP_TARGET} PRIVATE Boost::mpi MPI::MPI_CXX)
      endif()

      add_test(NAME ${GROUP_TARGET} COMMAND ${GROUP_TARGET})
      set_property(TEST ${GROUP_TARGET} PROPERTY LABELS "${GROUP_LABEL}")
      list(APPEND ALL_GROUP_TARGETS ${GROUP_TARGET})
    endif()

    # ---- Per-directory *_subtests meta target (recursive) ----
    # Pull the list of subtests for this directory key and create a meta target that depends on them.
    set(_varname "SUBTESTS_${_key}")
    if (DEFINED ${_varname})
      # Deduplicate the list for cleanliness.
      list(REMOVE_DUPLICATES ${_varname})

      if (_key STREQUAL "ROOT")
        set(SUBTESTS_TARGET "${folder_test_name}_subtests")
      else()
        set(SUBTESTS_TARGET "${folder_test_name}_${_key}_subtests")
      endif()

      if (NOT TARGET ${SUBTESTS_TARGET})
        add_custom_target(${SUBTESTS_TARGET})
        add_dependencies(${SUBTESTS_TARGET} ${${_varname}})
        # Optional: stash the list as a property for inspection.
        set_property(TARGET ${SUBTESTS_TARGET} PROPERTY SUBTESTS ${${_varname}})
      endif()
    endif()
  endforeach()

  # ---------------------------------------------
  # 3) Global convenience meta targets
  # ---------------------------------------------
  add_custom_target(${folder_test_name}_groups)
  add_dependencies(${folder_test_name}_groups ${ALL_GROUP_TARGETS})

  add_custom_target(${folder_test_name}_all)
  add_dependencies(${folder_test_name}_all ${folder_test_name}_groups)

endfunction()

# Entry point to start the recursive addition of tests
message(STATUS "Building tests...")
# Enable CTest to recognize the tests
enable_testing()
add_subdirectory(prob)
add_subdirectory(unit)
message(STATUS "Building tests: DONE")

if (STAN_TEST_HEADERS)
message(STATUS, "Building Test Headers")
# Detect OS
if(WIN32)
  set(DEV_NULL "NUL")
  set(IS_UCRT true) # Adjust this condition as necessary
  if(IS_UCRT)
    set(UCRT_NULL_FLAG "-S")
  endif()
else()
  set(DEV_NULL "/dev/null")
endif()

# Compile one test manually so that we can reuse precompile headers
add_executable(dae_header_pch ${CMAKE_CURRENT_SOURCE_DIR}/unit/math/rev/functor/pph_dae_typed_test.cpp)
target_include_directories(dae_header_pch PRIVATE ${CMAKE_SOURCE_DIR} ${CLCPP_SOURCE_DIR}/include)
target_precompile_headers(dae_header_pch PRIVATE
    "$<$<COMPILE_LANGUAGE:CXX>:<gtest/gtest.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:C>:<stddef.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<test/unit/math/test_ad.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/fwd.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/mix.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/rev.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/prim.hpp$<ANGLE-R>>")

target_link_libraries(dae_header_pch
  gtest_main benchmark::benchmark TBB::tbb
  Eigen3::Eigen sundials_kinsol_static sundials_cvodes_static
  sundials_idas_static sundials_nvecserial_static Boost::math Boost::numeric_odeint Boost::lexical_cast Boost::optional)
if (STAN_OPENCL)
  target_link_libraries(dae_header_pch OpenCL::OpenCL)
  target_include_directories(dae_header_pch PRIVATE ${CLCPP_SOURCE_DIR}/include)
endif()
if (STAN_MPI)
  target_link_libraries(dae_header_pch Boost::mpi MPI::MPI_CXX)
endif()
target_compile_options(dae_header_pch PUBLIC -O0 -Wunused-local-typedefs)



# Find header files
file(GLOB_RECURSE HEADER_FILES "${CMAKE_SOURCE_DIR}/stan/*.hpp")

# Create dummy file if it doesn't exist
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test)
set(DUMMY_CPP "${CMAKE_BINARY_DIR}/test/dummy.cpp")
if(NOT EXISTS ${DUMMY_CPP})
  file(WRITE ${DUMMY_CPP} "int main() {return 0;}\n")
endif()

# Custom target for testing headers
add_custom_target(test-headers)
# Create targets for each header file
foreach(HEADER ${HEADER_FILES})
  get_filename_component(HEADER_NAME ${HEADER} NAME_WE)
  set(HEADER_TEST "${HEADER_NAME}-test")
  string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_FILE ${HEADER})
  string(REPLACE "/" "_" TMP_TEST_TARGET_NAME ${RELATIVE_FILE})
  string(REPLACE ".hpp" "-test" SUB_TEST_TARGET_NAME ${TMP_TEST_TARGET_NAME})

  add_executable(${SUB_TEST_TARGET_NAME} ${HEADER} ${DUMMY_CPP})
  set_target_properties(${SUB_TEST_TARGET_NAME} PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE
  )
  target_compile_options(${SUB_TEST_TARGET_NAME} PUBLIC -O0 -Wunused-local-typedefs)
  target_include_directories(${SUB_TEST_TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
  target_precompile_headers(${SUB_TEST_TARGET_NAME} REUSE_FROM dae_header_pch)
  target_link_libraries(${SUB_TEST_TARGET_NAME}
    gtest_main benchmark::benchmark TBB::tbb
    Eigen3::Eigen sundials_kinsol_static sundials_cvodes_static
    sundials_idas_static sundials_nvecserial_static Boost::math Boost::numeric_odeint Boost::lexical_cast Boost::optional)
  if (STAN_OPENCL)
    target_link_libraries(${SUB_TEST_TARGET_NAME} OpenCL::OpenCL)
    target_include_directories(${SUB_TEST_TARGET_NAME} PRIVATE ${CLCPP_SOURCE_DIR}/include)
  endif()
  if (STAN_MPI)
    target_link_libraries(${SUB_TEST_TARGET_NAME} Boost::mpi MPI::MPI_CXX)
  endif()

  add_custom_command(
    TARGET ${SUB_TEST_TARGET_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Testing header ${HEADER}"
#    COMMAND rm -f $<TARGET_FILE:${SUB_TEST_TARGET_NAME}>
#    COMMENT "Deleting target output for ${SUB_TEST_TARGET_NAME}"
  )

  add_dependencies(test-headers ${SUB_TEST_TARGET_NAME})
endforeach()
endif()
