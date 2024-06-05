# Unit Tests

# Compile one test manually so that we can reuse precompile headers
add_executable(dae_test ${CMAKE_CURRENT_SOURCE_DIR}/unit/math/rev/functor/pph_dae_typed_test.cpp)
target_include_directories(dae_test PRIVATE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR} ${CLCPP_SOURCE_DIR}/include)
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
  Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
  sundials_idas sundials_nvecserial $<$<BOOL:STAN_OPENCL>:OpenCL::OpenCL>)


# Function to add a test target for each directory containing C++ source files
function(add_gtest_grouped_test test_directory folder_test_name target_pch)
  # Recursively find all directories under the specified test directory
  file(GLOB_RECURSE ALL_CPP_FILES "${test_directory}/*.cpp")

  # Get unique directories from file list
  set(UNIQUE_DIRS "")
  foreach(cpp_file IN LISTS ALL_CPP_FILES)
    get_filename_component(dir_path "${cpp_file}" DIRECTORY)
    list(APPEND UNIQUE_DIRS "${dir_path}")
  endforeach()
  list(REMOVE_DUPLICATES UNIQUE_DIRS)

  # List to accumulate all test files for the single target
  set(ALL_TEST_FILES "")
  set(ALL_TEST_TARGETS "")
  # Create a test target for each unique directory with CPP files
  foreach(dir IN LISTS UNIQUE_DIRS)
    file(GLOB CPP_FILES_IN_DIR "${dir}/*test.cpp")
    list(APPEND ALL_TEST_FILES ${CPP_FILES_IN_DIR})  # Add files to global list

    # Generate a safe test target name from the directory path
    string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_DIR ${dir})
    string(REPLACE "/" "_" TEST_TARGET_NAME ${RELATIVE_DIR})
    if(TEST_TARGET_NAME STREQUAL "")
      set(TEST_TARGET_NAME "root")
    endif()

    # Check if target already exists to prevent redefinition
    if (NOT TARGET ${TEST_TARGET_NAME})
      # Add executable and related commands only if there are source files
      if (CPP_FILES_IN_DIR)
        foreach(cpp_file IN LISTS CPP_FILES_IN_DIR)
          string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_FILE ${cpp_file})
          string(REPLACE "/" "_" TMP_TEST_TARGET_NAME ${RELATIVE_FILE})
          string(REPLACE ".cpp" "" SUB_TEST_TARGET_NAME ${TMP_TEST_TARGET_NAME})
          add_executable(${SUB_TEST_TARGET_NAME} ${cpp_file})
          target_include_directories(${SUB_TEST_TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}
            ${boost_SOURCE_DIR} ${CLCPP_SOURCE_DIR}/include)
          target_precompile_headers(${SUB_TEST_TARGET_NAME} REUSE_FROM ${target_pch})
          target_link_libraries(${SUB_TEST_TARGET_NAME}
            gtest_main benchmark::benchmark TBB::tbb
            Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
            sundials_idas sundials_nvecserial $<$<BOOL:STAN_OPENCL>:OpenCL::OpenCL>)
          # Register the test with CTest
          add_test(NAME ${SUB_TEST_TARGET_NAME} COMMAND ${SUB_TEST_TARGET_NAME})
        endforeach()
        message(STATUS "Adding grouped test for directory: ${dir} as ${TEST_TARGET_NAME}")
        add_executable(${TEST_TARGET_NAME} ${CPP_FILES_IN_DIR})
        # Configure target properties such as include directories and linked libraries
        target_include_directories(${TEST_TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}
         ${boost_SOURCE_DIR} ${CLCPP_SOURCE_DIR}/include)
        target_compile_options(${TEST_TARGET_NAME} PUBLIC -O1)
        target_precompile_headers(${TEST_TARGET_NAME} REUSE_FROM ${target_pch})
        target_link_libraries(${TEST_TARGET_NAME}
          gtest_main benchmark::benchmark TBB::tbb
          Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
          sundials_idas sundials_nvecserial $<$<BOOL:STAN_OPENCL>:OpenCL::OpenCL>)
          # Register the test with CTest
        add_test(NAME ${TEST_TARGET_NAME} COMMAND ${TEST_TARGET_NAME})
        list(APPEND ALL_TEST_TARGETS ${TEST_TARGET_NAME})
      endif()
    else()
      message(WARNING "Target already exists: ${TEST_TARGET_NAME}")
    endif()
  endforeach()

  # Create a single target for all tests in the test_directory
  if(ALL_TEST_FILES)
    message(STATUS "Adding single test target for all tests in: ${test_directory} as ${folder_test_name}_test")
    add_executable(${folder_test_name}_test ${ALL_TEST_FILES})
    target_compile_options(${folder_test_name}_test PUBLIC -O1)
    target_include_directories(${folder_test_name}_test PRIVATE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR} ${CLCPP_SOURCE_DIR}/include)
    target_precompile_headers(${folder_test_name}_test REUSE_FROM ${target_pch})
    target_link_libraries(${folder_test_name}_test
          gtest_main benchmark::benchmark TBB::tbb
          Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
          sundials_idas sundials_nvecserial $<$<BOOL:STAN_OPENCL>:OpenCL::OpenCL>)
    # Register the combined test with CTest
    add_test(${test_directory}_test ${SINGLE_TARGET_HPP} COMMAND ${folder_test_name}_test)
  endif()
  if (ALL_TEST_TARGETS)
  add_custom_target(${folder_test_name}_subtests)
  add_dependencies(${folder_test_name}_subtests ${ALL_TEST_TARGETS})
  message(STATUS "Adding ${folder_test_name}_subtests to build")
    message(STATUS "${ALL_TEST_TARGETS} ")
  endif()
endfunction()


# Entry point to start the recursive addition of tests
message(STATUS "Building tests...")
# Enable CTest to recognize the tests
enable_testing()
add_subdirectory(prob)
add_subdirectory(unit)
message(STATUS "Building tests: DONE")


