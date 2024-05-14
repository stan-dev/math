# Unit Tests
# Function to add a test target for each directory containing C++ source files
function(add_gtest_grouped_test test_directory target_pch)
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
        message(STATUS "Adding grouped test for directory: ${dir} as ${TEST_TARGET_NAME}")
        add_executable(${TEST_TARGET_NAME} ${CPP_FILES_IN_DIR})
        # Configure target properties such as include directories and linked libraries
        target_include_directories(${TEST_TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${boost_SOURCE_DIR})
        target_link_libraries(${TEST_TARGET_NAME} target_pch gtest_main benchmark::benchmark TBB::tbb Eigen3::Eigen sundials_kinsol sundials_cvodes)

        # Register the test with CTest
        add_test(NAME ${TEST_TARGET_NAME} COMMAND ${TEST_TARGET_NAME})
      endif()
    else()
      message(WARNING "Target already exists: ${TEST_TARGET_NAME}")
    endif()
  endforeach()

  # Create a single target for all tests in the test_directory
  string(REPLACE "/" "_" SINGLE_TARGET_NAME "${test_directory}_test")
  if(ALL_TEST_FILES)
    message(STATUS "Adding single test target for all tests in: ${test_directory} as ${SINGLE_TARGET_NAME}")
    add_executable(${SINGLE_TARGET_NAME} ${ALL_TEST_FILES})
    target_include_directories(${SINGLE_TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${boost_SOURCE_DIR})
    target_link_libraries(${SINGLE_TARGET_NAME} target_pch gtest_main benchmark::benchmark TBB::tbb Eigen3::Eigen sundials_kinsol sundials_cvodes)

    # Register the combined test with CTest
    add_test(${test_directory}_test ${SINGLE_TARGET_HPP} COMMAND ${SINGLE_TARGET_NAME})
  endif()
endfunction()


# Entry point to start the recursive addition of tests
message(STATUS "Building tests...")
# Enable CTest to recognize the tests
enable_testing()
#add_subdirectory(prob)
add_subdirectory(unit)
message(STATUS "Building tests: DONE")
