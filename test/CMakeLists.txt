# Prob tests
# Define the test generator executable
add_executable(generate_tests prob/generate_tests.cpp)


# Linking libraries, adjust according to actual dependencies required
target_include_directories(generate_tests PRIVATE ${boost_SOURCE_DIR})
# Add a custom target to run the test generator
# Define a custom command that always runs and uses the TIMESTAMP to force a re-run
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/always_rebuild
    COMMAND ${CMAKE_COMMAND} -E echo "Running generate_tests to produce distribution tests"
    COMMAND generate_tests
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/always_rebuild
    DEPENDS generate_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating distribution tests"
    VERBATIM
)

# Custom target that always builds when the project is built
add_custom_target(run_generate_tests ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/always_rebuild
    COMMENT "Running distribution test generator"
)

# Ensure generate_tests runs at build time
add_dependencies(run_generate_tests generate_tests)


# Discover all .hpp files that define tests and generate corresponding .cpp test files
file(GLOB_RECURSE HEADER_TEST_FILES "${CMAKE_CURRENT_SOURCE_DIR}/prob/*.hpp")
foreach(HEADER_FILE ${HEADER_TEST_FILES})
    # Replace slashes in the header file path to create a valid target name
    string(REPLACE "/" "_" SAFE_HEADER_NAME ${HEADER_FILE})
    string(REGEX REPLACE "/" "_" GENERATED_TEST_FILE ${SAFE_HEADER_NAME})

    # Generate each test file with a custom command
    add_custom_command(
        OUTPUT ${GENERATED_TEST_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "Generating test for ${HEADER_FILE}"
        COMMAND generate_tests ${HEADER_FILE} ${N_TESTS}
        DEPENDS generate_tests ${HEADER_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating ${GENERATED_TEST_FILE} from ${HEADER_FILE}"
        VERBATIM
    )

    # Define a custom target for each generated test file, using the safe name
    add_custom_target("${GENERATED_TEST_FILE}_target" ALL
        DEPENDS ${GENERATED_TEST_FILE}
        COMMENT "Generating test target for ${GENERATED_TEST_FILE}"
    )
endforeach()

# Unit Tests
add_library(test_ad_pch INTERFACE)
target_precompile_headers(test_ad_pch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/unit/math/test_ad.hpp)
target_include_directories(test_ad_pch INTERFACE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR})
target_link_libraries(test_ad_pch INTERFACE benchmark::benchmark TBB::tbb Eigen3::Eigen sundials_kinsol sundials_cvodes)
# Function to add a test target for each directory containing C++ source files
function(add_gtest_grouped_test test_directory)
  # Recursively find all directories under the specified test directory
  file(GLOB_RECURSE ALL_CPP_FILES "${test_directory}/*.cpp")

  # Get unique directories from file list
  set(UNIQUE_DIRS "")
  foreach(cpp_file IN LISTS ALL_CPP_FILES)
    get_filename_component(dir_path "${cpp_file}" DIRECTORY)
    list(APPEND UNIQUE_DIRS "${dir_path}")
  endforeach()
  list(REMOVE_DUPLICATES UNIQUE_DIRS)

  # Create a test target for each unique directory with CPP files
  foreach(dir IN LISTS UNIQUE_DIRS)
    file(GLOB CPP_FILES_IN_DIR "${dir}/*test.cpp")

    # Generate a safe test target name from the directory path
    string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_DIR ${dir})
    string(REPLACE "/" "_" TEST_TARGET_NAME ${RELATIVE_DIR})
    if(TEST_TARGET_ONLY_NAME STREQUAL "")
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
        target_link_libraries(${TEST_TARGET_NAME} test_ad_pch gtest_main benchmark::benchmark TBB::tbb Eigen3::Eigen sundials_kinsol sundials_cvodes)

        # Register the test with CTest
        add_test(NAME ${TEST_TARGET_NAME} COMMAND ${TEST_TARGET_NAME})
      endif()
    else()
      message(WARNING "Target already exists: ${TEST_TARGET_NAME}")
    endif()
  endforeach()
endfunction()

# Enable CTest to recognize the tests
enable_testing()

# Entry point to start the recursive addition of tests
message(STATUS "Building tests...")
add_gtest_grouped_test(${CMAKE_CURRENT_SOURCE_DIR}/unit)
add_gtest_grouped_test(${CMAKE_CURRENT_SOURCE_DIR}/prob)
message(STATUS "Building tests: DONE")
