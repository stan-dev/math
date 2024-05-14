
add_library(test_ad_pch2 INTERFACE)
target_precompile_headers(test_ad_pch2 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../unit/math/test_ad.hpp)
target_include_directories(test_ad_pch2 INTERFACE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR})
target_link_libraries(test_ad_pch2 INTERFACE benchmark::benchmark TBB::tbb Eigen3::Eigen sundials_kinsol sundials_cvodes)

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

add_gtest_grouped_test(${CMAKE_CURRENT_SOURCE_DIR} test_ad_pch2)
