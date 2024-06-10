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
  sundials_idas sundials_nvecserial)
if (STAN_OPENCL)
  target_link_libraries(dae_test OpenCL::OpenCL)
  target_include_directories(dae_test PRIVATE ${CLCPP_SOURCE_DIR}/include)
endif()
if (STAN_MPI)
  target_link_libraries(dae_test Boost::mpi)
endif()


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
            ${boost_SOURCE_DIR})
          target_precompile_headers(${SUB_TEST_TARGET_NAME} REUSE_FROM ${target_pch})
          target_link_libraries(${SUB_TEST_TARGET_NAME}
            gtest_main benchmark::benchmark TBB::tbb
            Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
            sundials_idas sundials_nvecserial)
          if (STAN_OPENCL)
            target_link_libraries(${SUB_TEST_TARGET_NAME} OpenCL::OpenCL)
            target_include_directories(${SUB_TEST_TARGET_NAME} PRIVATE ${CLCPP_SOURCE_DIR}/include)
          endif()
          if (STAN_MPI)
            target_link_libraries(${SUB_TEST_TARGET_NAME} Boost::mpi)
          endif()
          # Register the test with CTest
          add_test(NAME ${SUB_TEST_TARGET_NAME} COMMAND ${SUB_TEST_TARGET_NAME})
        endforeach()
        message(STATUS "Adding grouped test ${TEST_TARGET_NAME} for ${dir}")
        add_executable(${TEST_TARGET_NAME} ${CPP_FILES_IN_DIR})
        # Configure target properties such as include directories and linked libraries
        target_include_directories(${TEST_TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}
         ${boost_SOURCE_DIR})
        target_compile_options(${TEST_TARGET_NAME} PUBLIC -O1)
        target_precompile_headers(${TEST_TARGET_NAME} REUSE_FROM ${target_pch})
        target_link_libraries(${TEST_TARGET_NAME}
          gtest_main benchmark::benchmark TBB::tbb
          Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
          sundials_idas sundials_nvecserial)
        if (STAN_OPENCL)
          target_link_libraries(${TEST_TARGET_NAME} OpenCL::OpenCL)
          target_include_directories(${TEST_TARGET_NAME} PRIVATE ${CLCPP_SOURCE_DIR}/include)
        endif()
        if (STAN_MPI)
          target_link_libraries(${TEST_TARGET_NAME} Boost::mpi)
        endif()

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
    message(STATUS "Adding single test target for all tests ${folder_test_name}_test for ${test_directory}")
    add_executable(${folder_test_name}_test ${ALL_TEST_FILES})
    target_compile_options(${folder_test_name}_test PUBLIC -O1)
    target_include_directories(${folder_test_name}_test PRIVATE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR})
    target_precompile_headers(${folder_test_name}_test REUSE_FROM ${target_pch})
    target_link_libraries(${folder_test_name}_test
          gtest_main benchmark::benchmark TBB::tbb
          Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
          sundials_idas sundials_nvecserial)
    if (STAN_OPENCL)
      target_link_libraries(${folder_test_name}_test OpenCL::OpenCL)
      target_include_directories(${folder_test_name}_test PRIVATE ${CLCPP_SOURCE_DIR}/include)
    endif()
    if (STAN_MPI)
      target_link_libraries(${folder_test_name}_test Boost::mpi)
    endif()

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

if (STAN_TEST_HEADERS)
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
target_include_directories(dae_header_pch PRIVATE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR} ${CLCPP_SOURCE_DIR}/include)
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
  Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
  sundials_idas sundials_nvecserial)
if (STAN_OPENCL)
  target_link_libraries(dae_header_pch OpenCL::OpenCL)
  target_include_directories(dae_header_pch PRIVATE ${CLCPP_SOURCE_DIR}/include)
endif()
if (STAN_MPI)
  target_link_libraries(dae_header_pch Boost::mpi)
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

  target_include_directories(${SUB_TEST_TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}
    ${boost_SOURCE_DIR})
  target_precompile_headers(${SUB_TEST_TARGET_NAME} REUSE_FROM dae_header_pch)
  target_link_libraries(${SUB_TEST_TARGET_NAME}
    gtest_main benchmark::benchmark TBB::tbb
    Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
    sundials_idas sundials_nvecserial)
  if (STAN_OPENCL)
    target_link_libraries(${SUB_TEST_TARGET_NAME} OpenCL::OpenCL)
    target_include_directories(${SUB_TEST_TARGET_NAME} PRIVATE ${CLCPP_SOURCE_DIR}/include)
  endif()
  if (STAN_MPI)
    target_link_libraries(${SUB_TEST_TARGET_NAME} Boost::mpi)
  endif()

  target_compile_options(${SUB_TEST_TARGET_NAME} PUBLIC -O0 -Wunused-local-typedefs)
  add_custom_command(
    TARGET ${SUB_TEST_TARGET_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Testing header ${HEADER}" > ${DEV_NULL}
  )

  add_dependencies(test-headers ${SUB_TEST_TARGET_NAME})
endforeach()
endif()
