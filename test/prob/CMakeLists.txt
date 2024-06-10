# Compile one test manually so that we can reuse precompile headers
add_executable(prob_pch
  ${CMAKE_CURRENT_SOURCE_DIR}/dummy_precompile_target.cpp)
target_compile_options(prob_pch PRIVATE -O1)
target_include_directories(prob_pch PRIVATE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR})
target_precompile_headers(prob_pch PRIVATE
    "$<$<COMPILE_LANGUAGE:CXX>:<gtest/gtest.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:C>:<stddef.h$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<stan/math/mix.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<stan/math/rev.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/test_fixture_ccdf_log.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/test_fixture_cdf_log.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/test_fixture_cdf.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/test_fixture_distr.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/utility.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<gtest/gtest.h$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<stdexcept$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<tuple$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<type_traits$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_int_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_int_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_int_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_real_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_real_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_real_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_int_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_int_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_int_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_int_int_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_int_int_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_int_int_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_int_int_int_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_int_int_int_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_int_int_int_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_real_real_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_real_real_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_real_real_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_int_int_int_int_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_int_int_int_int_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_int_int_int_int_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_real_real_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_real_real_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_real_real_int_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_v_real_real_real_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_ffv_real_real_real_real_real_pch.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/prob/args/arg_generated_vv_real_real_real_real_real_pch.hpp$<ANGLE-R>>"
)

target_link_libraries(prob_pch
  gtest_main benchmark::benchmark TBB::tbb
  Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
  sundials_idas sundials_nvecserial)

if (STAN_OPENCL)
  target_link_libraries(prob_pch OpenCL::OpenCL)
  target_include_directories(prob_pch PRIVATE ${CLCPP_SOURCE_DIR}/include)
endif()
target_compile_options(prob_pch PRIVATE -O1)

# Prob tests
# Define the test generator executable
add_executable(generate_tests generate_tests.cpp)


# Linking libraries, adjust according to actual dependencies required
target_include_directories(generate_tests PRIVATE ${boost_SOURCE_DIR})
# Add a custom target to run the test generator
# Define a custom command that always runs and uses the TIMESTAMP to force a re-run
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/always_rebuild1
    COMMAND ${CMAKE_COMMAND} -E echo "Running generate_tests to produce distribution tests"
    COMMAND generate_tests ${CMAKE_CURRENT_SOURCE_DIR} 300
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/always_rebuild1
    DEPENDS generate_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating distribution tests"
    VERBATIM
)

# Custom target that always builds when the project is built
add_custom_target(run_generate_tests ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/always_rebuild1
    COMMENT "Running distribution test generator"
)

# Ensure generate_tests runs at build time
add_dependencies(run_generate_tests generate_tests)


add_gtest_grouped_test(${CMAKE_CURRENT_SOURCE_DIR} prob_math prob_pch)
