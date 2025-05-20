
# Common utilities
add_library(boost_redis_project_options INTERFACE)
target_link_libraries(boost_redis_project_options INTERFACE boost_redis)
if (MSVC)
  target_compile_options(boost_redis_project_options INTERFACE /bigobj)
  target_compile_definitions(boost_redis_project_options INTERFACE _WIN32_WINNT=0x0601)
endif()

add_library(boost_redis_src STATIC boost_redis.cpp)
target_compile_features(boost_redis_src PRIVATE cxx_std_17)
target_link_libraries(boost_redis_src PRIVATE boost_redis_project_options)

# Test utils
add_library(boost_redis_tests_common STATIC common.cpp)
target_compile_features(boost_redis_tests_common PRIVATE cxx_std_17)
target_link_libraries(boost_redis_tests_common PRIVATE boost_redis_project_options)

macro(make_test TEST_NAME STANDARD)
  set(EXE_NAME "boost_redis_${TEST_NAME}")
  add_executable(${EXE_NAME} ${TEST_NAME}.cpp)
  target_link_libraries(${EXE_NAME} PRIVATE
    boost_redis_src
    boost_redis_tests_common
    boost_redis_project_options
    Boost::unit_test_framework
  )
  target_compile_features(${EXE_NAME} PRIVATE cxx_std_${STANDARD})
  add_test(${EXE_NAME} ${EXE_NAME})
endmacro()

make_test(test_conn_quit 17)
# TODO: Configure a Redis server with TLS in the CI and reenable this test.
#make_test(test_conn_tls 17)
make_test(test_low_level 17)
make_test(test_conn_exec_retry 17)
make_test(test_conn_exec_error 17)
make_test(test_request 17)
make_test(test_run 17)
make_test(test_low_level_sync_sans_io 17)
make_test(test_conn_check_health 17)

make_test(test_conn_exec 20)
make_test(test_conn_push 20)
make_test(test_conn_reconnect 20)
make_test(test_conn_exec_cancel 20)
make_test(test_conn_exec_cancel2 20)
make_test(test_conn_echo_stress 20)
make_test(test_conn_run_cancel 20)
make_test(test_issue_50 20)
make_test(test_issue_181 17)

# Coverage
set(
  COVERAGE_TRACE_COMMAND
  lcov --capture
  -output-file "${PROJECT_BINARY_DIR}/coverage.info"
  --directory "${PROJECT_BINARY_DIR}"
  --include "${PROJECT_SOURCE_DIR}/include/*"
)

set(
  COVERAGE_HTML_COMMAND
  genhtml --legend -f -q
  "${PROJECT_BINARY_DIR}/coverage.info"
  --prefix "${PROJECT_SOURCE_DIR}"
  --output-directory "${PROJECT_BINARY_DIR}/coverage_html"
)

add_custom_target(
  coverage
  COMMAND ${COVERAGE_TRACE_COMMAND}
  COMMAND ${COVERAGE_HTML_COMMAND}
  COMMENT "Generating coverage report"
  VERBATIM
)
