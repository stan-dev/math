include_directories(${CMAKE_HOME_DIRECTORY})

include(CTest)

enable_testing()

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -j4 -C ${CMAKE_CFG_INTDIR})

##################################################
# Parser tests
##################################################
add_executable(parser_ parser.cpp)
set_property(TARGET parser_ PROPERTY CXX_STANDARD ${CXX_STD})
target_link_libraries(parser_ parser boost ${link_flags})
if (MSVC)
    target_compile_options(parser_ PRIVATE /source-charset:utf-8 /bigobj)
elseif (USE_ASAN OR USE_UBSAN)
    target_compile_options(parser_ PRIVATE ${compile_flags})
endif ()
add_test(NAME parser_ COMMAND parser_)

add_executable(parser_api parser_api.cpp)
set_property(TARGET parser_api PROPERTY CXX_STANDARD ${CXX_STD})
target_link_libraries(parser_api parser boost ${link_flags})
if (MSVC)
    target_compile_options(parser_api PRIVATE /source-charset:utf-8 /bigobj)
elseif (USE_ASAN OR USE_UBSAN)
    target_compile_options(parser_api PRIVATE ${compile_flags})
endif ()
add_test(NAME parser_api COMMAND parser_api)

add_executable(
    compile_tests
    compile_tests_main.cpp
    compile_attribute.cpp
    compile_seq_attribute.cpp
    compile_or_attribute.cpp
    compile_combining_groups.cpp
    compile_all_t.cpp
)
set_property(TARGET compile_tests PROPERTY CXX_STANDARD ${CXX_STD})
target_link_libraries(compile_tests parser boost)

macro(add_test_executable name)
    add_executable(${name} ${name}.cpp)
    set_property(TARGET ${name} PROPERTY CXX_STANDARD ${CXX_STD})
    target_link_libraries(${name} parser boost ${link_flags})
    if (MSVC)
        target_compile_options(${name} PRIVATE /source-charset:utf-8 /bigobj)
    elseif (USE_ASAN OR USE_UBSAN)
        target_compile_options(${name} PRIVATE ${compile_flags})
    endif ()
    add_test(NAME ${name} COMMAND ${name})
endmacro()

add_test_executable(github_issues)
add_test_executable(search)
add_test_executable(split)
add_test_executable(replace)
add_test_executable(transform_replace)
add_test_executable(hl)
add_test_executable(aggr_tuple_assignment)
add_test_executable(parser_perm)
add_test_executable(parser_attributes)
add_test_executable(parser_lazy_params)
add_test_executable(parser_if_switch)
add_test_executable(parser_rule)
add_test_executable(parser_rule_with_params)
add_test_executable(parser_action)
add_test_executable(parser_action_with_params)
add_test_executable(parser_symbol_table)
add_test_executable(parser_quoted_string)
add_test_executable(tracing)
add_test_executable(parse_empty)
add_test_executable(tuple_aggregate)
add_test_executable(class_type)
add_test_executable(case_fold_generated)
add_test_executable(no_case)
add_test_executable(merge_separate)
add_test_executable(parse_coords_new)
add_test_executable(parser_seq_permutations_1)
add_test_executable(parser_seq_permutations_2)
add_test_executable(parser_or_permutations_1)
add_test_executable(parser_or_permutations_2)

if (MSVC)
    add_executable(vs_output_tracing tracing.cpp)
    set_property(TARGET vs_output_tracing PROPERTY CXX_STANDARD ${CXX_STD})
    target_link_libraries(vs_output_tracing parser boost ${link_flags})
    target_compile_options(vs_output_tracing PRIVATE /source-charset:utf-8 /bigobj)
    target_compile_definitions(vs_output_tracing PRIVATE BOOST_PARSER_TRACE_TO_VS_OUTPUT)
    add_test(NAME vs_output_tracing COMMAND vs_output_tracing)
endif()

set(COVERAGE_TRACE_COMMAND
    lcov --capture
    -output-file ${PROJECT_BINARY_DIR}/coverage.info
    --directory ${PROJECT_BINARY_DIR}
    --include ${PROJECT_SOURCE_DIR}/include/*
)

set(COVERAGE_HTML_COMMAND
    genhtml --legend -f -q
    ${PROJECT_BINARY_DIR}/coverage.info
    --prefix ${PROJECT_SOURCE_DIR}
    --output-directory ${PROJECT_BINARY_DIR}/coverage_html
)

add_custom_target(
    coverage
        COMMAND ${COVERAGE_TRACE_COMMAND}
        COMMAND ${COVERAGE_HTML_COMMAND}
        COMMENT "Generating coverage report"
        VERBATIM
)
