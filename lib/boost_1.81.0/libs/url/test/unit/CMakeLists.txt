#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 Dmitry Arkhipov (grisumbras@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

if(NOT TARGET tests)
    add_custom_target(tests)
    set_property(TARGET tests PROPERTY FOLDER Dependencies)
endif()

set(BOOST_URL_TESTS_FILES
    CMakeLists.txt
    Jamfile
    test_rule.hpp
    authority_view.cpp
    doc_3_urls.cpp
    doc_grammar.cpp
    decode_view.cpp
    encode.cpp
    encoding_opts.cpp
    error.cpp
    error_types.cpp
    grammar.cpp
    host_type.cpp
    ignore_case.cpp
    ipv4_address.cpp
    ipv6_address.cpp
    natvis.cpp
    optional.cpp
    param.cpp
    params_base.cpp
    params_encoded_view.cpp
    params_view.cpp
    params_encoded_base.cpp
    params_encoded_ref.cpp
    params_ref.cpp
    parse.cpp
    parse_path.cpp
    parse_query.cpp
    pct_string_view.cpp
    scheme.cpp
    segments_base.cpp
    segments_encoded_base.cpp
    segments_encoded_ref.cpp
    segments_encoded_view.cpp
    segments_ref.cpp
    segments_view.cpp
    snippets.cpp
    static_url.cpp
    string_view.cpp
    url.cpp
    url_base.cpp
    url_view.cpp
    url_view_base.cpp
    urls.cpp
    variant.cpp
    grammar/alnum_chars.cpp
    grammar/alpha_chars.cpp
    grammar/charset.cpp
    grammar/ci_string.cpp
    grammar/dec_octet_rule.cpp
    grammar/delim_rule.cpp
    grammar/digit_chars.cpp
    grammar/grammar_error.cpp
    grammar/grammar_parse.cpp
    grammar/hexdig_chars.cpp
    grammar/literal_rule.cpp
    grammar/lut_chars.cpp
    grammar/not_empty_rule.cpp
    grammar/optional_rule.cpp
    grammar/range_rule.cpp
    grammar/recycled.cpp
    grammar/string_token.cpp
    grammar/string_view_base.cpp
    grammar/token_rule.cpp
    grammar/tuple_rule.cpp
    grammar/type_traits.cpp
    grammar/unsigned_rule.cpp
    grammar/variant_rule.cpp
    grammar/vchars.cpp
    rfc/absolute_uri_rule.cpp
    rfc/authority_rule.cpp
    rfc/gen_delim_chars.cpp
    rfc/ipv4_address_rule.cpp
    rfc/ipv6_address_rule.cpp
    rfc/origin_form_rule.cpp
    rfc/pchars.cpp
    rfc/pct_encoded_rule.cpp
    rfc/query_rule.cpp
    rfc/relative_ref_rule.cpp
    rfc/reserved_chars.cpp
    rfc/sub_delim_chars.cpp
    rfc/unreserved_chars.cpp
    rfc/uri_rule.cpp
    rfc/uri_reference_rule.cpp
    )

set(SUITE_FILES ../../extra/test_main.cpp ../../extra/test_suite.hpp)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${BOOST_URL_TESTS_FILES})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../../extra PREFIX "_extra" FILES ${SUITE_FILES})
add_executable(boost_url_tests ${BOOST_URL_TESTS_FILES} ${SUITE_FILES})
target_include_directories(boost_url_tests PRIVATE . ../../extra)
# The include dependencies are found in the CMakeLists.txt
# of the root project directory.
# See: BOOST_URL_UNIT_TEST_LIBRARIES
target_link_libraries(boost_url_tests PRIVATE
    Boost::url
    Boost::container
    Boost::filesystem
    Boost::unordered)
add_test(NAME boost_url_tests COMMAND boost_url_tests)
add_dependencies(boost_url_all_tests boost_url_tests)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(boost_url_tests PUBLIC -Wno-error=unused-but-set-variable)
    if (TARGET boost_filesystem AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        get_target_property(FS_IS_IMPORTED boost_filesystem IMPORTED)
        if (FS_IS_IMPORTED)
            target_compile_options(boost_filesystem PUBLIC -Wno-error=restrict)
        endif()
    endif()
endif()