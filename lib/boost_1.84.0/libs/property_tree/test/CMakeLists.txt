#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 Richard Hodges (hodges.r@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

file(GLOB_RECURSE BOOST_PROPERTY_TREE_TESTS_FILES CONFIGURE_DEPENDS Jamfile.v2 *.cpp *.hpp)
list(FILTER BOOST_PROPERTY_TREE_TESTS_FILES EXCLUDE REGEX "^${CMAKE_CURRENT_SOURCE_DIR}/cmake_install_test/.*$")

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "" FILES ${BOOST_PROPERTY_TREE_TESTS_FILES})

if(NOT TARGET tests)
    add_custom_target(tests)
endif()

macro(PTREE_TEST name)
    #message(STATUS "ptree_test ${PROJECT_NAME}-${name} SRCS ${ARGN}")
    add_executable("${PROJECT_NAME}-${name}" ${ARGN})
    target_include_directories("${PROJECT_NAME}-${name}" PRIVATE .)
    target_link_libraries("${PROJECT_NAME}-${name}" PRIVATE Boost::property_tree Boost::serialization)
    add_dependencies(tests "${PROJECT_NAME}-${name}")
    add_test(NAME "${PROJECT_NAME}-${name}" COMMAND "${PROJECT_NAME}-${name}")
endmacro()

PTREE_TEST(test-multi test_multi_module1.cpp test_multi_module2.cpp)
PTREE_TEST(test-property-tree test_property_tree.cpp test_property_tree.hpp)
PTREE_TEST(test-info-parser test_info_parser.cpp)
PTREE_TEST(test-rapidxml test_rapidxml.cpp)
PTREE_TEST(test-json-parser test_json_parser.cpp)
PTREE_TEST(test-json-parser2 test_json_parser2.cpp)
PTREE_TEST(test-ini-parser test_ini_parser.cpp)
PTREE_TEST(test-xml-parser-rapidxml test_xml_parser_rapidxml.cpp)
