#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 Dmitry Arkhipov (grisumbras@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

# Files
file(GLOB_RECURSE BOOST_URL_TESTS_FILES CONFIGURE_DEPENDS *.cpp *.hpp)
if (DEFINED BOOST_URL_TEST_FLAGS AND NOT BOOST_URL_TEST_FLAGS STREQUAL "")
    set_source_files_properties(${BOOST_URL_TESTS_FILES} PROPERTIES COMPILE_FLAGS ${BOOST_URL_TEST_FLAGS})
    set_target_properties(boost_url PROPERTIES COMPILE_FLAGS ${BOOST_URL_TEST_FLAGS})
endif()
set_property(SOURCE doc_grammar.cpp PROPERTY COMPILE_FLAGS "")
set_property(SOURCE doc_3_urls.cpp PROPERTY COMPILE_FLAGS "")
list(APPEND BOOST_URL_TESTS_FILES CMakeLists.txt Jamfile)
set(EXAMPLE_FILES ../../example/router/impl/matches.cpp ../../example/router/detail/impl/router.cpp)

# Test target
add_executable(boost_url_unit_tests EXCLUDE_FROM_ALL ${BOOST_URL_TESTS_FILES} ${SUITE_FILES} ${EXAMPLE_FILES})
target_include_directories(boost_url_unit_tests PRIVATE . ../../extra ../../example/router)
target_link_libraries(boost_url_unit_tests PUBLIC Boost::url)
foreach (BOOST_URL_UNIT_TEST_LIBRARY ${BOOST_URL_UNIT_TEST_LIBRARIES})
    target_link_libraries(boost_url_unit_tests PUBLIC Boost::${BOOST_URL_UNIT_TEST_LIBRARY})
endforeach ()

# Compile options
target_compile_options(boost_url_unit_tests PUBLIC $<$<CXX_COMPILER_ID:GNU>:-Wno-unused-but-set-variable>)
target_compile_options(boost_url_unit_tests PUBLIC $<$<CXX_COMPILER_ID:GNU>:-Wno-unused-function>)
if (TARGET boost_filesystem AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_target_property(FS_IS_IMPORTED boost_filesystem IMPORTED)
    target_compile_options(boost_filesystem PUBLIC $<$<BOOL:${FS_IS_IMPORTED}>:-Wno-error=restrict>)
endif()

# Folders
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${BOOST_URL_TESTS_FILES})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../../extra PREFIX "_extra" FILES ${SUITE_FILES})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../../example/router PREFIX "_router" FILES ${EXAMPLE_FILES})

# CTest target
add_test(NAME boost_url_unit_tests COMMAND boost_url_unit_tests)
add_dependencies(tests boost_url_unit_tests)

