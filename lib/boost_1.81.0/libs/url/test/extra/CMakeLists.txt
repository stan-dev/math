#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 DMitry Arkhipov (grisumbras@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

set(TEST_FILES
    Jamfile
    test_suite.cpp
    )

file(GLOB_RECURSE SUITE_FILES CONFIGURE_DEPENDS
    ../../extra/*.cpp
    ../../extra/*.hpp
    ../../extra/*.ipp
)

#set(SUITE_FILES ../../extra/test_main.cpp ../../extra/test_suite.hpp)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${TEST_FILES})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../../extra PREFIX "_extra" FILES ${SUITE_FILES})
add_executable(boost_url_extra ${TEST_FILES} ${SUITE_FILES})
target_include_directories(boost_url_extra PRIVATE . ../../extra)
if (BOOST_URL_FIND_PACKAGE_BOOST)
    target_link_libraries(boost_url_extra PRIVATE Boost::headers)
else()
    target_link_libraries(boost_url_extra PRIVATE
        Boost::align
        Boost::config
        Boost::core
        Boost::optional
        Boost::type_traits
        Boost::system
        Boost::variant2)
endif()
add_test(NAME boost_url_extra COMMAND boost_url_extra)
add_dependencies(boost_url_all_tests boost_url_extra)
