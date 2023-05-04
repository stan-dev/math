#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# This may be called from add_subdirectory, too
if (BOOST_MYSQL_USE_FIND_PACKAGE)
    find_package(Boost ${BOOST_MYSQL_VERSION} REQUIRED COMPONENTS unit_test_framework)
endif()

# Test common
add_library(
    boost_mysql_testing
    INTERFACE
)
target_include_directories(
    boost_mysql_testing
    INTERFACE
    common
)

if (BOOST_MYSQL_USE_FIND_PACKAGE)
    set(_TEST_LIBRARIES
        Boost::headers
        Boost::unit_test_framework
        Boost::disable_autolinking
    )
else()
    set(_TEST_LIBRARIES
        Boost::any
        Boost::optional
        Boost::unit_test_framework
        Boost::type_index
    )
endif()

target_link_libraries(
    boost_mysql_testing
    INTERFACE
    ${_TEST_LIBRARIES}
    Boost::mysql
    boost_mysql_asio
)

# Custom target tests; required by the Boost superproject
if(NOT TARGET tests)
    add_custom_target(tests)
endif()

# Unit testing
add_subdirectory(unit)

# Integration testing
if(BOOST_MYSQL_INTEGRATION_TESTS)
    add_subdirectory(integration)
endif()
