#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/static_string
#

# Custom target used by the boost super-project
if(NOT TARGET tests)
        add_custom_target(tests)
        set_property(TARGET tests PROPERTY FOLDER _deps)
endif()

set(BOOST_STATIC_STRING_TESTS_FILES
        CMakeLists.txt
        Jamfile
        constexpr_tests.hpp
        compile_fail.hpp
        static_string.cpp
        )

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${BOOST_STATIC_STRING_TESTS_FILES})
add_executable(boost_static_string_tests ${BOOST_STATIC_STRING_TESTS_FILES})

# The include dependencies are found in the CMakeLists.txt
# of the root project directory.
# See: BOOST_STATIC_STRING_UNIT_TEST_LIBRARIES
target_link_libraries(boost_static_string_tests PRIVATE Boost::static_string)
add_test(NAME boost_static_string_tests COMMAND boost_static_string_tests)
add_dependencies(tests boost_static_string_tests)
