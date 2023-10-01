#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 DMitry Arkhipov (grisumbras@gmail.com)
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

set(F
    CMakeLists.txt
    main.cpp
    )

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "" FILES ${F})
add_executable(boost_wpt_tests ${F})
target_include_directories(boost_wpt_tests PRIVATE . ../../extra)
target_link_libraries(boost_wpt_tests PRIVATE Boost::url Boost::json)
add_test(NAME boost_url_tests COMMAND boost_url_tests)
add_dependencies(tests boost_url_tests)
