# Copyright 2019 Mike Dev
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt
#
# NOTE: CMake support for Boost.Chrono is currently experimental at best
#       and the interface is likely to change in the future

cmake_minimum_required( VERSION 3.5 )
project( BoostChrono LANGUAGES CXX)

add_library( boost_chrono
    src/chrono.cpp
    src/process_cpu_clocks.cpp
    src/thread_clock.cpp
)
add_library( Boost::chrono ALIAS boost_chrono )

target_include_directories( boost_chrono PUBLIC include )
target_compile_definitions( boost_chrono
    PUBLIC
        # NOTE:
        # We deactivate autolinking because cmake based builds don't need it
        # and we don't implement name mangling for the library file anyway.
        # Ususally the parent CMakeLists.txt file should already have globally defined BOOST_ALL_NO_LIB
        BOOST_CHRONO_NO_LIB
        $<$<STREQUAL:$<TARGET_PROPERTY:boost_chrono,TYPE>,SHARED_LIBRARY>:BOOST_CHRONO_DYN_LINK=1>
        $<$<STREQUAL:$<TARGET_PROPERTY:boost_chrono,TYPE>,STATIC_LIBRARY>:BOOST_CHRONO_STATIC_LINK=1>
)

target_link_libraries( boost_chrono
    PUBLIC
        Boost::assert
        Boost::config
        Boost::core
        Boost::integer
        Boost::move
        Boost::mpl
        Boost::predef
        Boost::ratio
        Boost::static_assert
        Boost::system
        Boost::throw_exception
        Boost::type_traits
        Boost::typeof
        Boost::utility
     PRIVATE
        $<$<STREQUAL:Windows,$<PLATFORM_ID>>:Boost::winapi>
)
