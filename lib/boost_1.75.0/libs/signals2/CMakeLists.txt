# Copyright 2019 Mike Dev
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt
#
# NOTE: CMake support for Boost.Signals2 is currently experimental at best
#       and the interface is likely to change in the future

cmake_minimum_required( VERSION 3.5...3.16 )
project( BoostSignals2 LANGUAGES CXX)

option( BOOST_SIGNALS2_INCLUDE_EXAMPLES "Add boost signals2 examples" OFF )

add_library( boost_signals2 INTERFACE )
add_library( Boost::signals2 ALIAS boost_signals2 )

target_include_directories( boost_signals2 INTERFACE include )

target_link_libraries( boost_signals2
    INTERFACE
        Boost::assert
        Boost::bind
        Boost::config
        Boost::core
        Boost::function
        Boost::iterator
        Boost::mpl
        Boost::optional
        Boost::parameter
        Boost::predef
        Boost::preprocessor
        Boost::smart_ptr
        Boost::throw_exception
        Boost::tuple
        Boost::type_traits
        Boost::variant
)

if( BOOST_SIGNALS2_INCLUDE_EXAMPLES )
    add_subdirectory( example )
endif()

