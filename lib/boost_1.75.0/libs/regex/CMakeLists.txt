# Copyright 2018-2019 Mike Dev
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt
#
# NOTE: CMake support for Boost.Regex is currently experimental at best
#       and the interface is likely to change in the future


##### How-To:
#
# If you have a cmake project that wants to use and compile
# boost_regex, as part of a single build system run, do the following:
# 1) clone the boost project and all its sub-projects:
#
#   git clone --branch develop --depth 1 --recursive --shallow-submodules https://github.com/boostorg/boost.git boost-root
#
# 2) add to your cmake script:
#
#   add_subdirectory( <path-to-boost-root> [<build-dir-for-boost-libs>])
#   target_link_libraries( <my-exec> PUBLIC Boost::regex)
#
# 3) run your cmake build as usual
#
# ## Explanation:
#
# Currently this file does not work standalone. It is expected to be
# invoked from a parent script via add_subdirectory. That parent script
# is responsible for providing targets for direct and indirect dependencies,
# such as Boost::assert,  Boost::concept_check, e.g. by also adding those
# libraries via add_submodule (order doesn't matter).
# The parent script can be your own cmake script, but it is easier to just
# use add the CMakeLists in the root of the boost super project, which
# will in turn add all boost libraries usable with the add_subdirectory
# Workflow.
#
# Note: You don't need to actually clone all boost libraries. E.g. look
# into the travis ci file to see on which libraries boost_regex actually
# depends or use boostdep https://github.com/boostorg/boostdep


##### Current Limitations:
#
# - Doesn't compile or run tests
#

cmake_minimum_required( VERSION 3.5...3.16 )
project( boost_regex VERSION "${BOOST_SUPERPROJECT_VERSION}" LANGUAGES CXX )

option( BOOST_REGEX_INCLUDE_EXAMPLES "Also build (some) boost regex examples" OFF )
option( BOOST_REGEX_USE_ICU "Enable ICU support in boost regex" OFF )

file( GLOB BOOST_REGEX_SRC ./src/*.cpp )

add_library( boost_regex ${BOOST_REGEX_SRC} )
add_library( Boost::regex ALIAS boost_regex )

target_include_directories( boost_regex PUBLIC include )

target_compile_definitions( boost_regex
    PUBLIC
        # No need for autolink
        BOOST_REGEX_NO_LIB
        $<$<STREQUAL:$<TARGET_PROPERTY:boost_regex,TYPE>,SHARED_LIBRARY>:BOOST_REGEX_DYN_LINK=1>
        $<$<STREQUAL:$<TARGET_PROPERTY:boost_regex,TYPE>,STATIC_LIBRARY>:BOOST_REGEX_STATIC_LINK=1>
)

# Specify dependencies (including header-only libraries)
target_link_libraries( boost_regex
    PUBLIC
        Boost::assert
        Boost::concept_check
        Boost::config
        Boost::container_hash
        Boost::core
        Boost::integer
        Boost::iterator
        Boost::mpl
        Boost::predef
        Boost::smart_ptr
        Boost::static_assert
        Boost::throw_exception
        Boost::type_traits
)

if( BOOST_REGEX_USE_ICU )
    # ICU Targets could be provided by parent project,
    # if not, look for them ourselves
    if( NOT TARGET ICU::dt )
        # components need to be listed explicitly
        find_package( ICU COMPONENTS dt in uc REQUIRED )
    endif()

    target_link_libraries( boost_regex
        PRIVATE
            ICU::dt ICU::in ICU::uc
    )
    target_compile_definitions( boost_regex PRIVATE BOOST_HAS_ICU=1 )
endif()

if( BOOST_REGEX_INCLUDE_EXAMPLES )
    add_subdirectory( example/snippets )
endif()

