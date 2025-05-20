#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# Sets _WIN32_WINNT on Windows
function(boost_mysql_set_windows_version TARGET_NAME)
    if(MSVC)
        if(WIN32 AND CMAKE_SYSTEM_VERSION)
            set(WINNT_VERSION ${CMAKE_SYSTEM_VERSION})
            string(REPLACE "." "" WINNT_VERSION ${WINNT_VERSION})
            string(REGEX REPLACE "([0-9])" "0\\1" WINNT_VERSION ${WINNT_VERSION})

            set(WINNT_VERSION "0x${WINNT_VERSION}")
        else()
            set(WINNT_VERSION "0x0601")
        endif()

        target_compile_definitions(
            ${TARGET_NAME}
            PUBLIC
            _WIN32_WINNT=${WINNT_VERSION} # Silence warnings in Windows
        )
    endif()
endfunction()

# Utility function to set warnings and other compile properties of
# our test targets
function(boost_mysql_common_target_settings TARGET_NAME)
    boost_mysql_set_windows_version(${TARGET_NAME})

    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        target_compile_definitions(
            ${TARGET_NAME}
            PUBLIC
            _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING # Warnings in C++17 for Asio
        )
        target_compile_options(${TARGET_NAME} PUBLIC /bigobj) # Prevent failures on Windows
    else()
        # gcc-13 doesn't understand view types
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13.0)
            target_compile_options(${TARGET_NAME} PUBLIC -Wno-dangling-reference -Wno-array-bounds)
        endif()
        target_compile_options(${TARGET_NAME} PUBLIC -Wall -Wextra)
    endif()

    set_target_properties(${TARGET_NAME} PROPERTIES CXX_EXTENSIONS OFF) # disable extensions

    # Follow the Boost convention: don't build test targets by default,
    # and only when explicitly requested by building target tests
    set_target_properties(${TARGET_NAME} PROPERTIES EXCLUDE_FROM_ALL ON)
    add_dependencies(tests ${TARGET_NAME})
endfunction()
