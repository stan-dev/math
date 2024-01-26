#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
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

    if(MSVC)
        target_compile_definitions(
            ${TARGET_NAME}
            PUBLIC
            _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING # Warnings in C++17 for Asio
        )
        target_compile_options(${TARGET_NAME} PUBLIC /bigobj) # Prevent failures on Windows
    else()
        target_compile_options(${TARGET_NAME} PUBLIC -Wall -Wextra -pedantic -Werror)
    endif()

    set_target_properties(${TARGET_NAME} PROPERTIES CXX_EXTENSIONS OFF) # disable extensions

    # Valgrind
    if(BOOST_MYSQL_VALGRIND_TESTS)
        target_include_directories(${TARGET_NAME} PUBLIC ${VALGRIND_INCLUDE_DIR})
        target_compile_definitions(${TARGET_NAME} PUBLIC BOOST_MYSQL_VALGRIND_TESTS)
    endif()

    # Coverage
    if(BOOST_MYSQL_COVERAGE)
        target_compile_options(${TARGET_NAME} PUBLIC --coverage)
        target_link_options(${TARGET_NAME} PUBLIC --coverage)
    endif()
endfunction()

# Valgrind stuff
if(BOOST_MYSQL_VALGRIND_TESTS)
    # Locate executable
    find_program(VALGRIND_EXECUTABLE valgrind)

    if(NOT VALGRIND_EXECUTABLE)
        message(FATAL_ERROR "Cannot locate valgrind executable")
    endif()

    # Locate includes
    find_path(VALGRIND_INCLUDE_DIR "valgrind/memcheck.h")

    if(NOT VALGRIND_INCLUDE_DIR)
        message(FATAL_ERROR "Cannot locate valgrind include files")
    endif()

    # Path to suppressions. Don't move inside any function
    set(_SUPPRESSIONS_FILE "${CMAKE_CURRENT_LIST_DIR}/../tools/valgrind_suppressions.txt")

    # Helper to define tests
    function(add_memcheck_test)
        set(options "")
        set(oneValueArgs NAME TARGET)
        set(multiValueArgs ARGUMENTS)
        cmake_parse_arguments(
            AddMemcheckTest
            "${options}"
            "${oneValueArgs}"
            "${multiValueArgs}"
            ${ARGN}
        )

        add_test(
            NAME ${AddMemcheckTest_NAME}
            COMMAND
            ${VALGRIND_EXECUTABLE}
            --leak-check=full
            --error-limit=yes
            --suppressions=${_SUPPRESSIONS_FILE}
            --error-exitcode=1
            --gen-suppressions=all
            $<TARGET_FILE:${AddMemcheckTest_TARGET}>
            ${AddMemcheckTest_ARGUMENTS}
        )
    endfunction()
endif()
