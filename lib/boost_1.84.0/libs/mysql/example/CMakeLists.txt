#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# Note: examples count as integration tests. This is only processed
# when we're the main project and BOOST_MYSQL_INTEGRATION_TESTS is on
find_package(Boost ${BOOST_MYSQL_VERSION} REQUIRED COMPONENTS context)

# Get the MySQL hostname to use for examples
if(DEFINED ENV{BOOST_MYSQL_SERVER_HOST})
    set(SERVER_HOST $ENV{BOOST_MYSQL_SERVER_HOST})
else()
    set(SERVER_HOST "localhost")
endif()

add_library(boost_mysql_examples_common INTERFACE)
target_link_libraries(
    boost_mysql_examples_common
    INTERFACE
    Boost::disable_autolinking
    boost_mysql_compiled
)

# Declare an example target
function (add_example_target EXAMPLE_NAME EXAMPLE_PATH)
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_PATH})
    target_link_libraries(${EXAMPLE_NAME} PRIVATE boost_mysql_examples_common ${ARGN})
    boost_mysql_common_target_settings(${EXAMPLE_NAME})
endfunction()

# Build and run an example
function (add_example EXAMPLE_NAME EXAMPLE_PATH)
    add_example_target(${EXAMPLE_NAME} ${EXAMPLE_PATH})

    if (BOOST_MYSQL_VALGRIND_TESTS)
        add_memcheck_test(
            NAME "${EXAMPLE_NAME}_memcheck"
            TARGET ${EXAMPLE_NAME}
            ARGUMENTS ${ARGN}
        )
    else()
        add_test(
            NAME ${EXAMPLE_NAME}
            COMMAND ${EXAMPLE_NAME} ${ARGN}
        )
    endif()
endfunction()

# Coroutines needs Boost.Context and shouldn't be memchecked
function (add_example_coroutines)
    set(EXECUTABLE_NAME boost_mysql_example_async_coroutines)
    add_example_target(${EXECUTABLE_NAME} async_coroutines.cpp Boost::context)
    add_test(
        NAME ${EXECUTABLE_NAME}
        COMMAND ${EXECUTABLE_NAME} example_user example_password ${SERVER_HOST}
    )
endfunction()

# The order management examples must be run several times through a Python script
function (add_example_order_management EXAMPLE_NAME EXAMPLE_PATH)
    add_example_target(${EXAMPLE_NAME} ${EXAMPLE_PATH})
    add_test(
        NAME ${EXAMPLE_NAME}
        COMMAND
            python
            ${CMAKE_CURRENT_SOURCE_DIR}/private/run_stored_procedures.py
            $<TARGET_FILE:${EXAMPLE_NAME}>
            ${SERVER_HOST}
    )
endfunction()

# Regular examples are the ones that require no extra linking libs and can be run
# with example_user example_password
set (REGULAR_EXAMPLES
    snippets
    tutorial
    async_callbacks
    async_coroutinescpp20
    async_futures
    metadata
    ssl
    timeouts
)
foreach (FILE_NAME ${REGULAR_EXAMPLES})
    add_example(
        "boost_mysql_example_${FILE_NAME}"
        "${FILE_NAME}.cpp"
        example_user example_password ${SERVER_HOST}
    )
endforeach()

# Order management examples must be run several times through a Python script
set(ORDER_EXAMPLES
    prepared_statements_cpp11
    prepared_statements_cpp14
    stored_procedures_cpp11
    stored_procedures_cpp14
)
foreach (FILE_NAME ${ORDER_EXAMPLES})
    add_example_order_management(
        "boost_mysql_example_${FILE_NAME}"
        "order_management/${FILE_NAME}.cpp"
    )
endforeach()

# Rest of the examples
add_example_coroutines()
if ("$ENV{BOOST_MYSQL_NO_UNIX_SOCKET_TESTS}" STREQUAL "")
    add_example(
        boost_mysql_example_unix_socket
        unix_socket.cpp
        example_user example_password
    )
endif()
add_example(
    boost_mysql_example_source_script
    source_script.cpp
    example_user example_password ${SERVER_HOST} ${CMAKE_CURRENT_SOURCE_DIR}/private/test_script.sql
)
