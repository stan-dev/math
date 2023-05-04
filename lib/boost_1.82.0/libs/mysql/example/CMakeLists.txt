#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# Note: examples count as integration tests. This is only processed
# when we're the main project and BOOST_MYSQL_INTEGRATION_TESTS is on
find_package(Boost ${BOOST_MYSQL_VERSION} REQUIRED COMPONENTS context)

# Compile the example
function (build_example EXECUTABLE_NAME CPPFILE)
    add_executable(
        ${EXECUTABLE_NAME}
        ${CPPFILE}
    )
    target_link_libraries(
        ${EXECUTABLE_NAME}
        PRIVATE
        Boost::disable_autolinking
        Boost::mysql
        Boost::context
        boost_mysql_asio
    )
    common_target_settings(${EXECUTABLE_NAME})
endfunction()

# Run it as a test
function (test_example EXECUTABLE_NAME CONNECTION_ARG)
    add_test(
        NAME ${EXECUTABLE_NAME}
        COMMAND ${EXECUTABLE_NAME} example_user example_password ${CONNECTION_ARG}
    )
endfunction()

# Run it as a test using Valgrind
function (memcheck_example EXECUTABLE_NAME CONNECTION_ARG)
    set(TEST_NAME "${EXECUTABLE_NAME}_memcheck")
    add_memcheck_test(
        NAME ${TEST_NAME}
        TARGET ${EXECUTABLE_NAME}
        ARGUMENTS example_user example_password ${CONNECTION_ARG}
    )
endfunction()

# Build and run it as either of them
function (add_example EXAMPLE_NAME DO_MEMCHECK CONNECTION_ARG)
    set(EXECUTABLE_NAME "boost_mysql_example_${EXAMPLE_NAME}")

    build_example(${EXECUTABLE_NAME} "${EXAMPLE_NAME}.cpp")

    if (BOOST_MYSQL_VALGRIND_TESTS AND DO_MEMCHECK)
        memcheck_example(${EXECUTABLE_NAME} ${CONNECTION_ARG})
    else()
        test_example(${EXECUTABLE_NAME} ${CONNECTION_ARG})
    endif()
endfunction()

# Get the MySQL hostname to use for examples
if(DEFINED ENV{BOOST_MYSQL_SERVER_HOST})
    set(SERVER_HOST $ENV{BOOST_MYSQL_SERVER_HOST})
else()
    set(SERVER_HOST "localhost")
endif()

# Build and run all the examples
add_example(snippets                    TRUE    ${SERVER_HOST})
add_example(tutorial                    TRUE    ${SERVER_HOST})
add_example(text_queries                TRUE    ${SERVER_HOST})
add_example(prepared_statements         TRUE    ${SERVER_HOST})
add_example(async_callbacks             TRUE    ${SERVER_HOST})
add_example(async_coroutines            FALSE   ${SERVER_HOST})
add_example(async_coroutinescpp20       TRUE    ${SERVER_HOST})
add_example(async_futures               TRUE    ${SERVER_HOST})
add_example(default_completion_tokens   TRUE    ${SERVER_HOST})
add_example(metadata                    TRUE    ${SERVER_HOST})
add_example(ssl                         TRUE    ${SERVER_HOST})
add_example(timeouts                    TRUE    ${SERVER_HOST})
if ("$ENV{BOOST_MYSQL_NO_UNIX_SOCKET_TESTS}" STREQUAL "")
    add_example(unix_socket             TRUE    "/var/run/mysqld/mysqld.sock")
endif()

