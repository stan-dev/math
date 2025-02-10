#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# Note: examples count as integration tests. This is only processed
# when BOOST_MYSQL_INTEGRATION_TESTS is on

# Get the MySQL hostname to use for examples
if(DEFINED ENV{BOOST_MYSQL_SERVER_HOST})
    set(SERVER_HOST $ENV{BOOST_MYSQL_SERVER_HOST})
else()
    set(SERVER_HOST "127.0.0.1")
endif()

add_library(boost_mysql_examples_common INTERFACE)
target_link_libraries(
    boost_mysql_examples_common
    INTERFACE
    boost_mysql_compiled
)


function(add_example EXAMPLE_NAME)
    # Parse the arguments
    set(ONE_VALUE_ARGS PYTHON_RUNNER)
    set(MULTI_VALUE_ARGS SOURCES LIBS ARGS)
    cmake_parse_arguments(ADD_EXAMPLE "" "${ONE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}" ${ARGN})

    # Create the target
    set(TARGET_NAME "boost_mysql_example_${EXAMPLE_NAME}")
    add_executable(${TARGET_NAME} ${ADD_EXAMPLE_SOURCES})
    target_link_libraries(${TARGET_NAME} PRIVATE boost_mysql_examples_common)
    boost_mysql_common_target_settings(${TARGET_NAME})
    target_link_libraries(${TARGET_NAME} PRIVATE ${ADD_EXAMPLE_LIBS})

    # Add it as a test
    if (ADD_EXAMPLE_PYTHON_RUNNER)
        add_test(
            NAME ${TARGET_NAME}
            COMMAND
            python
            ${CMAKE_CURRENT_SOURCE_DIR}/private/${ADD_EXAMPLE_PYTHON_RUNNER}
            $<TARGET_FILE:${TARGET_NAME}>
            ${ADD_EXAMPLE_ARGS}
        )
    else()
        add_test(
            NAME ${TARGET_NAME}
            COMMAND ${TARGET_NAME} ${ADD_EXAMPLE_ARGS}
        )
    endif()
endfunction()

function(add_tutorial EXAMPLE_NAME EXAMPLE_PATH)
    add_example(${EXAMPLE_NAME} SOURCES "1_tutorial/${EXAMPLE_PATH}" ${ARGN})
endfunction()

function(add_simple_example EXAMPLE_NAME)
    add_example(${EXAMPLE_NAME} SOURCES "2_simple/${EXAMPLE_NAME}.cpp" ${ARGN})
endfunction()

set(REGULAR_ARGS example_user example_password ${SERVER_HOST})

# Tutorials
add_tutorial(tutorial_sync                 1_sync.cpp                  ARGS ${REGULAR_ARGS})
add_tutorial(tutorial_async                2_async.cpp                 ARGS ${REGULAR_ARGS})
add_tutorial(tutorial_with_params          3_with_params.cpp           ARGS ${REGULAR_ARGS} 1)
add_tutorial(tutorial_static_interface     4_static_interface.cpp      ARGS ${REGULAR_ARGS} 1        LIBS Boost::pfr)
add_tutorial(tutorial_updates_transactions 5_updates_transactions.cpp  ARGS ${REGULAR_ARGS} 1 "John" LIBS Boost::pfr)
add_tutorial(tutorial_connection_pool      6_connection_pool.cpp       ARGS ${SERVER_HOST}           LIBS Boost::pfr
    PYTHON_RUNNER run_tutorial_connection_pool.py)
add_tutorial(tutorial_error_handling       7_error_handling.cpp        ARGS ${SERVER_HOST} --test-errors  LIBS Boost::pfr
    PYTHON_RUNNER run_tutorial_connection_pool.py)

# Simple
add_simple_example(inserts                      ARGS ${REGULAR_ARGS} "John" "Doe" "HGS")
add_simple_example(deletes                      ARGS ${REGULAR_ARGS} 20)
add_simple_example(callbacks                    ARGS ${REGULAR_ARGS})
add_simple_example(coroutines_cpp11             ARGS ${REGULAR_ARGS} LIBS Boost::context)
add_simple_example(batch_inserts                ARGS ${SERVER_HOST}  PYTHON_RUNNER run_batch_inserts.py LIBS Boost::json)
add_simple_example(batch_inserts_generic        ARGS ${SERVER_HOST}  PYTHON_RUNNER run_batch_inserts.py LIBS Boost::json)
add_simple_example(patch_updates                ARGS ${SERVER_HOST}  PYTHON_RUNNER run_patch_updates.py)
add_simple_example(dynamic_filters              ARGS ${SERVER_HOST}  PYTHON_RUNNER run_dynamic_filters.py)
add_simple_example(disable_tls                  ARGS ${REGULAR_ARGS})
add_simple_example(tls_certificate_verification ARGS ${REGULAR_ARGS})
add_simple_example(metadata                     ARGS ${REGULAR_ARGS})
add_simple_example(prepared_statements          ARGS ${REGULAR_ARGS} "HGS")
add_simple_example(pipeline                     ARGS ${REGULAR_ARGS} "HGS")
add_simple_example(multi_function               ARGS ${REGULAR_ARGS})
add_simple_example(source_script                ARGS ${REGULAR_ARGS} ${CMAKE_CURRENT_SOURCE_DIR}/private/test_script.sql)

# UNIX sockets. Don't run the example on Windows machines
if (NOT WIN32)
    add_simple_example(unix_socket ARGS example_user example_password)
endif()

# Advanced
add_example(
    connection_pool
    SOURCES
        3_advanced/connection_pool/repository.cpp
        3_advanced/connection_pool/handle_request.cpp
        3_advanced/connection_pool/server.cpp
        3_advanced/connection_pool/main.cpp
    LIBS
        Boost::context
        Boost::json
        Boost::url
        Boost::beast
    PYTHON_RUNNER run_connection_pool.py
    ARGS ${SERVER_HOST}
)
