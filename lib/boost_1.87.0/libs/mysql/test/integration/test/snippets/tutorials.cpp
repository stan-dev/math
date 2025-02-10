//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/pfr.hpp>

#include <boost/asio/awaitable.hpp>
#if defined(BOOST_ASIO_HAS_CO_AWAIT) && BOOST_PFR_CORE_NAME_ENABLED

#include <boost/mysql/any_address.hpp>
#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/pfr.hpp>
#include <boost/mysql/pool_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/resultset_view.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/static_results.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/asio/as_tuple.hpp>
#include <boost/asio/cancel_after.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/read.hpp>
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/this_coro.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/system/error_code.hpp>
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <string>
#include <tuple>

#include "test_common/ci_server.hpp"
#include "test_integration/run_coro.hpp"
#include "test_integration/snippets/credentials.hpp"
#include "test_integration/snippets/snippets_fixture.hpp"

namespace mysql = boost::mysql;
namespace asio = boost::asio;
using namespace mysql::test;

// Common
inline namespace tutorials {
struct employee
{
    std::string first_name;
    std::string last_name;
};
}  // namespace tutorials

namespace {

//
// Tutorial 4: static interface
//
BOOST_FIXTURE_TEST_CASE(section_tutorial_static_interface, snippets_fixture)
{
    mysql::results result;
    conn.execute("SELECT first_name, last_name FROM employee WHERE id = 1", result);
    auto print_employee = [](mysql::string_view, mysql::string_view) {};

    //[tutorial_static_casts
    mysql::row_view employee = result.rows().at(0);
    print_employee(employee.at(0).as_string(), employee.at(1).as_string());
    //]
}

//
// Tutorial 5: updates and txns
//
asio::awaitable<void> tutorial_updates_transactions(mysql::any_connection& conn)
{
    const char* new_first_name = "John";
    int employee_id = 1;

    {
        //[tutorial_updates_transactions_update
        // Run an UPDATE. We can use with_params to compose it, too
        // If new_first_name contains 'John' and employee_id contains 42, this will run:
        //    UPDATE employee SET first_name = 'John' WHERE id = 42
        // result contains an empty resultset: it has no rows
        mysql::results result;
        co_await conn.async_execute(
            mysql::with_params(
                "UPDATE employee SET first_name = {} WHERE id = {}",
                new_first_name,
                employee_id
            ),
            result
        );
        //]
        //[tutorial_updates_transactions_select
        // Retrieve the newly created employee.
        // As we will see, this is a potential race condition
        // that can be avoided with transactions.
        co_await conn.async_execute(
            mysql::with_params("SELECT first_name, last_name FROM employee WHERE id = {}", employee_id),
            result
        );

        if (result.rows().empty())
        {
            std::cout << "No employee with ID = " << employee_id << std::endl;
        }
        else
        {
            std::cout << "Updated: " << result.rows().at(0).at(0) << " " << result.rows().at(0).at(1)
                      << std::endl;
        }
        //]
    }
    {
        //[tutorial_updates_transactions_txn
        mysql::results empty_result, select_result;

        // Start a transaction block. Subsequent statements will belong
        // to the transaction block, until a COMMIT or ROLLBACK is encountered,
        // or the connection is closed.
        // START TRANSACTION returns no rows.
        co_await conn.async_execute("START TRANSACTION", empty_result);

        // Run the UPDATE as we did before
        co_await conn.async_execute(
            mysql::with_params(
                "UPDATE employee SET first_name = {} WHERE id = {}",
                new_first_name,
                employee_id
            ),
            empty_result
        );

        // Run the SELECT. If a row is returned here, it is the one
        // that we modified.
        co_await conn.async_execute(
            mysql::with_params("SELECT first_name, last_name FROM employee WHERE id = {}", employee_id),
            select_result
        );

        // Commit the transaction. This makes the updated row visible
        // to other transactions and releases any locked rows.
        co_await conn.async_execute("COMMIT", empty_result);

        // Process the retrieved rows
        if (select_result.rows().empty())
        {
            std::cout << "No employee with ID = " << employee_id << std::endl;
        }
        else
        {
            std::cout << "Updated: " << select_result.rows().at(0).at(0) << " "
                      << select_result.rows().at(0).at(1) << std::endl;
        }
        //]
    }
    {
        //[tutorial_updates_transactions_multi_queries
        // Run the 4 statements in a single round-trip.
        // If an error is encountered, successive statements won't be executed
        // and the transaction won't be committed.
        mysql::results result;
        co_await conn.async_execute(
            mysql::with_params(
                "START TRANSACTION;"
                "UPDATE employee SET first_name = {} WHERE id = {};"
                "SELECT first_name, last_name FROM employee WHERE id = {};"
                "COMMIT",
                new_first_name,
                employee_id,
                employee_id
            ),
            result
        );
        //]

        //[tutorial_updates_transactions_dynamic_results
        // Get the 3rd resultset. resultset_view API is similar to results
        mysql::resultset_view select_result = result.at(2);
        if (select_result.rows().empty())
        {
            std::cout << "No employee with ID = " << employee_id << std::endl;
        }
        else
        {
            std::cout << "Updated: " << select_result.rows().at(0).at(0) << " "
                      << select_result.rows().at(0).at(1) << std::endl;
        }
        //]
    }
    {
        //[tutorial_updates_transactions_manual_indices
        // {0} will be replaced by the first format arg, {1} by the second
        mysql::results result;
        co_await conn.async_execute(
            mysql::with_params(
                "START TRANSACTION;"
                "UPDATE employee SET first_name = {0} WHERE id = {1};"
                "SELECT first_name, last_name FROM employee WHERE id = {1};"
                "COMMIT",
                new_first_name,
                employee_id
            ),
            result
        );
        //]
    }
}

BOOST_FIXTURE_TEST_CASE(section_tutorial_updates_transactions, snippets_fixture)
{
    run_coro(ctx, [&]() { return tutorial_updates_transactions(conn); });
}

//
// Tutorial 6: connection pool
//
mysql::pool_params create_pool_params()
{
    mysql::pool_params res;
    res.server_address.emplace_host_and_port(get_hostname());
    res.username = mysql_username;
    res.password = mysql_password;
    res.database = "boost_mysql_examples";
    return res;
}

BOOST_FIXTURE_TEST_CASE(section_tutorial_connection_pool, snippets_fixture)
{
    run_coro(ctx, [&]() -> asio::awaitable<void> {
        mysql::connection_pool pool(ctx, create_pool_params());
        pool.async_run(asio::detached);

        //[tutorial_connection_pool_get_connection
        // Get a connection from the pool.
        // This will wait until a healthy connection is ready to be used.
        // pooled_connection grants us exclusive access to the connection until
        // the object is destroyed
        mysql::pooled_connection conn = co_await pool.async_get_connection();
        //]
    });
}

//
// Tutorial 7: error handling
//
void log_error(const char*, boost::system::error_code) {}

//[tutorial_error_handling_db_nodiag
asio::awaitable<std::string> get_employee_details(mysql::connection_pool& pool, std::int64_t employee_id)
{
    // Get a connection from the pool.
    // This will wait until a healthy connection is ready to be used.
    // ec is an error code, conn is a pooled_connection
    auto [ec, conn] = co_await pool.async_get_connection(asio::as_tuple);
    if (ec)
    {
        // A connection couldn't be obtained.
        // This may be because a timeout happened.
        log_error("Error in async_get_connection", ec);
        co_return "ERROR";
    }

    // Use the connection normally to query the database.
    mysql::static_results<mysql::pfr_by_name<employee>> result;
    auto [ec2] = co_await conn->async_execute(
        mysql::with_params("SELECT first_name, last_name FROM employee WHERE id = {}", employee_id),
        result,
        asio::as_tuple
    );
    if (ec2)
    {
        log_error("Error running query", ec);
        co_return "ERROR";
    }

    // Handle the result as we did in the previous tutorial
    //<-
    co_return "";
    //->
}
//]

[[maybe_unused]]
//[tutorial_error_handling_session_as_tuple
asio::awaitable<void> handle_session(mysql::connection_pool& pool, asio::ip::tcp::socket client_socket)
{
    // Read the request from the client.
    unsigned char message[8]{};
    auto [ec1, bytes_read] = co_await asio::async_read(client_socket, asio::buffer(message), asio::as_tuple);
    if (ec1)
    {
        log_error("Error reading from the socket", ec1);
        co_return;
    }

    // Process the request as before (omitted)
    //<-
    boost::ignore_unused(pool);
    std::string response;
    //->

    // Write the response back to the client.
    auto [ec2, bytes_written] = co_await asio::async_write(
        client_socket,
        asio::buffer(response),
        asio::as_tuple
    );
    if (ec2)
    {
        log_error("Error writing to the socket", ec2);
        co_return;
    }
}
//]

asio::awaitable<void> tutorial_error_handling()
{
    // Setup
    mysql::connection_pool pool(co_await asio::this_coro::executor, create_pool_params());
    asio::steady_timer cv(co_await asio::this_coro::executor, (std::chrono::steady_clock::time_point::max)());
    pool.async_run([&](mysql::error_code) { cv.cancel(); });

    {
        //[tutorial_error_handling_callbacks
        // Function to call when async_get_connection completes
        auto on_available_connection = [](boost::system::error_code ec, mysql::pooled_connection conn) {
            // Do something useful with the connection
            //<-
            BOOST_TEST(ec == mysql::error_code());
            BOOST_TEST(conn.valid());
            //->
        };

        // Start the operation. on_available_connection will be called when the operation
        // completes. on_available_connection is the completion token.
        // When a callback is passed, async_get_connection returns void,
        // so we can't use co_await with it.
        pool.async_get_connection(on_available_connection);
        //]
    }

    {
        //[tutorial_error_handling_default_tokens
        // These two lines are equivalent.
        // Both of them can be read as "I want to use C++20 coroutines as my completion style"
        auto conn1 = co_await pool.async_get_connection();
        auto conn2 = co_await pool.async_get_connection(mysql::with_diagnostics(asio::deferred));
        //]

        BOOST_TEST(conn1.valid());
        BOOST_TEST(conn2.valid());
    }
    {
        //[tutorial_error_handling_adapter_tokens
        // Enable the use of the "s" suffix for std::chrono::seconds
        using namespace std::chrono_literals;

        // The following two lines are equivalent.
        // Both get a connection, waiting no more than 20s before cancelling the operation.
        // If no token is passed to cancel_after, the default one will be used,
        // which transforms the operation into an awaitable.
        // asio::cancel_after(20s) is usually termed "partial completion token"
        auto conn1 = co_await pool.async_get_connection(asio::cancel_after(20s));
        auto conn2 = co_await pool.async_get_connection(
            asio::cancel_after(20s, mysql::with_diagnostics(asio::deferred))
        );
        //]

        BOOST_TEST(conn1.valid());
        BOOST_TEST(conn2.valid());
    }

    {
        //[tutorial_error_handling_as_tuple
        // Passing asio::as_tuple transforms the operation's handler signature:
        //    Original:    void(error_code, mysql::pooled_connection)
        //    Transformed: void(std::tuple<error_code, mysql::pooled_connection>)
        // The transformed signature no longer has an error_code as first parameter,
        // so no automatic error code to exception transformation happens.
        std::tuple<boost::system::error_code, mysql::pooled_connection>
            res = co_await pool.async_get_connection(asio::as_tuple);
        //]

        BOOST_TEST(std::get<0>(res) == mysql::error_code());
    }

    {
        //[tutorial_error_handling_as_tuple_structured_bindings
        // ec is an error_code, conn is the mysql::pooled_connection.
        // If the operation fails, ec will be non-empty.
        auto [ec, conn] = co_await pool.async_get_connection(asio::as_tuple);
        //]

        BOOST_TEST(ec == mysql::error_code());
        BOOST_TEST(conn.valid());
    }

    {
        //[tutorial_error_handling_as_tuple_default_tokens
        // The following two lines are equivalent.
        // Both of them produce an awaitable that produces a tuple when awaited.
        auto [ec1, conn1] = co_await pool.async_get_connection(asio::as_tuple);
        auto [ec2, conn2] = co_await pool.async_get_connection(
            asio::as_tuple(mysql::with_diagnostics(asio::deferred))
        );
        //]

        BOOST_TEST(ec1 == mysql::error_code());
        BOOST_TEST(ec2 == mysql::error_code());
        BOOST_TEST(conn1.valid());
        BOOST_TEST(conn2.valid());
    }

    {
        using namespace std::chrono_literals;

        //[tutorial_error_handling_as_tuple_cancel_after
        // ec is an error_code, conn is the mysql::pooled_connection
        // Apply a timeout and don't throw on error
        auto [ec, conn] = co_await pool.async_get_connection(asio::cancel_after(20s, asio::as_tuple));
        //]

        BOOST_TEST(ec == mysql::error_code());
        BOOST_TEST(conn.valid());
    }

    // Call the functions requiring a pool
    co_await get_employee_details(pool, 1);

    // Cancel the pool and wait run to return, so no work is left in the io_context
    pool.cancel();
    boost::ignore_unused(co_await cv.async_wait(asio::as_tuple));
}

BOOST_FIXTURE_TEST_CASE(section_tutorial_error_handling, io_context_fixture)
{
    run_coro(ctx, &tutorial_error_handling);
}

}  // namespace

#endif
