//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/awaitable.hpp>
#ifdef BOOST_ASIO_HAS_CO_AWAIT

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/pfr.hpp>
#include <boost/mysql/pool_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/static_results.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/experimental/cancellation_condition.hpp>
#include <boost/asio/experimental/parallel_group.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/redirect_error.hpp>
#include <boost/asio/this_coro.hpp>
#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <string>

#include "test_common/ci_server.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/run_coro.hpp"
#include "test_integration/snippets/credentials.hpp"

namespace mysql = boost::mysql;
namespace asio = boost::asio;
using namespace boost::mysql::test;

// Defined outside the namespace to prevent unused warnings
#if !defined(BOOST_ASIO_USE_TS_EXECUTOR_AS_DEFAULT)
asio::awaitable<void> dont_run(mysql::any_connection& conn)
{
    //[overview_async_dont
    // Coroutine body
    // DO NOT DO THIS!!!!
    mysql::results result1, result2;
    co_await asio::experimental::make_parallel_group(
        conn.async_execute("SELECT 1", result1, asio::deferred),
        conn.async_execute("SELECT 2", result2, asio::deferred)
    )
        .async_wait(asio::experimental::wait_for_all(), asio::deferred);
    //]
}
#endif

inline namespace overview {
//[overview_static_struct
// This must be placed at namespace scope.
// Should contain a member for each field of interest present in our query.
// Declaration order doesn't need to match field order in the query.
// Field names should match the ones in our query
struct employee
{
    std::int64_t id;
    std::string first_name;
    std::string last_name;
};
//]
}  // namespace overview

namespace {

asio::awaitable<void> overview_coro(mysql::any_connection& conn)
{
    std::string server_hostname = get_hostname();
    int employee_id = 1;
    const char* new_name = "John";

    auto& ctx = static_cast<asio::io_context&>((co_await asio::this_coro::executor).context());

    {
        //[overview_connect
        // The hostname, username, password and database to use.
        mysql::connect_params params;
        params.server_address.emplace_host_and_port(server_hostname);  // hostname
        params.username = mysql_username;
        params.password = mysql_password;
        params.database = "boost_mysql_examples";

        // Connect to the server
        co_await conn.async_connect(params);
        //]
    }
    {
        //[overview_text_query
        // Executes 'SELECT 1' and reads the resulting rows into memory
        mysql::results result;
        co_await conn.async_execute("SELECT 1", result);
        //]
    }

    {
        //[overview_with_params
        // If employee_id is 42, executes 'SELECT first_name FROM employee WHERE id = 42'
        mysql::results result;
        co_await conn.async_execute(
            mysql::with_params("SELECT first_name FROM employee WHERE id = {}", employee_id),
            result
        );
        //]
    }

    {
        //[overview_statement
        // First prepare the statement. Parsing happens server-side.
        mysql::statement stmt = co_await conn.async_prepare_statement(
            "SELECT first_name FROM employee WHERE company_id = ?"
        );

        // Now execute it. Parameter substitution happens server-side.
        mysql::results result;
        co_await conn.async_execute(stmt.bind(employee_id), result);
        //]
    }
    {
        //[overview_ifaces_dynamic
        // Passing a results to async_execute selects the dynamic interface
        mysql::results result;
        co_await conn.async_execute("SELECT id, first_name, last_name FROM employee", result);

        // Every employee is a collection of fields, which are variant-like objects
        // that represent data. We use as_string() to cast them to the appropriate type
        for (mysql::row_view emp : result.rows())
        {
            std::cout << "First name: " << emp.at(1).as_string() << ", last name: " << emp.at(2).as_string()
                      << std::endl;
        }
        //]
    }
#if BOOST_PFR_CORE_NAME_ENABLED
    {
        // The struct definition is included above this
        //[overview_ifaces_static
        //
        // This must be placed inside your function or method:
        //

        // Passing a static_results to async_execute selects the static interface
        mysql::static_results<mysql::pfr_by_name<employee>> result;
        co_await conn.async_execute("SELECT id, first_name, last_name FROM employee", result);

        // Query results are parsed directly into your own type
        for (const employee& emp : result.rows())
        {
            std::cout << "First name: " << emp.first_name << ", last name: " << emp.last_name << std::endl;
        }
        //]
    }
#endif
    {
        //[overview_update
        mysql::results result;
        co_await conn.async_execute(
            mysql::with_params("UPDATE employee SET first_name = {} WHERE id = {}", new_name, employee_id),
            result
        );
        //]
    }
    {
        //[overview_no_exceptions
        mysql::error_code ec;
        mysql::diagnostics diag;
        mysql::results result;

        // The provided SQL is invalid. The server will return an error.
        // ec will be set to a non-zero value, and diag will be populated
        co_await conn.async_execute("this is not SQL!", result, diag, asio::redirect_error(ec));
        if (ec)
        {
            // The error code will likely report a syntax error
            std::cout << "Operation failed with error code: " << ec << '\n';

            // diag.server_message() will contain the classic phrase
            // "You have an error in your SQL syntax; check the manual..."
            // Bear in mind that server_message() may contain user input, so treat it with caution
            std::cout << "Server diagnostics: " << diag.server_message() << std::endl;
        }
        //]
        BOOST_TEST(ec != mysql::error_code());
    }
    {
        //[overview_multifn
        // execution_state stores state about our operation, and must be passed to all functions
        mysql::execution_state st;

        // Writes the query request and reads the server response, but not the rows
        co_await conn.async_start_execution("SELECT first_name, last_name FROM employee", st);

        // Reads all the returned rows, in batches.
        // st.should_read_rows() returns false once there are no more rows to read
        while (st.should_read_rows())
        {
            // row_batch will be valid until conn performs the next network operation
            mysql::rows_view row_batch = co_await conn.async_read_some_rows(st);

            for (mysql::row_view emp : row_batch)
            {
                // Process the employee as required
                std::cout << "Name:" << emp.at(0) << " " << emp.at(1) << std::endl;
            }
        }
        //]
    }
    {
        //[overview_pool_create
        // pool_params contains configuration for the pool.
        // You must specify enough information to establish a connection,
        // including the server address and credentials.
        // You can configure a lot of other things, like pool limits
        mysql::pool_params params;
        params.server_address.emplace_host_and_port(server_hostname);
        params.username = mysql_username;
        params.password = mysql_password;
        params.database = "boost_mysql_examples";

        // Construct a pool of connections. The execution context will be used internally
        // to create the connections and other I/O objects
        mysql::connection_pool pool(ctx, std::move(params));

        // You need to call async_run on the pool before doing anything useful with it.
        // async_run creates connections and keeps them healthy. It must be called
        // only once per pool.
        // The detached completion token means that we don't want to be notified when
        // the operation ends. It's similar to a no-op callback.
        pool.async_run(asio::detached);
        //]

        // If we don't use the pool, we may leave unfinished work in the context
        co_await pool.async_get_connection();
    }
}

BOOST_FIXTURE_TEST_CASE(section_overview, any_connection_fixture)
{
    run_coro(ctx, [&]() { return overview_coro(conn); });
}

}  // namespace

#endif