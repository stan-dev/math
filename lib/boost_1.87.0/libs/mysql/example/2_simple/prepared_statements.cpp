//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/awaitable.hpp>
#ifdef BOOST_ASIO_HAS_CO_AWAIT

//[example_prepared_statements

/**
 * This example demonstrates how to prepare, execute
 * and deallocate prepared statements. This program retrieves
 * all employees in a company, given its ID.
 *
 * It uses C++20 coroutines. If you need, you can backport
 * it to C++11 by using callbacks, asio::yield_context
 * or sync functions instead of coroutines.
 */

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/this_coro.hpp>

#include <iostream>
#include <string_view>

namespace mysql = boost::mysql;
namespace asio = boost::asio;

void print_employee(mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

// The main coroutine
asio::awaitable<void> coro_main(
    std::string_view server_hostname,
    std::string_view username,
    std::string_view password,
    std::string_view company_id
)
{
    // Create a connection. It will use the same executor as our coroutine
    mysql::any_connection conn(co_await asio::this_coro::executor);

    // The hostname, username, password and database to use
    mysql::connect_params params;
    params.server_address.emplace_host_and_port(std::string(server_hostname));
    params.username = username;
    params.password = password;
    params.database = "boost_mysql_examples";

    // Connect to the server
    co_await conn.async_connect(params);

    // Prepared statements can be used to execute queries with untrusted
    // parameters securely. They are an option to mysql::with_params,
    // but work server-side.
    // They are more complex but can yield more efficiency when retrieving
    // lots of numeric data, or when executing the same query several times with the same parameters.
    // Ask the server to prepare a statement and retrieve its handle
    mysql::statement stmt = co_await conn.async_prepare_statement(
        "SELECT first_name, last_name, salary FROM employee WHERE company_id = ?"
    );

    // Execute the statement. bind() must be passed as many parameters (number of ?)
    // as the statement has. bind() packages the statement handle with the parameters,
    // and async_execute sends them to the server
    mysql::results result;
    co_await conn.async_execute(stmt.bind(company_id), result);
    for (mysql::row_view employee : result.rows())
        print_employee(employee);

    // We can execute stmt as many times as we want, potentially with different
    // parameters, without the need to re-prepare it.

    // Once we're done with a statement, we can close it, to deallocate it from the server.
    // Closing the connection will also deallocate active statements, so this is not
    // strictly required here, but it's shown for completeness.
    // This can be relevant if you're using long-lived sessions.
    // Note that statement's destructor does NOT close the statement.
    co_await conn.async_close_statement(stmt);

    // Notify the MySQL server we want to quit, then close the underlying connection.
    co_await conn.async_close();
}

void main_impl(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> <company-id>\n";
        exit(1);
    }

    // Create an I/O context, required by all I/O objects
    asio::io_context ctx;

    // Launch our coroutine
    asio::co_spawn(
        ctx,
        [=] { return coro_main(argv[3], argv[1], argv[2], argv[4]); },
        // If any exception is thrown in the coroutine body, rethrow it.
        [](std::exception_ptr ptr) {
            if (ptr)
            {
                std::rethrow_exception(ptr);
            }
        }
    );

    // Calling run will actually execute the coroutine until completion
    ctx.run();

    std::cout << "Done\n";
}

int main(int argc, char** argv)
{
    try
    {
        main_impl(argc, argv);
    }
    catch (const boost::mysql::error_with_diagnostics& err)
    {
        // Some errors include additional diagnostics, like server-provided error messages.
        // Security note: diagnostics::server_message may contain user-supplied values (e.g. the
        // field value that caused the error) and is encoded using to the connection's character set
        // (UTF-8 by default). Treat is as untrusted input.
        std::cerr << "Error: " << err.what() << ", error code: " << err.code() << '\n'
                  << "Server diagnostics: " << err.get_diagnostics().server_message() << std::endl;
        return 1;
    }
    catch (const std::exception& err)
    {
        std::cerr << "Error: " << err.what() << std::endl;
        return 1;
    }
}

//]

#else

#include <iostream>

int main()
{
    std::cout << "Sorry, your compiler doesn't have the required capabilities to run this example"
              << std::endl;
}

#endif