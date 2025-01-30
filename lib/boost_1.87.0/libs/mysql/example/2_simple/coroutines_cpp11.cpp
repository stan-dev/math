//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_coroutines_cpp11

/**
 * This example demonstrates how to use stackful coroutines when using async functions.
 * This can be a good choice when targeting a standard lower than C++20.
 * This example uses the 'boost_mysql_examples' database, which you
 * can get by running db_setup.sql.
 * You need to link your program to Boost.Context to use asio::spawn.
 *
 * This example uses the 'boost_mysql_examples' database, which you
 * can get by running db_setup.sql.
 */

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/with_diagnostics.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/spawn.hpp>

#include <iostream>

namespace mysql = boost::mysql;
namespace asio = boost::asio;

void print_employee(mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

/**
 * The main coroutine. It will suspend every time we call one of the asynchronous functions, saving
 * all information it needs for resuming. When the asynchronous operation completes,
 * the coroutine will resume in the point it was left.
 * We need to pass the yield object to async functions for this to work.
 */
void coro_main(
    const char* server_hostname,
    const char* username,
    const char* password,
    const char* company_id,
    asio::yield_context yield
)
{
    // Represents a connection to the MySQL server.
    // The connection will use the same executor as the coroutine
    mysql::any_connection conn(yield.get_executor());

    // The hostname, username, password and database to use
    mysql::connect_params conn_params;
    conn_params.server_address.emplace_host_and_port(server_hostname);
    conn_params.username = username;
    conn_params.password = password;
    conn_params.database = "boost_mysql_examples";

    // Connect to server. with_diagnostics turns thrown exceptions
    // into error_with_diagnostics, which contain more info than regular exceptions
    conn.async_connect(conn_params, mysql::with_diagnostics(yield));

    // Initiate the query execution. company_id is an untrusted value.
    // with_params will securely compose a SQL query and send it to the server for execution.
    // Returned rows will be read into result.
    mysql::results result;
    conn.async_execute(
        mysql::with_params(
            "SELECT first_name, last_name, salary FROM employee WHERE company_id = {}",
            company_id
        ),
        result,
        mysql::with_diagnostics(yield)
    );

    // Print the employees
    for (boost::mysql::row_view employee : result.rows())
    {
        print_employee(employee);
    }

    // Notify the MySQL server we want to quit, then close the underlying connection.
    conn.async_close(mysql::with_diagnostics(yield));
}

void main_impl(int argc, char** argv)
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> [company-id]\n";
        exit(1);
    }

    // The company_id whose employees we will be listing. This
    // is user-supplied input, and should be treated as untrusted.
    const char* company_id = argc == 5 ? argv[4] : "HGS";

    // The execution context, required to run I/O operations.
    asio::io_context ctx;

    // Launch the coroutine
    asio::spawn(
        ctx,
        [argv, company_id](asio::yield_context yield) {
            coro_main(argv[3], argv[1], argv[2], company_id, yield);
        },
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
}

int main(int argc, char** argv)
{
    try
    {
        main_impl(argc, argv);
    }
    catch (const boost::mysql::error_with_diagnostics& err)
    {
        // You will only get this type of exceptions if you use with_diagnostics.
        // Some errors include additional diagnostics, like server-provided error messages.
        // Security note: diagnostics::server_message may contain user-supplied values (e.g. the
        // field value that caused the error) and is encoded using to the connection's character set
        // (UTF-8 by default). Treat is as untrusted input.
        std::cerr << "Error: " << err.what() << '\n'
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
