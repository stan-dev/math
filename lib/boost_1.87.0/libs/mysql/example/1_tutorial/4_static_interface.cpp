//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/pfr.hpp>

#include <boost/asio/awaitable.hpp>
#if defined(BOOST_ASIO_HAS_CO_AWAIT) && BOOST_PFR_CORE_NAME_ENABLED

//[example_tutorial_static_interface

/**
 * This example shows how to use the static interface to parse
 * the results of a query into a C++ struct.
 * Like the previous tutorial, given an employee ID,
 * it prints their full name.
 *
 * It uses Boost.Pfr for reflection, which requires C++20.
 * You can backport it to C++14 if you need by using Boost.Describe.
 *
 * This example uses the 'boost_mysql_examples' database, which you
 * can get by running db_setup.sql.
 */

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/pfr.hpp>
#include <boost/mysql/static_results.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/this_coro.hpp>

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <string_view>

namespace mysql = boost::mysql;
namespace asio = boost::asio;

//[tutorial_static_fn
void print_employee(std::string_view first_name, std::string_view last_name)
{
    std::cout << "Employee's name is: " << first_name << ' ' << last_name << std::endl;
}
//]

//[tutorial_static_struct
// Should contain a member for each field of interest present in our query.
// Declaration order doesn't need to match field order in the query.
// Field names should match the ones in our query
struct employee
{
    std::string first_name;
    std::string last_name;
};
//]

asio::awaitable<void> coro_main(
    std::string_view server_hostname,
    std::string_view username,
    std::string_view password,
    std::int64_t employee_id
)
{
    // Represents a connection to the MySQL server.
    // The connection will use the same executor as the coroutine
    mysql::any_connection conn(co_await asio::this_coro::executor);

    // The hostname, username, password and database to use.
    mysql::connect_params params;
    params.server_address.emplace_host_and_port(std::string(server_hostname));
    params.username = username;
    params.password = password;
    params.database = "boost_mysql_examples";

    // Connect to the server
    co_await conn.async_connect(params);

    //[tutorial_static_execute
    // Using static_results will parse the result of our query
    // into instances of the employee type. Fields will be matched
    // by name, instead of by position.
    // pfr_by_name tells the library to use Boost.Pfr for reflection,
    // and to match fields by name.
    mysql::static_results<mysql::pfr_by_name<employee>> result;

    // Execute the query with the given parameters, performing the required
    // escaping to prevent SQL injection.
    co_await conn.async_execute(
        mysql::with_params("SELECT first_name, last_name FROM employee WHERE id = {}", employee_id),
        result
    );
    //]

    //[tutorial_static_results
    // Did we find an employee with that ID?
    if (result.rows().empty())
    {
        std::cout << "Employee not found" << std::endl;
    }
    else
    {
        // Print the retrieved details
        const employee& emp = result.rows()[0];
        print_employee(emp.first_name, emp.last_name);
    }
    //]

    // Close the connection
    co_await conn.async_close();
}

void main_impl(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> <employee-id>\n";
        exit(1);
    }

    // The execution context, required to run I/O operations.
    asio::io_context ctx;

    // Enqueue the coroutine for execution.
    asio::co_spawn(
        // The execution context where the coroutine will run
        ctx,

        // The coroutine to run. This must be a function taking no arguments
        // and returning an asio::awaitable<T>
        [argv] { return coro_main(argv[3], argv[1], argv[2], std::stoi(argv[4])); },

        // Callback to run when the coroutine completes.
        // If any exception is thrown in the coroutine body, propagate it to terminate the program.
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
    catch (const mysql::error_with_diagnostics& err)
    {
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

#else

#include <iostream>

int main()
{
    std::cout << "Sorry, your compiler doesn't have the required capabilities to run this example"
              << std::endl;
}

#endif