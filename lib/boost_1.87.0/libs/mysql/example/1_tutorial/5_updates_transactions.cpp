//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/pfr.hpp>

#include <boost/asio/awaitable.hpp>
#if defined(BOOST_ASIO_HAS_CO_AWAIT) && BOOST_PFR_CORE_NAME_ENABLED

//[example_tutorial_updates_transactions

/**
 * This example demonstrates how to use UPDATE statements,
 * transactions and semicolon-separated queries.
 *
 * The program updates the first name of an employee given their ID
 * and prints their full details.
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
#include <boost/mysql/results.hpp>
#include <boost/mysql/resultset_view.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/static_results.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>

#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <tuple>

namespace mysql = boost::mysql;
namespace asio = boost::asio;

// As in the previous tutorial, this struct models
// the data returned by our SELECT query. It should contain a member
// for each field of interest, with a matching name.
struct employee
{
    std::string first_name;
    std::string last_name;
};

// The main coroutine
asio::awaitable<void> coro_main(
    std::string_view server_hostname,
    std::string_view username,
    std::string_view password,
    std::int64_t employee_id,
    std::string_view new_first_name
)
{
    // Create a connection.
    // Will use the same executor as the coroutine.
    mysql::any_connection conn(co_await asio::this_coro::executor);

    //[section_connection_establishment_multi_queries
    //[tutorial_updates_transactions_connect
    // The server host, username, password and database to use.
    // Setting multi_queries to true makes it possible to run several
    // semicolon-separated queries with async_execute.
    mysql::connect_params params;
    params.server_address.emplace_host_and_port(std::string(server_hostname));
    params.username = std::move(username);
    params.password = std::move(password);
    params.database = "boost_mysql_examples";
    params.multi_queries = true;

    // Connect to the server
    co_await conn.async_connect(params);
    //]
    //]

    // Perform the update and retrieve the results:
    //   1. Begin a transaction block. Further updates won't be visible to
    //      other transactions until this one commits.
    //   2. Perform the update.
    //   3. Retrieve the employee we just updated. Since we're in a transaction,
    //      this will be the employee we just updated (if any),
    //      without the possibility of other transactions interfering.
    //   4. Commit the transaction and make everything visible to other transactions.
    //      If any of the previous steps fail, the commit won't be run, and the
    //      transaction will be rolled back when the connection is closed.
    //[tutorial_updates_transactions_static
    // MySQL returns one resultset for each query, so we pass 4 params to static_results
    //<-
    // clang-format off
    //->
    mysql::static_results<
        std::tuple<>,                  // START TRANSACTION doesn't generate rows
        std::tuple<>,                  // The UPDATE doesn't generate rows
        mysql::pfr_by_name<employee>,  // The SELECT generates employees
        std::tuple<>                   // The COMMIT doesn't generate rows
    > result;
    //<-
    // clang-format on
    //->

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

    // We've run 4 SQL queries, so MySQL has returned us 4 resultsets.
    // The SELECT is the 3rd resultset. Retrieve the generated rows.
    // employees is a span<const employee>
    auto employees = result.rows<2>();
    if (employees.empty())
    {
        std::cout << "No employee with ID = " << employee_id << std::endl;
    }
    else
    {
        const employee& emp = employees[0];
        std::cout << "Updated: employee is now " << emp.first_name << " " << emp.last_name << std::endl;
    }
    //]

    // Notify the MySQL server we want to quit, then close the underlying connection.
    co_await conn.async_close();
}

void main_impl(int argc, char** argv)
{
    if (argc != 6)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <username> <password> <server-hostname> <employee-id> <new-first-name>\n";
        exit(1);
    }

    // Create an I/O context, required by all I/O objects
    asio::io_context ctx;

    // Launch our coroutine
    asio::co_spawn(
        ctx,
        [=] { return coro_main(argv[3], argv[1], argv[2], std::stoi(argv[4]), argv[5]); },
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
