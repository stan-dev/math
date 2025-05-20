//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/awaitable.hpp>
#ifdef BOOST_ASIO_HAS_CO_AWAIT

//[example_metadata

/**
 * This example shows how to obtain metadata from SQL queries,
 * including field and table names.
 *
 * This example uses C++20 coroutines. If you need, you can backport
 * it to C++11 by using callbacks, asio::yield_context
 * or sync functions instead of coroutines.
 *
 * This example uses the 'boost_mysql_examples' database, which you
 * can get by running db_setup.sql.
 */

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>

#include <iostream>

namespace asio = boost::asio;
namespace mysql = boost::mysql;

// The main coroutine
asio::awaitable<void> coro_main(
    std::string_view server_hostname,
    std::string_view username,
    std::string_view password
)
{
    // Create a connection.
    // Will use the same executor as the coroutine.
    mysql::any_connection conn(co_await asio::this_coro::executor);

    // By default, string metadata (like column names) won't be retained.
    // This is for efficiency reasons. You can change this setting by calling
    // connection::set_meta_mode. It will affect any subsequent queries and statement executions.
    conn.set_meta_mode(mysql::metadata_mode::full);

    // The socket path, username, password and database to use.
    mysql::connect_params params;
    params.server_address.emplace_host_and_port(std::string(server_hostname));
    params.username = username;
    params.password = password;
    params.database = "boost_mysql_examples";

    // Connect to the server
    co_await conn.async_connect(params);

    // Issue the query
    constexpr const char* sql = R"(
        SELECT comp.name AS company_name, emp.id AS employee_id
        FROM employee emp
        JOIN company comp ON (comp.id = emp.company_id)
    )";
    mysql::results result;
    co_await conn.async_execute(sql, result);

    /**
     * results objects allow you to access metadata about the columns in the query
     * using the meta() function, which returns span-like object containing metadata objects
     * (one per column in the query, and in the same order as in the query).
     * You can retrieve the column name, type, number of decimals,
     * suggested display width, whether the column is part of a key...
     * These metadata objects are owned by the results object.
     */
    assert(result.meta().size() == 2);

    // <-
    // clang-format off
    // ->
    const mysql::metadata& company_name = result.meta()[0];
    assert(company_name.database() == "boost_mysql_examples");  // database name
    assert(company_name.table() == "comp");  // the alias we assigned to the table in the query
    assert(company_name.original_table() == "company");   // the original table name
    assert(company_name.column_name() == "company_name");  // the name of the column in the query
    assert(company_name.original_column_name() == "name");  // the name of the physical column in the table
    assert(company_name.type() == boost::mysql::column_type::varchar);  // we created the column as a VARCHAR
    assert(!company_name.is_primary_key());     // column is not a primary key
    assert(!company_name.is_auto_increment());  // column is not AUTO_INCREMENT
    assert(company_name.is_not_null());         // column may not be NULL

    const mysql::metadata& employee_id = result.meta()[1];
    assert(employee_id.database() == "boost_mysql_examples");  // database name
    assert(employee_id.table() == "emp");  // the alias we assigned to the table in the query
    assert(employee_id.original_table() == "employee");  // the original table name
    assert(employee_id.column_name() == "employee_id");   // the name of the column in the query
    assert(employee_id.original_column_name() == "id");  // the name of the physical column in the table
    assert(employee_id.type() == boost::mysql::column_type::int_);  // we created the column as INT
    assert(employee_id.is_primary_key()); // column is a primary key
    assert(employee_id.is_auto_increment()); // we declared the column as AUTO_INCREMENT
    assert(employee_id.is_not_null()); // column cannot be NULL
    // <-
    // clang-format on
    // avoid warnings in release mode
    static_cast<void>(company_name);
    static_cast<void>(employee_id);
    // ->

    // Notify the MySQL server we want to quit, then close the underlying connection.
    co_await conn.async_close();
}

void main_impl(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname>\n";
        exit(1);
    }

    // Create an I/O context, required by all I/O objects
    asio::io_context ctx;

    // Launch our coroutine
    asio::co_spawn(
        ctx,
        [=] { return coro_main(argv[3], argv[1], argv[2]); },
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
