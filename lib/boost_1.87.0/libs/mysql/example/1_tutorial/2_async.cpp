//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/awaitable.hpp>
#ifdef BOOST_ASIO_HAS_CO_AWAIT

//[example_tutorial_async

/**
 * This example is analogous to the synchronous tutorial, but uses async functions
 * with C++20 coroutines, instead.
 */

#include <boost/mysql/any_address.hpp>
#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/results.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/this_coro.hpp>

#include <exception>
#include <iostream>

namespace mysql = boost::mysql;
namespace asio = boost::asio;

/**
 * The main coroutine.
 * It must have a return type of asio::awaitable<T>.
 * Our coroutine does not communicate any result back, so T=void.
 *
 * The coroutine will suspend every time we call one of the asynchronous functions, saving
 * all information it needs for resuming. When the asynchronous operation completes,
 * the coroutine will resume in the point it was left.
 * We use the same program structure as in the sync world, replacing
 * sync functions by their async equivalents and adding co_await in front of them.
 */
//[tutorial_async_coro
asio::awaitable<void> coro_main(
    mysql::any_connection& conn,
    std::string_view server_hostname,
    std::string_view username,
    std::string_view password
)
{
    // The hostname, username, password and database to use.
    // TLS is used by default.
    mysql::connect_params params;
    params.server_address.emplace_host_and_port(std::string(server_hostname));
    params.username = username;
    params.password = password;

    // Connect to the server
    co_await conn.async_connect(params);

    // Issue the SQL query to the server
    const char* sql = "SELECT 'Hello world!'";
    mysql::results result;
    co_await conn.async_execute(sql, result);

    // Print the first field in the first row
    std::cout << result.rows().at(0).at(0) << std::endl;

    // Close the connection
    co_await conn.async_close();
}
//]

void main_impl(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname>\n";
        exit(1);
    }

    //[tutorial_async_connection
    // The execution context, required to run I/O operations.
    asio::io_context ctx;

    // Represents a connection to the MySQL server.
    mysql::any_connection conn(ctx);
    //]

    //[tutorial_async_co_spawn
    // Enqueue the coroutine for execution.
    // This does not wait for the coroutine to finish.
    asio::co_spawn(
        // The execution context where the coroutine will run
        ctx,

        // The coroutine to run. This must be a function taking no arguments
        // and returning an asio::awaitable<T>
        [&conn, argv] { return coro_main(conn, argv[3], argv[1], argv[2]); },

        // Callback to run when the coroutine completes.
        // If any exception is thrown in the coroutine body, propagate it to terminate the program.
        [](std::exception_ptr ptr) {
            if (ptr)
            {
                std::rethrow_exception(ptr);
            }
        }
    );
    //]

    //[tutorial_async_run
    // Calling run will actually execute the coroutine until completion
    ctx.run();
    //]
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