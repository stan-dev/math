//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/awaitable.hpp>
#include <boost/asio/local/stream_protocol.hpp>
#if defined(BOOST_ASIO_HAS_CO_AWAIT) && defined(BOOST_ASIO_HAS_LOCAL_SOCKETS)

//[example_unix_socket

/**
 * This example demonstrates how to connect to MySQL using a UNIX socket.
 *
 * It uses C++20 coroutines. If you need, you can backport
 * it to C++11 by using callbacks, asio::yield_context
 * or sync functions instead of coroutines.
 */

#include <boost/mysql/any_address.hpp>
#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/results.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>

#include <iostream>
#include <string_view>

namespace mysql = boost::mysql;
namespace asio = boost::asio;

// The main coroutine
asio::awaitable<void> coro_main(
    std::string_view unix_socket_path,
    std::string_view username,
    std::string_view password
)
{
    //[section_connection_establishment_unix_socket
    // Create a connection.
    // Will use the same executor as the coroutine.
    mysql::any_connection conn(co_await asio::this_coro::executor);

    // The socket path, username, password and database to use.
    // server_address is a variant-like type. Using emplace_unix_path,
    // we can specify a UNIX socket path, instead of a hostname and a port.
    // UNIX socket connections never use TLS.
    mysql::connect_params params;
    params.server_address.emplace_unix_path(std::string(unix_socket_path));
    params.username = username;
    params.password = password;
    params.database = "boost_mysql_examples";

    // Connect to the server
    co_await conn.async_connect(params);
    //]

    // The connection can now be used normally
    mysql::results result;
    co_await conn.async_execute("SELECT 'Hello world!'", result);
    std::cout << result.rows().at(0).at(0) << std::endl;

    // Notify the MySQL server we want to quit, then close the underlying connection.
    co_await conn.async_close();
}

void main_impl(int argc, char** argv)
{
    if (argc != 3 && argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> [<socket-path>]\n";
        exit(1);
    }

    // If not provided, use the default UNIX socket path,
    // compatible with most UNIX systems.
    const char* socket_path = argc >= 4 ? argv[3] : "/var/run/mysqld/mysqld.sock";

    // Create an I/O context, required by all I/O objects
    asio::io_context ctx;

    // Launch our coroutine
    asio::co_spawn(
        ctx,
        [=] { return coro_main(socket_path, argv[1], argv[2]); },
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
    std::cout << "Sorry, your compiler/system doesn't have the required capabilities to run this example"
              << std::endl;
}

#endif
