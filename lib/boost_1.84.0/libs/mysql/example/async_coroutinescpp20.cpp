//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_async_coroutinescpp20

#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/tcp_ssl.hpp>
#include <boost/mysql/throw_on_error.hpp>

#include <boost/asio/as_tuple.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/use_awaitable.hpp>

#include <exception>
#include <iostream>

using boost::mysql::error_code;

#ifdef BOOST_ASIO_HAS_CO_AWAIT

void print_employee(boost::mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

// Using this completion token instead of plain use_awaitable prevents
// co_await from throwing exceptions. Instead, co_await will return a std::tuple<error_code>
// with a non-zero code on error. We will then use boost::mysql::throw_on_error
// to throw exceptions with embedded diagnostics, if available. If you
// employ plain use_awaitable, you will get boost::system::system_error exceptions
// instead of boost::mysql::error_with_diagnostics exceptions. This is a limitation of use_awaitable.
constexpr auto tuple_awaitable = boost::asio::as_tuple(boost::asio::use_awaitable);

/**
 * Our coroutine. It must have a return type of boost::asio::awaitable<T>.
 * Our coroutine does not communicate any result back, so T=void.
 * Remember that you do not have to explicitly create any awaitable<void> in
 * your function. Instead, the return type is fed to std::coroutine_traits
 * to determine the semantics of the coroutine, like the promise type.
 * Asio already takes care of all this for us.
 *
 * The coroutine will suspend every time we call one of the asynchronous functions, saving
 * all information it needs for resuming. When the asynchronous operation completes,
 * the coroutine will resume in the point it was left.
 *
 * The return type of an asynchronous operation that uses use_awaitable
 * as completion token is a boost::asio::awaitable<T>, where T
 * is the second argument to the handler signature for the asynchronous operation.
 * If any of the asynchronous operations fail, an exception will be raised
 * within the coroutine.
 */
boost::asio::awaitable<void> coro_main(
    boost::mysql::tcp_ssl_connection& conn,
    boost::asio::ip::tcp::resolver& resolver,
    const boost::mysql::handshake_params& params,
    const char* hostname,
    const char* company_id
)
{
    error_code ec;
    boost::mysql::diagnostics diag;

    // Resolve hostname. We may use use_awaitable here, as hostname resolution
    // never produces any diagnostics.
    auto endpoints = co_await resolver.async_resolve(
        hostname,
        boost::mysql::default_port_string,
        boost::asio::use_awaitable
    );

    // Connect to server
    std::tie(ec) = co_await conn.async_connect(*endpoints.begin(), params, diag, tuple_awaitable);
    boost::mysql::throw_on_error(ec, diag);

    // We will be using company_id, which is untrusted user input, so we will use a prepared
    // statement.
    boost::mysql::statement stmt;
    std::tie(ec, stmt) = co_await conn.async_prepare_statement(
        "SELECT first_name, last_name, salary FROM employee WHERE company_id = ?",
        diag,
        tuple_awaitable
    );
    boost::mysql::throw_on_error(ec, diag);

    // Execute the statement
    boost::mysql::results result;
    std::tie(ec) = co_await conn.async_execute(stmt.bind(company_id), result, diag, tuple_awaitable);
    boost::mysql::throw_on_error(ec, diag);

    // Print all employees
    for (boost::mysql::row_view employee : result.rows())
    {
        print_employee(employee);
    }

    // Notify the MySQL server we want to quit, then close the underlying connection.
    std::tie(ec) = co_await conn.async_close(diag, tuple_awaitable);
    boost::mysql::throw_on_error(ec, diag);
}

void main_impl(int argc, char** argv)
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> [company-id]\n";
        exit(1);
    }

    const char* hostname = argv[3];

    // The company_id whose employees we will be listing. This
    // is user-supplied input, and should be treated as untrusted.
    const char* company_id = argc == 5 ? argv[4] : "HGS";

    // I/O context and connection. We use SSL because MySQL 8+ default settings require it.
    boost::asio::io_context ctx;
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    boost::mysql::tcp_ssl_connection conn(ctx, ssl_ctx);

    // Connection parameters
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database to use; leave empty or omit the parameter for no
                                // database
    );

    // Resolver for hostname resolution
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());

    // The entry point. We pass in a function returning
    // boost::asio::awaitable<void>, as required.
    boost::asio::co_spawn(
        ctx.get_executor(),
        [&conn, &resolver, params, hostname, company_id] {
            return coro_main(conn, resolver, params, hostname, company_id);
        },
        // If any exception is thrown in the coroutine body, rethrow it.
        [](std::exception_ptr ptr) {
            if (ptr)
            {
                std::rethrow_exception(ptr);
            }
        }
    );

    // Calling run will execute the requested operations.
    ctx.run();
}

#else

void main_impl(int, char**)
{
    std::cout << "Sorry, your compiler does not support C++20 coroutines" << std::endl;
}

#endif

int main(int argc, char** argv)
{
    try
    {
        main_impl(argc, argv);
    }
    catch (const boost::mysql::error_with_diagnostics& err)
    {
        // You will only get this type of exceptions if you use throw_on_error.
        // Some errors include additional diagnostics, like server-provided error messages.
        // Security note: diagnostics::server_message may contain user-supplied values (e.g. the
        // field value that caused the error) and is encoded using to the connection's encoding
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
