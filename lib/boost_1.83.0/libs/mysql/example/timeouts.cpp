//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_timeouts

#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
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
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/use_awaitable.hpp>

#include <chrono>
#include <exception>
#include <iostream>
#include <stdexcept>

#ifdef BOOST_ASIO_HAS_CO_AWAIT

#include <boost/asio/experimental/awaitable_operators.hpp>

using namespace boost::asio::experimental::awaitable_operators;
using boost::asio::use_awaitable;
using boost::mysql::error_code;

constexpr std::chrono::milliseconds TIMEOUT(8000);

void print_employee(boost::mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

/**
 * Helper functions to check whether an async operation, launched in parallel with
 * a timer, was successful, resulted in an error or timed out. The timer is always the first operation.
 * If the variant holds the first alternative, the timer fired before
 * the async operation completed, which means a timeout. We'll be using as_tuple with use_awaitable to be able
 * to use boost::mysql::throw_on_error and include server diagnostics in the thrown exceptions.
 */
template <class T>
T check_error(
    std::variant<std::monostate, std::tuple<error_code, T>>&& op_result,
    const boost::mysql::diagnostics& diag = {}
)
{
    if (op_result.index() == 0)
    {
        throw std::runtime_error("Operation timed out");
    }
    auto [ec, res] = std::get<1>(std::move(op_result));
    boost::mysql::throw_on_error(ec, diag);
    return res;
}

void check_error(
    const std::variant<std::monostate, std::tuple<error_code>>& op_result,
    const boost::mysql::diagnostics& diag
)
{
    if (op_result.index() == 0)
    {
        throw std::runtime_error("Operation timed out");
    }
    auto [ec] = std::get<1>(op_result);
    boost::mysql::throw_on_error(ec, diag);
}

// Using this completion token instead of plain use_awaitable prevents
// co_await from throwing exceptions. Instead, co_await will return a std::tuple<error_code>
// with a non-zero code on error. We will then use boost::mysql::throw_on_error
// to throw exceptions with embedded diagnostics, if available. If you
// employ plain use_awaitable, you will get boost::system::system_error exceptions
// instead of boost::mysql::error_with_diagnostics exceptions. This is a limitation of use_awaitable.
constexpr auto tuple_awaitable = boost::asio::as_tuple(boost::asio::use_awaitable);

/**
 * We use Boost.Asio's cancellation capabilities to implement timeouts for our
 * asynchronous operations. This is not something specific to Boost.MySQL, and
 * can be used with any other asynchronous operation that follows Asio's model.
 *
 * Each time we invoke an asynchronous operation, we also call timer_type::async_wait.
 * We then use Asio's overload for operator || to run the timer wait and the async operation
 * in parallel. Once the first of them finishes, the other operation is cancelled
 * (the behavior is similar to JavaScripts's Promise.race).
 * If we co_await the awaitable returned by operator ||, we get a std::variant<std::monostate, T>,
 * where T is the async operation's result type. If the timer wait finishes first (we have a
 * timeout), the variant will hold the std::monostate at index 0; otherwise, it will have the async
 * operation's result at index 1. The function check_error throws an exception in the case of
 * timeout and extracts the operation's result otherwise.
 *
 * If any of the MySQL specific operations result in a timeout, the connection is left
 * in an unspecified state. You should close it and re-open it to get it working again.
 */
boost::asio::awaitable<void> coro_main(
    boost::mysql::tcp_ssl_connection& conn,
    boost::asio::ip::tcp::resolver& resolver,
    boost::asio::steady_timer& timer,
    const boost::mysql::handshake_params& params,
    const char* hostname,
    const char* company_id
)
{
    boost::mysql::diagnostics diag;

    // Resolve hostname
    timer.expires_after(TIMEOUT);
    auto endpoints = check_error(co_await (
        timer.async_wait(use_awaitable) ||
        resolver.async_resolve(hostname, boost::mysql::default_port_string, tuple_awaitable)
    ));

    // Connect to server. Note that we need to reset the timer before using it again.
    timer.expires_after(TIMEOUT);
    auto op_result = co_await (
        timer.async_wait(use_awaitable) ||
        conn.async_connect(*endpoints.begin(), params, diag, tuple_awaitable)
    );
    check_error(op_result, diag);

    // We will be using company_id, which is untrusted user input, so we will use a prepared
    // statement.
    auto stmt_op_result = co_await (
        timer.async_wait(use_awaitable) ||
        conn.async_prepare_statement(
            "SELECT first_name, last_name, salary FROM employee WHERE company_id = ?",
            diag,
            tuple_awaitable
        )
    );
    boost::mysql::statement stmt = check_error(std::move(stmt_op_result), diag);

    // Execute the statement
    boost::mysql::results result;
    timer.expires_after(TIMEOUT);
    op_result = co_await (
        timer.async_wait(use_awaitable) ||
        conn.async_execute(stmt.bind(company_id), result, diag, tuple_awaitable)
    );
    check_error(op_result, diag);

    // Print all the obtained rows
    for (boost::mysql::row_view employee : result.rows())
    {
        print_employee(employee);
    }

    // Notify the MySQL server we want to quit, then close the underlying connection.
    op_result = co_await (timer.async_wait(use_awaitable) || conn.async_close(diag, tuple_awaitable));
    check_error(op_result, diag);
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
    boost::asio::steady_timer timer(ctx.get_executor());

    // Connection parameters
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database to use; leave empty or omit for no database
    );

    // Resolver for hostname resolution
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());

    // The entry point. We pass in a function returning a boost::asio::awaitable<void>, as required.
    boost::asio::co_spawn(
        ctx.get_executor(),
        [&conn, &resolver, &timer, params, hostname, company_id] {
            return coro_main(conn, resolver, timer, params, hostname, company_id);
        },
        // If any exception is thrown in the coroutine body, rethrow it.
        [](std::exception_ptr ptr) {
            if (ptr)
            {
                std::rethrow_exception(ptr);
            }
        }
    );

    // Calling run will actually start the requested operations.
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
