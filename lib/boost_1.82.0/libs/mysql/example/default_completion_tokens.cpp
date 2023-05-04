//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_default_completion_tokens]

#include <boost/mysql.hpp>

#include <boost/asio/as_tuple.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/use_awaitable.hpp>

#include <iostream>

using boost::mysql::error_code;

#ifdef BOOST_ASIO_HAS_CO_AWAIT

void print_employee(boost::mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

// In Boost.Asio, default completion tokens are associated to executors.
// Instead of using the usual I/O object types (like tcp_ssl_connection), we need
// instantiations of the base template that use an executor with a default
// completion token. We can achieve this using CompletionToken::as_default_on_t<IOObject>,
// where CompletionToken is the desired default token and IOObject is the usual I/O object type.
// as_default_on_t requires the I/O object to implement rebind_executor.
using tuple_awaitable_t = boost::asio::as_tuple_t<boost::asio::use_awaitable_t<>>;
using resolver_type = tuple_awaitable_t::as_default_on_t<boost::asio::ip::tcp::resolver>;
using connection_type = tuple_awaitable_t::as_default_on_t<boost::mysql::tcp_ssl_connection>;

// Our coroutine
boost::asio::awaitable<void> coro_main(
    connection_type& conn,
    resolver_type& resolver,
    const char* hostname,
    const boost::mysql::handshake_params& params,
    const char* company_id
)
{
    boost::mysql::diagnostics diag;

    // Resolve hostname
    auto [ec, endpoints] = co_await resolver.async_resolve(hostname, boost::mysql::default_port_string);
    boost::mysql::throw_on_error(ec);

    // Connect to server
    std::tie(ec) = co_await conn.async_connect(*endpoints.begin(), params, diag);
    boost::mysql::throw_on_error(ec, diag);

    // Prepare an statement
    boost::mysql::statement stmt;
    std::tie(ec, stmt) = co_await conn.async_prepare_statement(
        "SELECT first_name, last_name, salary FROM employee WHERE company_id = ?",
        diag
    );
    boost::mysql::throw_on_error(ec, diag);

    // Execute it
    boost::mysql::results result;
    std::tie(ec) = co_await conn.async_execute_statement(stmt, std::make_tuple(company_id), result, diag);
    boost::mysql::throw_on_error(ec, diag);

    // Use the received rows
    for (boost::mysql::row_view employee : result.rows())
    {
        print_employee(employee);
    }

    // Notify the MySQL server we want to quit, then close the underlying connection.
    // This will also deallocate the statement from the server.
    std::tie(ec) = co_await conn.async_close(diag);
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
    const char* company_id = argc == 5 ? argv[4] : "HGS";

    // I/O context and connection. We use SSL because MySQL 8+ default settings require it.
    boost::asio::io_context ctx;
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    connection_type conn(ctx.get_executor(), ssl_ctx);

    // Resolver, for hostname resolution
    resolver_type resolver(ctx.get_executor());

    // Connection parameters
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database to use; leave empty or omit for no database
    );

    // Spawn the coroutine
    boost::asio::co_spawn(
        ctx.get_executor(),
        [&conn, &resolver, hostname, params, company_id] {
            return coro_main(conn, resolver, hostname, params, company_id);
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
