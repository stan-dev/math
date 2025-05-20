//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/static_results.hpp>

#include "log_error.hpp"

#ifdef BOOST_MYSQL_CXX14

//[example_connection_pool_server_cpp
//
// File: server.cpp
//

#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/ip/address.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/spawn.hpp>
#include <boost/asio/strand.hpp>
#include <boost/beast/core/flat_buffer.hpp>
#include <boost/beast/http/error.hpp>
#include <boost/beast/http/message.hpp>
#include <boost/beast/http/parser.hpp>
#include <boost/beast/http/read.hpp>
#include <boost/beast/http/string_body.hpp>
#include <boost/beast/http/write.hpp>

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include "handle_request.hpp"
#include "repository.hpp"
#include "server.hpp"
#include "types.hpp"

// This file contains all the boilerplate code to implement a HTTP
// server. Functions here end up invoking handle_request.

namespace asio = boost::asio;
namespace http = boost::beast::http;
using boost::mysql::error_code;
using namespace notes;

namespace {

static void run_http_session(
    boost::asio::ip::tcp::socket sock,
    std::shared_ptr<shared_state> st,
    boost::asio::yield_context yield
)
{
    error_code ec;

    // A buffer to read incoming client requests
    boost::beast::flat_buffer buff;

    while (true)
    {
        // Construct a new parser for each message
        http::request_parser<http::string_body> parser;

        // Apply a reasonable limit to the allowed size
        // of the body in bytes to prevent abuse.
        parser.body_limit(10000);

        // Read a request
        http::async_read(sock, buff, parser.get(), yield[ec]);

        if (ec)
        {
            if (ec == http::error::end_of_stream)
            {
                // This means they closed the connection
                sock.shutdown(asio::ip::tcp::socket::shutdown_send, ec);
            }
            else
            {
                // An unknown error happened
                log_error("Error reading HTTP request: ", ec);
            }
            return;
        }

        // Process the request to generate a response.
        // This invokes the business logic, which will need to access MySQL data
        auto response = handle_request(parser.get(), note_repository(st->pool), yield);

        // Determine if we should close the connection
        bool keep_alive = response.keep_alive();

        // Send the response
        http::async_write(sock, response, yield[ec]);
        if (ec)
            return log_error("Error writing HTTP response: ", ec);

        // This means we should close the connection, usually because
        // the response indicated the "Connection: close" semantic.
        if (!keep_alive)
        {
            sock.shutdown(asio::ip::tcp::socket::shutdown_send, ec);
            return;
        }
    }
}

// Implements the server's accept loop. The server will
// listen for connections until stopped.
static void do_accept(
    asio::any_io_executor executor,  // The original executor (without strands)
    std::shared_ptr<asio::ip::tcp::acceptor> acceptor,
    std::shared_ptr<shared_state> st
)
{
    acceptor->async_accept([executor, st, acceptor](error_code ec, asio::ip::tcp::socket sock) {
        // If there was an error accepting the connection, exit our loop
        if (ec)
            return log_error("Error while accepting connection", ec);

        // Launch a new session for this connection. Each session gets its
        // own stackful coroutine, so we can get back to listening for new connections.
        boost::asio::spawn(
            // Every session gets its own strand. This prevents data races.
            asio::make_strand(executor),

            // The actual coroutine
            [st, socket = std::move(sock)](boost::asio::yield_context yield) mutable {
                run_http_session(std::move(socket), std::move(st), yield);
            },

            // All errors in the session are handled via error codes or by catching
            // exceptions explicitly. An unhandled exception here means an error.
            // Rethrowing it will propagate the exception, making io_context::run()
            // to throw and terminate the program.
            [](std::exception_ptr ex) {
                if (ex)
                    std::rethrow_exception(ex);
            }
        );

        // Accept a new connection
        do_accept(executor, acceptor, st);
    });
}

}  // namespace

error_code notes::launch_server(
    boost::asio::any_io_executor ex,
    std::shared_ptr<shared_state> st,
    unsigned short port
)
{
    error_code ec;

    // An object that allows us to accept incoming TCP connections.
    // Since we're in a multi-threaded environment, we create a strand for the acceptor,
    // so all accept handlers are run serialized
    auto acceptor = std::make_shared<asio::ip::tcp::acceptor>(asio::make_strand(ex));

    // The endpoint where the server will listen. Edit this if you want to
    // change the address or port we bind to.
    boost::asio::ip::tcp::endpoint listening_endpoint(boost::asio::ip::make_address("0.0.0.0"), port);

    // Open the acceptor
    acceptor->open(listening_endpoint.protocol(), ec);
    if (ec)
        return ec;

    // Allow address reuse
    acceptor->set_option(asio::socket_base::reuse_address(true), ec);
    if (ec)
        return ec;

    // Bind to the server address
    acceptor->bind(listening_endpoint, ec);
    if (ec)
        return ec;

    // Start listening for connections
    acceptor->listen(asio::socket_base::max_listen_connections, ec);
    if (ec)
        return ec;

    std::cout << "Server listening at " << acceptor->local_endpoint() << std::endl;

    // Launch the acceptor loop
    do_accept(std::move(ex), std::move(acceptor), std::move(st));

    // Done
    return error_code();
}

//]

#endif
