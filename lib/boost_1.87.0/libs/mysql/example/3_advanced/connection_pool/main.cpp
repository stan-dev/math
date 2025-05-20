//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/static_results.hpp>

#include "log_error.hpp"

#ifdef BOOST_MYSQL_CXX14

//[example_connection_pool_main_cpp

/**
 * This example demonstrates how to use a connection_pool.
 * It implements a minimal REST API to manage notes.
 * A note is a simple object containing a user-defined title and content.
 * The REST API offers CRUD operations on such objects:
 *    POST   /notes        Creates a new note.
 *    GET    /notes        Retrieves all notes.
 *    GET    /notes/<id>   Retrieves a single note.
 *    PUT    /notes/<id>   Replaces a note, changing its title and content.
 *    DELETE /notes/<id>   Deletes a note.
 *
 * Notes are stored in MySQL. The note_repository class encapsulates
 *   access to MySQL, offering friendly functions to manipulate notes.
 * server.cpp encapsulates all the boilerplate to launch an HTTP server,
 *   match URLs to API endpoints, and invoke the relevant note_repository functions.
 * All communication happens asynchronously. We use stackful coroutines to simplify
 *   development, using boost::asio::spawn and boost::asio::yield_context.
 * This example requires linking to Boost::context, Boost::json and Boost::url.
 */

#include <boost/mysql/any_address.hpp>
#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/pool_params.hpp>

#include <boost/asio/detached.hpp>
#include <boost/asio/signal_set.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/system/error_code.hpp>

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include "server.hpp"

using namespace notes;

// The number of threads to use
static constexpr std::size_t num_threads = 5;

int main(int argc, char* argv[])
{
    // Check command line arguments.
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <mysql-hostname> <port>\n";
        return EXIT_FAILURE;
    }

    // Application config
    const char* mysql_username = argv[1];
    const char* mysql_password = argv[2];
    const char* mysql_hostname = argv[3];
    auto port = static_cast<unsigned short>(std::stoi(argv[4]));

    // An event loop, where the application will run.
    // We will use the main thread to run the pool, too, so we use
    // one thread less than configured
    boost::asio::thread_pool th_pool(num_threads - 1);

    // Configuration for the connection pool
    boost::mysql::pool_params pool_prms{
        // Connect using TCP, to the given hostname and using the default port
        boost::mysql::host_and_port{mysql_hostname},

        // Authenticate using the given username
        mysql_username,

        // Password for the above username
        mysql_password,

        // Database to use when connecting
        "boost_mysql_examples",
    };

    // Using thread_safe will make the pool thread-safe by internally
    // creating and using a strand.
    // This allows us to share the pool between sessions, which may run
    // concurrently, on different threads.
    pool_prms.thread_safe = true;

    // Create the connection pool
    auto shared_st = std::make_shared<shared_state>(
        boost::mysql::connection_pool(th_pool, std::move(pool_prms))
    );

    // A signal_set allows us to intercept SIGINT and SIGTERM and
    // exit gracefully
    boost::asio::signal_set signals{th_pool.get_executor(), SIGINT, SIGTERM};

    // Launch the MySQL pool
    shared_st->pool.async_run(boost::asio::detached);

    // Start listening for HTTP connections. This will run until the context is stopped
    auto ec = launch_server(th_pool.get_executor(), shared_st, port);
    if (ec)
    {
        log_error("Error launching server: ", ec);
        exit(EXIT_FAILURE);
    }

    // Capture SIGINT and SIGTERM to perform a clean shutdown
    signals.async_wait([shared_st, &th_pool](boost::system::error_code, int) {
        // Cancel the pool. This will cause async_run to complete.
        shared_st->pool.cancel();

        // Stop the execution context. This will cause main to exit
        th_pool.stop();
    });

    // Attach the current thread to the thread pool. This will block
    // until stop() is called
    th_pool.attach();

    // Wait until all threads have exited
    th_pool.join();

    std::cout << "Server exiting" << std::endl;

    // (If we get here, it means we got a SIGINT or SIGTERM)
    return EXIT_SUCCESS;
}

//]

#else

#include <iostream>

int main() { std::cout << "Sorry, your compiler doesn't support C++14\n"; }

#endif
