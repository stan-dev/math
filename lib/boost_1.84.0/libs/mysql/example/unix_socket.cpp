//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_unix_socket

#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/unix.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/local/stream_protocol.hpp>

#include <iostream>
#include <tuple>

void print_employee(boost::mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

#define ASSERT(expr)                                          \
    if (!(expr))                                              \
    {                                                         \
        std::cerr << "Assertion failed: " #expr << std::endl; \
        exit(1);                                              \
    }

/**
 * UNIX sockets are only available on, er, UNIX systems. Typedefs for
 * UNIX socket-based connections are only available in UNIX systems.
 * Check for BOOST_ASIO_HAS_LOCAL_SOCKETS to know if UNIX socket
 * typedefs are available in your system.
 */
#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS

void main_impl(int argc, char** argv)
{
    if (argc != 3 && argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> [<socket-path>] [<company-id>]\n";
        exit(1);
    }

    const char* socket_path = argc >= 4 ? argv[3] : "/var/run/mysqld/mysqld.sock";

    // The company_id whose employees we will be listing. This
    // is user-supplied input, and should be treated as untrusted.
    const char* company_id = argc == 5 ? argv[4] : "HGS";

    // Connection parameters that tell us where and how to connect to the MySQL server.
    // There are two types of parameters:
    //   - UNIX-level connection parameters, identifying the UNIX socket to connect to.
    //   - MySQL level parameters: database credentials and schema to use.
    boost::asio::local::stream_protocol::endpoint ep(socket_path);
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database to use; leave empty or omit the parameter for no
                                // database
    );

    boost::asio::io_context ctx;

    // Connection to the MySQL server, over a UNIX socket. Note that we don't need
    // to use SSL when using UNIX sockets because it's restricted to the local machine,
    // so MySQL considers it secure, even if it's not encrypted.
    boost::mysql::unix_connection conn(ctx);
    conn.connect(ep, params);  // UNIX socket connect and MySQL handshake

    // We will be using company_id, which is untrusted user input, so we will use a prepared
    // statement.
    boost::mysql::statement stmt = conn.prepare_statement(
        "SELECT first_name, last_name, salary FROM employee WHERE company_id = ?"
    );

    // Execute the statement
    boost::mysql::results result;
    conn.execute(stmt.bind(company_id), result);

    // Print employees
    for (boost::mysql::row_view employee : result.rows())
    {
        print_employee(employee);
    }

    // Notify the MySQL server we want to quit, then close the underlying connection.
    conn.close();
}

#else

void main_impl(int, char**) { std::cout << "Sorry, your system does not support UNIX sockets" << std::endl; }

#endif

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
        // field value that caused the error) and is encoded using to the connection's encoding
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
