//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_async_callbacks

#include <boost/mysql/connection.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/tcp_ssl.hpp>
#include <boost/mysql/throw_on_error.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl/context.hpp>

#include <iostream>

using boost::mysql::error_code;

void print_employee(boost::mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

class application
{
    boost::asio::ip::tcp::resolver::results_type eps;  // Physical endpoint(s) to connect to
    boost::mysql::handshake_params conn_params;        // MySQL credentials and other connection config
    boost::asio::io_context ctx;                       // boost::asio context
    boost::asio::ip::tcp::resolver resolver;           // To perform hostname resolution
    boost::asio::ssl::context ssl_ctx;                 // MySQL 8+ default settings require SSL
    boost::mysql::tcp_ssl_connection conn;             // Represents the connection to the MySQL server
    boost::mysql::statement stmt;                      // A prepared statement
    boost::mysql::results result;                      // A result from a query
    boost::mysql::diagnostics diag;                    // Will be populated with info about server errors
    const char* company_id;  // The ID of the company whose employees we want to list. Untrusted.
public:
    application(const char* username, const char* password, const char* company_id)
        : conn_params(username, password, "boost_mysql_examples"),
          resolver(ctx.get_executor()),
          ssl_ctx(boost::asio::ssl::context::tls_client),
          conn(ctx, ssl_ctx),
          company_id(company_id)
    {
    }

    void start(const char* hostname) { resolve_hostname(hostname); }

    void resolve_hostname(const char* hostname)
    {
        resolver.async_resolve(
            hostname,
            boost::mysql::default_port_string,
            [this](error_code err, boost::asio::ip::tcp::resolver::results_type results) {
                boost::mysql::throw_on_error(err);
                eps = std::move(results);
                connect();
            }
        );
    }

    void connect()
    {
        conn.async_connect(*eps.begin(), conn_params, diag, [this](error_code err) {
            boost::mysql::throw_on_error(err, diag);
            prepare_statement();
        });
    }

    void prepare_statement()
    {
        // We will be using company_id, which is untrusted user input, so we will use a prepared
        // statement.
        conn.async_prepare_statement(
            "SELECT first_name, last_name, salary FROM employee WHERE company_id = ?",
            diag,
            [this](error_code err, boost::mysql::statement temp_stmt) {
                boost::mysql::throw_on_error(err, diag);
                stmt = temp_stmt;
                query_employees();
            }
        );
    }

    void query_employees()
    {
        conn.async_execute(stmt.bind(company_id), result, diag, [this](error_code err) {
            boost::mysql::throw_on_error(err, diag);
            for (boost::mysql::row_view employee : result.rows())
            {
                print_employee(employee);
            }
            close();
        });
    }

    void close()
    {
        // Notify the MySQL server we want to quit and then close the socket
        conn.async_close(diag, [this](error_code err) { boost::mysql::throw_on_error(err, diag); });
    }

    void run() { ctx.run(); }
};

void main_impl(int argc, char** argv)
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> [company-id]\n";
        exit(1);
    }

    // The company_id whose employees we will be listing. This
    // is user-supplied input, and should be treated as untrusted.
    const char* company_id = argc == 5 ? argv[4] : "HGS";

    application app(argv[1], argv[2], company_id);
    app.start(argv[3]);  // starts the async chain
    app.run();           // run the asio::io_context until the async chain finishes
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
