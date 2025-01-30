//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_callbacks

/**
 * This example demonstrates how to use callbacks when using async functions.
 * This can be a good choice when targeting a standard lower than C++20.
 * This example uses the 'boost_mysql_examples' database, which you
 * can get by running db_setup.sql.
 */

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/system/error_code.hpp>

#include <functional>
#include <iostream>
#include <memory>

// When using callbacks, we usually employ error codes instead of exceptions.
using boost::system::error_code;

namespace mysql = boost::mysql;
namespace asio = boost::asio;

// Prints a database employee to stdout
void print_employee(mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

// A session object, containing all variables that need to be kept alive for our session.
// We will use a shared_ptr to ensure that all these variables are kept alive
// until the last callback is executed
class session : public std::enable_shared_from_this<session>
{
    mysql::connect_params conn_params;  // MySQL credentials and other connection config
    mysql::any_connection conn;         // Represents the connection to the MySQL server
    mysql::results result;              // A result from a query
    mysql::error_code final_error;      // Will be set in case of error
    mysql::diagnostics diag;            // Will be populated with info about server errors
    const char* company_id;             // The ID of the company whose employees we want to list. Untrusted.
public:
    session(
        asio::io_context& ctx,
        const char* server_hostname,
        const char* username,
        const char* password,
        const char* company_id
    )
        : conn(ctx), company_id(company_id)
    {
        conn_params.server_address.emplace_host_and_port(server_hostname);
        conn_params.username = username;
        conn_params.password = password;
        conn_params.database = "boost_mysql_examples";
    }

    // Accessor for error information, so main can access it
    error_code get_error() const { return final_error; }
    const boost::mysql::diagnostics& get_diagnostics() const { return diag; }

    // Initiates the callback chain
    void start()
    {
        // Will call on_connect when the connect operation completes.
        // The session object is kept alive with the shared_ptr that shared_from_this produces
        conn.async_connect(
            conn_params,
            diag,
            std::bind(&session::on_connect, shared_from_this(), std::placeholders::_1)
        );
    }

    void on_connect(error_code ec)
    {
        // If there was an error, stop the callback chain
        if (ec)
        {
            final_error = ec;
            return;
        }

        // Initiate the query execution. company_id is an untrusted value.
        // with_params will securely compose a SQL query and send it to the server for execution.
        // Returned rows will be read into result.
        // We use the callback chain + shared_ptr technique again
        conn.async_execute(
            mysql::with_params(
                "SELECT first_name, last_name, salary FROM employee WHERE company_id = {}",
                company_id
            ),
            result,
            diag,
            std::bind(&session::on_execute, shared_from_this(), std::placeholders::_1)
        );
    }

    void on_execute(error_code ec)
    {
        // If there was an error, stop the callback chain
        if (ec)
        {
            final_error = ec;
            return;
        }

        // Print the rows returned by the query
        for (boost::mysql::row_view employee : result.rows())
        {
            print_employee(employee);
        }

        // Notify the MySQL server we want to quit and then close the socket
        conn.async_close(diag, std::bind(&session::finish, shared_from_this(), std::placeholders::_1));
    }

    void finish(error_code err) { final_error = err; }
};

void main_impl(int argc, char** argv)
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> [company-id]\n";
        exit(1);
    }

    // The execution context, required to run I/O operations.
    boost::asio::io_context ctx;

    // The company_id whose employees we will be listing. This
    // is user-supplied input, and should be treated as untrusted.
    const char* company_id = argc == 5 ? argv[4] : "HGS";

    // Create the session object and launch it
    auto sess = std::make_shared<session>(ctx, argv[3], argv[1], argv[2], company_id);
    sess->start();

    // Run the callback chain until it completes
    ctx.run();

    // Check for errors
    if (error_code ec = sess->get_error())
    {
        std::cerr << "Error: " << ec << ": " << ec.message() << '\n'
                  << "Server diagnostics: " << sess->get_diagnostics().server_message() << std::endl;
        exit(1);
    }
}

int main(int argc, char** argv)
{
    try
    {
        main_impl(argc, argv);
    }
    catch (const std::exception& err)
    {
        std::cerr << "Error: " << err.what() << std::endl;
        return 1;
    }
}

//]
