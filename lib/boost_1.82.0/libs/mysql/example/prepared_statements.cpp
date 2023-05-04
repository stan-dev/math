//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_prepared_statements

#include <boost/mysql.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ssl/context.hpp>

#include <iostream>
#include <random>
#include <tuple>

#define ASSERT(expr)                                          \
    if (!(expr))                                              \
    {                                                         \
        std::cerr << "Assertion failed: " #expr << std::endl; \
        exit(1);                                              \
    }

double generate_random_payrise()
{
    std::random_device dev;
    std::default_random_engine eng(dev());
    std::uniform_real_distribution<> dist(500.0, 3000.0);
    return dist(eng);
}

void main_impl(int argc, char** argv)
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <username> <password> <server-hostname> [employee-first-name]\n";
        exit(1);
    }

    // The first_name of the employee we will be working with. This
    // is user-supplied input, and should be treated as untrusted.
    const char* first_name = argc == 5 ? argv[4] : "Efficient";

    // I/O context and connection. We use SSL because MySQL 8+ default settings require it.
    boost::asio::io_context ctx;
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    boost::mysql::tcp_ssl_connection conn(ctx, ssl_ctx);

    // Resolver for hostname resolution
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());

    // Connection params
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database to use; leave empty or omit the parameter for no
                                // database
    );

    // Hostname resolution
    auto endpoints = resolver.resolve(argv[3], boost::mysql::default_port_string);

    // TCP and MySQL level connect
    conn.connect(*endpoints.begin(), params);

    /**
     * We can tell MySQL to prepare a statement using connection::prepare_statement.
     * We provide a string SQL statement, which can include any number of parameters,
     * identified by question marks.
     *
     * Prepared statements are stored in the server on a per-connection basis.
     * Once a connection is closed, all prepared statements for that connection are deallocated.
     *
     * The result of prepare_statement is a boost::mysql::statement object, which
     * is a lightweight handle for the server-side statement.
     *
     * We prepare two statements, a SELECT and an UPDATE.
     */
    //[prepared_statements_prepare
    boost::mysql::statement salary_getter = conn.prepare_statement(
        "SELECT salary FROM employee WHERE first_name = ?"
    );
    //]

    // num_params() returns the number of parameters (question marks)
    ASSERT(salary_getter.num_params() == 1);

    boost::mysql::statement salary_updater = conn.prepare_statement(
        "UPDATE employee SET salary = salary + ? WHERE first_name = ?"
    );
    ASSERT(salary_updater.num_params() == 2);

    /*
     * Once a statement has been prepared, it can be executed by calling
     * connection::execute_statement(). Parameter actual values are provided
     * as a std::tuple. Executing a statement yields a results object.
     */
    //[prepared_statements_execute
    boost::mysql::results select_result, update_result;
    conn.execute_statement(salary_getter, std::make_tuple(first_name), select_result);
    //]

    // First row, first column, cast to double
    double old_salary = select_result.rows().at(0).at(0).as_double();
    std::cout << "The salary before the payrise was: " << old_salary << std::endl;

    // Run the update. In this case, we must pass in two parameters.
    double payrise = generate_random_payrise();
    conn.execute_statement(salary_updater, std::make_tuple(payrise, first_name), update_result);
    ASSERT(update_result.rows().empty());  // an UPDATE never returns rows

    /**
     * Execute the select again. We can execute a prepared statement
     * as many times as we want. We do NOT need to call
     * connection::prepare_statement() again.
     */
    conn.execute_statement(salary_getter, std::make_tuple(first_name), select_result);
    double new_salary = select_result.rows().at(0).at(0).as_double();
    ASSERT(new_salary > old_salary);  // Our update took place
    std::cout << "The salary after the payrise was: " << new_salary << std::endl;

    /**
     * Close the statements. Closing a statement deallocates it from the server.
     * Closing statements implies communicating with the server and can thus fail.
     *
     * Statements are automatically deallocated once the connection is closed.
     * If you are re-using connection objects and preparing statements over time,
     * you should close your statements to prevent excessive resource usage.
     * If you are not re-using the connections, or are preparing your statements
     * just once at application startup, there is no need to perform this step.
     */
    conn.close_statement(salary_updater);
    conn.close_statement(salary_getter);

    // Close the connection
    conn.close();
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
