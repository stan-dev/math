//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_text_queries

#include <boost/mysql.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ssl/context.hpp>

#include <iostream>

#define ASSERT(expr)                                          \
    if (!(expr))                                              \
    {                                                         \
        std::cerr << "Assertion failed: " #expr << std::endl; \
        exit(1);                                              \
    }

/**
 * Prints an employee to std::cout. An employee here is a boost::mysql::row_view,
 * which represents a row returned by a SQL query. row_view objects are an ordered
 * collection of SQL fields, representing each value returned by the query.
 *
 * Indexing a row_view yields a boost::mysql::field_view, which is a variant-like
 * type representing a single value returned by MySQL.
 */
void print_employee(boost::mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

void main_impl(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname>\n";
        exit(1);
    }

    // The I/O context to perform all operations.
    boost::asio::io_context ctx;

    /**
     * Connection parameters that tell us how to connect to the MySQL server:
     * database credentials and schema to use.
     */
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database to use; leave empty or omit for no database
    );

    /* We will use SSL in all our examples. To enable SSL, use boost::mysql::tcp_ssl_connection.
     * MySQL 8+ default is to use an authentication method that requires SSL, so we encourage
     * you to use SSL connections if you can.
     */
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);

    // Represents a single connection over TCP to a MySQL server.
    boost::mysql::tcp_ssl_connection conn(ctx, ssl_ctx);

    // To establish the connection, we need a TCP endpoint. We have a hostname,
    // so we need to perform hostname resolution. We create a resolver for this.
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());

    // Invoke the resolver's appropriate function to perform the resolution.
    const char* hostname = argv[3];
    auto endpoints = resolver.resolve(hostname, boost::mysql::default_port_string);

    /**
     * Before using the connection, we have to connect to the server by:
     *  - Establishing the TCP-level session.
     *  - Authenticating to the MySQL server. The SSL handshake is performed as part of this.
     *    connection::connect takes care of both.
     */
    conn.connect(*endpoints.begin(), params);

    /**
     * To issue a SQL query to the database server, use tcp_ssl_connection::query, which takes
     * the SQL to be executed as parameter and returns a results object by lvalue reference.
     * Resultset objects contain the retrieved rows, among other info.
     * We will get all employees working for 'High Growth Startup'.
     */
    const char* sql = "SELECT first_name, last_name, salary FROM employee WHERE company_id = 'HGS'";
    boost::mysql::results result;
    conn.query(sql, result);

    // We can access the rows using results::rows
    for (boost::mysql::row_view employee : result.rows())
    {
        print_employee(employee);
    }

    // We can issue any SQL statement, not only SELECTs. In this case, the returned
    // results will have no fields and no rows
    sql = "UPDATE employee SET salary = 10000 WHERE first_name = 'Underpaid'";
    conn.query(sql, result);
    ASSERT(result.rows().empty());  // UPDATEs don't retrieve rows

    // Check we have updated our poor intern salary
    conn.query("SELECT salary FROM employee WHERE first_name = 'Underpaid'", result);
    double salary = result.rows().at(0).at(0).as_double();
    ASSERT(salary == 10000.0);

    // Close the connection. This notifies the MySQL we want to log out
    // and then closes the underlying socket. This operation implies a network
    // transfer and thus can fail
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
