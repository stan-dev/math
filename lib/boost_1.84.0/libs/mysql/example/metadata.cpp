//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_metadata

#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl/context.hpp>

#include <iostream>

#define ASSERT(expr)                                          \
    if (!(expr))                                              \
    {                                                         \
        std::cerr << "Assertion failed: " #expr << std::endl; \
        exit(1);                                              \
    }

void main_impl(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname>\n";
        exit(1);
    }

    // I/O context and connection. We use SSL because MySQL 8+ default settings require it.
    boost::asio::io_context ctx;
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    boost::mysql::tcp_ssl_connection conn(ctx, ssl_ctx);

    // By default, string metadata (like column names) won't be retained.
    // This is for efficiency reasons. You can change this setting by calling
    // connection::set_meta_mode. It will affect any subsequent queries and statement executions.
    conn.set_meta_mode(boost::mysql::metadata_mode::full);

    // Hostname resolution
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());
    auto endpoints = resolver.resolve(argv[3], boost::mysql::default_port_string);

    // TCP and MySQL level connect
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database to use; leave empty or omit the parameter for no
                                // database
    );
    conn.connect(*endpoints.begin(), params);

    // Issue the query
    const char* sql = R"(
        SELECT comp.name AS company_name, emp.id AS employee_id
        FROM employee emp
        JOIN company comp ON (comp.id = emp.company_id)
    )";
    boost::mysql::results result;
    conn.execute(sql, result);

    /**
     * results objects allow you to access metadata about the columns in the query
     * using the meta() function, which returns span-like object containing metadata objects
     * (one per column in the query, and in the same order as in the query).
     * You can retrieve the column name, type, number of decimals,
     * suggested display width, whether the column is part of a key...
     * These metadata objects are owned by the results object.
     */
    ASSERT(result.meta().size() == 2);

    // clang-format off
    const boost::mysql::metadata& company_name = result.meta()[0];
    ASSERT(company_name.database() == "boost_mysql_examples");  // database name
    ASSERT(company_name.table() == "comp");  // the alias we assigned to the table in the query
    ASSERT(company_name.original_table() == "company");   // the original table name
    ASSERT(company_name.column_name() == "company_name");  // the name of the column in the query
    ASSERT(company_name.original_column_name() == "name");  // the name of the physical column in the table
    ASSERT(company_name.type() == boost::mysql::column_type::varchar);  // we created the column as a VARCHAR
    ASSERT(!company_name.is_primary_key());     // column is not a primary key
    ASSERT(!company_name.is_auto_increment());  // column is not AUTO_INCREMENT
    ASSERT(company_name.is_not_null());         // column may not be NULL

    const boost::mysql::metadata& employee_id = result.meta()[1];
    ASSERT(employee_id.database() == "boost_mysql_examples");  // database name
    ASSERT(employee_id.table() == "emp");  // the alias we assigned to the table in the query
    ASSERT(employee_id.original_table() == "employee");  // the original table name
    ASSERT(employee_id.column_name() == "employee_id");   // the name of the column in the query
    ASSERT(employee_id.original_column_name() == "id");  // the name of the physical column in the table
    ASSERT(employee_id.type() == boost::mysql::column_type::int_);  // we created the column as INT
    ASSERT(employee_id.is_primary_key()); // column is a primary key
    ASSERT(employee_id.is_auto_increment()); // we declared the column as AUTO_INCREMENT
    ASSERT(employee_id.is_not_null()); // column cannot be NULL
    // clang-format on

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
