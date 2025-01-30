//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_source_script

#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/metadata_collection_view.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/system/system_error.hpp>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

/**
 * This example runs all command in a .sql file, using multi-queries.
 * Note that special commands that are handled by the mysql command line tool
 * (like DELIMITER) won't work.
 *
 * For this example, we will be using the 'boost_mysql_examples' database.
 * You can get this database by running db_setup.sql.
 * This example assumes you are connecting to a localhost MySQL server.
 *
 * This example uses synchronous functions and handles errors using exceptions.
 */

// Reads a file into memory
std::string read_file(const char* file_name)
{
    std::ifstream ifs(file_name);
    if (!ifs)
        throw std::runtime_error("Cannot open file: " + std::string(file_name));
    return std::string(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
}

void print_column_names(boost::mysql::metadata_collection_view meta_collection)
{
    if (meta_collection.empty())
        return;

    bool is_first = true;
    for (auto meta : meta_collection)
    {
        if (!is_first)
        {
            std::cout << " | ";
        }
        is_first = false;
        std::cout << meta.column_name();
    }
    std::cout << "\n-----------------\n";
}

void print_row(boost::mysql::row_view row)
{
    bool is_first = true;
    for (auto field : row)
    {
        if (!is_first)
        {
            std::cout << " | ";
        }
        is_first = false;
        std::cout << field;
    }
    std::cout << '\n';
}

void print_ok(const boost::mysql::execution_state& st)
{
    std::cout << "Affected rows: " << st.affected_rows()
              << "\n"
                 "Last insert ID: "
              << st.last_insert_id()
              << "\n"
                 "Warnings: "
              << st.warning_count() << "\n\n"
              << std::flush;
}

void main_impl(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> <path-to-script>\n";
        exit(1);
    }

    // Read the script file into memory
    std::string script_contents = read_file(argv[4]);

    // Set up the io_context, SSL context and connection required to
    // connect to the server.
    boost::asio::io_context ctx;
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    boost::mysql::tcp_ssl_connection conn(ctx.get_executor(), ssl_ctx);

    // Resolve the server hostname to get a collection of endpoints
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());
    auto endpoints = resolver.resolve(argv[3], boost::mysql::default_port_string);

    // The username, password and database to use
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database
    );

    // We're going to use multi-queries, which enables passing the server
    // a set of semicolon-separated queries. We need to explicitly enable support for it.
    params.set_multi_queries(true);

    // We'll be using metadata strings to print column names, so we need to enable support for it
    conn.set_meta_mode(boost::mysql::metadata_mode::full);

    // Connect to the server using the first endpoint returned by the resolver
    conn.connect(*endpoints.begin(), params);

    // The executed commands may generate a lot of output, so we're going to
    // use multi-function operations (i.e. start_execution) to read it in batches.
    boost::mysql::execution_state st;
    conn.start_execution(script_contents, st);

    // The main read loop. Each executed command will yield a resultset.
    // st.comoplete() returns true once all resultsets have been read.
    for (std::size_t resultset_number = 0; !st.complete(); ++resultset_number)
    {
        // Advance to next resultset, if required
        if (st.should_read_head())
        {
            conn.read_resultset_head(st);
        }

        // Print the name of the fields
        std::cout << "Resultset number " << resultset_number << "\n";
        print_column_names(st.meta());

        // Read the rows and print them
        while (st.should_read_rows())
        {
            boost::mysql::rows_view batch = conn.read_some_rows(st);
            for (auto row : batch)
            {
                print_row(row);
            }
        }

        // Print OK packet data
        print_ok(st);
    }

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
        // field value that caused the error) and is encoded using to the connection's character set
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
