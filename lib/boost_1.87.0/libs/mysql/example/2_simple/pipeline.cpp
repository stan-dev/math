//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/awaitable.hpp>
#ifdef BOOST_ASIO_HAS_CO_AWAIT

//[example_pipeline

/**
 * (EXPERIMENTAL)
 * This example demonstrates how to use the pipeline API to prepare,
 * execute and close statements in batch.
 * Pipelines are a experimental API.
 *
 * This example uses C++20 coroutines. If you need, you can backport
 * it to C++11 by using callbacks or asio::yield_context.
 * Timeouts can't be used with sync functions.
 *
 * This example uses the 'boost_mysql_examples' database, which you
 * can get by running db_setup.sql.
 */

#include <boost/mysql/any_address.hpp>
#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/pipeline.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>

#include <array>
#include <iostream>
#include <span>
#include <vector>

namespace asio = boost::asio;
namespace mysql = boost::mysql;

// Prepare several statements in batch.
// This is faster than preparing them one by one, as it saves round-trips to the server.
asio::awaitable<std::vector<mysql::statement>> batch_prepare(
    mysql::any_connection& conn,
    std::span<const std::string_view> statements
)
{
    // Construct a pipeline request describing the work to be performed.
    // There must be one prepare_statement_stage per statement to prepare
    mysql::pipeline_request req;
    for (auto stmt_sql : statements)
        req.add_prepare_statement(stmt_sql);

    // Run the pipeline.
    // stage_response is a variant-like type that can hold the response of any stage type.
    std::vector<mysql::stage_response> pipe_res;
    co_await conn.async_run_pipeline(req, pipe_res);

    // If we got here, all statements were prepared successfully.
    // pipe_res contains as many elements as statements.size(), holding statement objects
    // Extract them into a vector
    std::vector<mysql::statement> res;
    res.reserve(statements.size());
    for (const auto& stage_res : pipe_res)
        res.push_back(stage_res.get_statement());
    co_return res;
}

// The main coroutine
asio::awaitable<void> coro_main(
    std::string_view server_hostname,
    std::string_view username,
    std::string_view password,
    std::string_view company_id
)
{
    // Create a connection.
    // Will use the same executor as the coroutine.
    mysql::any_connection conn(co_await asio::this_coro::executor);

    // The hostname, username, password and database to use
    mysql::connect_params params;
    params.server_address.emplace_host_and_port(std::string(server_hostname));
    params.username = username;
    params.password = password;
    params.database = "boost_mysql_examples";

    // Connect to server
    co_await conn.async_connect(params);

    // Prepare the statements using the batch prepare function that we previously defined
    const std::array<std::string_view, 2> stmt_sql{
        "INSERT INTO employee (company_id, first_name, last_name) VALUES (?, ?, ?)",
        "INSERT INTO audit_log (msg) VALUES (?)"
    };
    std::vector<mysql::statement> stmts = co_await batch_prepare(conn, stmt_sql);

    // Create a pipeline request to execute them.
    // Warning: do NOT include the COMMIT statement in this pipeline.
    // COMMIT must only be executed if all the previous statements succeeded.
    // In a pipeline, all stages get executed, regardless of the outcome of previous stages.
    // We say that COMMIT has a dependency on the result of previous stages.
    mysql::pipeline_request req;
    req.add_execute("START TRANSACTION")
        .add_execute(stmts.at(0), company_id, "Juan", "Lopez")
        .add_execute(stmts.at(0), company_id, "Pepito", "Rodriguez")
        .add_execute(stmts.at(0), company_id, "Someone", "Random")
        .add_execute(stmts.at(1), "Inserted 3 new emplyees");
    std::vector<mysql::stage_response> res;

    // Execute the pipeline
    co_await conn.async_run_pipeline(req, res);

    // If we got here, all stages executed successfully.
    // Since they were execution stages, the response contains a results object.
    // Get the IDs of the newly created employees
    auto id1 = res.at(1).as_results().last_insert_id();
    auto id2 = res.at(2).as_results().last_insert_id();
    auto id3 = res.at(3).as_results().last_insert_id();

    // We can now commit our transaction and close the statements.
    // Clear the request and populate it again
    req.clear();
    req.add_execute("COMMIT").add_close_statement(stmts.at(0)).add_close_statement(stmts.at(1));

    // Run it
    co_await conn.async_run_pipeline(req, res);

    // If we got here, our insertions got committed.
    std::cout << "Inserted employees: " << id1 << ", " << id2 << ", " << id3 << std::endl;

    // Notify the MySQL server we want to quit, then close the underlying connection.
    co_await conn.async_close();
}

void main_impl(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> <company-id>\n";
        exit(1);
    }

    // Create an I/O context, required by all I/O objects
    asio::io_context ctx;

    // Launch our coroutine
    asio::co_spawn(
        ctx,
        [=] { return coro_main(argv[3], argv[1], argv[2], argv[4]); },
        // If any exception is thrown in the coroutine body, rethrow it.
        [](std::exception_ptr ptr) {
            if (ptr)
            {
                std::rethrow_exception(ptr);
            }
        }
    );

    // Calling run will actually execute the coroutine until completion
    ctx.run();
}

int main(int argc, char** argv)
{
    try
    {
        main_impl(argc, argv);
    }
    catch (const mysql::error_with_diagnostics& err)
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

#else

#include <iostream>

int main()
{
    std::cout << "Sorry, your compiler doesn't have the required capabilities to run this example"
              << std::endl;
}

#endif
