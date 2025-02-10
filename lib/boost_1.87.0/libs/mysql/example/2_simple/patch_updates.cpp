//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/awaitable.hpp>
#ifdef BOOST_ASIO_HAS_CO_AWAIT

//[example_patch_updates

/**
 * This example demonstrates how to implement dynamic updates
 * with PATCH-like semantics using client-side SQL formatting.
 *
 * The program updates an employee by ID, modifying fields
 * as provided by command-line arguments, and leaving all other
 * fields unmodified.
 *
 * This example uses C++20 coroutines. If you need, you can backport
 * it to C++14 (required by Boost.Describe) by using callbacks, asio::yield_context
 * or sync functions instead of coroutines.
 *
 * This example uses the 'boost_mysql_examples' database, which you
 * can get by running db_setup.sql.
 */

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/sequence.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>

#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>

namespace mysql = boost::mysql;
namespace asio = boost::asio;

/**
 * Represents a single update as a name, value pair.
 * The idea is to use command-line arguments to compose
 * a std::vector<update_field> with the fields to be updated,
 * and use mysql::sequence() to join these with commas
 */
struct update_field
{
    // The field name to set (i.e. the column name)
    std::string_view field_name;

    // The value to set the field to. Recall that field_view is
    // a variant-like type that can hold all types that MySQL supports.
    mysql::field_view field_value;
};

// Contains the parsed command-line arguments
struct cmdline_args
{
    // MySQL username to use during authentication.
    std::string_view username;

    // MySQL password to use during authentication.
    std::string_view password;

    // Hostname where the MySQL server is listening.
    std::string_view server_hostname;

    // The ID of the employee we want to update.
    std::int64_t employee_id{};

    // A list of name, value pairs containing the employee fields to update.
    std::vector<update_field> updates;
};

// Parses the command line arguments, calling exit on failure.
static cmdline_args parse_cmdline_args(int argc, char** argv)
{
    // Available options
    constexpr std::string_view company_id_prefix = "--company-id=";
    constexpr std::string_view first_name_prefix = "--first-name=";
    constexpr std::string_view last_name_prefix = "--last-name=";
    constexpr std::string_view salary_prefix = "--salary=";

    // Helper function to print the usage message and exit
    auto print_usage_and_exit = [argv]() {
        std::cerr << "Usage: " << argv[0]
                  << " <username> <password> <server-hostname> employee_id [updates]\n";
        exit(1);
    };

    // Check number of arguments
    if (argc <= 5)
        print_usage_and_exit();

    // Parse the required arguments
    cmdline_args res;
    res.username = argv[1];
    res.password = argv[2];
    res.server_hostname = argv[3];
    res.employee_id = std::stoll(argv[4]);

    // Parse the requested updates
    for (int i = 5; i < argc; ++i)
    {
        // Get the argument
        std::string_view arg = argv[i];

        // Attempt to match it with the options we have
        if (arg.starts_with(company_id_prefix))
        {
            std::string_view new_value = arg.substr(company_id_prefix.size());
            res.updates.push_back(update_field{"company_id", mysql::field_view(new_value)});
        }
        else if (arg.starts_with(first_name_prefix))
        {
            std::string_view new_value = arg.substr(first_name_prefix.size());
            res.updates.push_back(update_field{"first_name", mysql::field_view(new_value)});
        }
        else if (arg.starts_with(last_name_prefix))
        {
            std::string_view new_value = arg.substr(last_name_prefix.size());
            res.updates.push_back(update_field{"last_name", mysql::field_view(new_value)});
        }
        else if (arg.starts_with(salary_prefix))
        {
            double new_value = std::stod(std::string(arg.substr(salary_prefix.size())));
            res.updates.push_back(update_field{"salary", mysql::field_view(new_value)});
        }
        else
        {
            std::cerr << "Unrecognized option: " << arg << std::endl;
            print_usage_and_exit();
        }
    }

    // There should be one update, at least
    if (res.updates.empty())
    {
        std::cerr << "There should be one update, at least\n";
        print_usage_and_exit();
    }

    return res;
}

// The main coroutine
asio::awaitable<void> coro_main(const cmdline_args& args)
{
    // Create a connection.
    // Will use the same executor as the coroutine.
    mysql::any_connection conn(co_await asio::this_coro::executor);

    // The hostname, username, password and database to use
    mysql::connect_params params;
    params.server_address.emplace_host_and_port(std::string(args.server_hostname));
    params.username = args.username;
    params.password = std::string(args.password);
    params.database = "boost_mysql_examples";
    params.multi_queries = true;

    // Connect to the server
    co_await conn.async_connect(params);

    // Formats an individual update. Used by sequence().
    // For update_field{"first_name", "John"}, it generates the string
    // "`first_name` = 'John'"
    // Format contexts can build a query string incrementally, and are used by sequence() internally
    auto update_format_fn = [](update_field upd, mysql::format_context_base& ctx) {
        mysql::format_sql_to(ctx, "{:i} = {}", upd.field_name, upd.field_value);
    };

    // Compose and execute the query. with_params will expand placeholders
    // before sending the query to the server.
    // We use sequence() to output the update list separated by commas.
    // We want to update the employee and then retrieve it. MySQL doesn't support
    // the UPDATE ... RETURNING statement to update and retrieve data atomically,
    // so we will use a transaction to guarantee consistency.
    // Instead of running every statement separately, we activated params.multi_queries,
    // which allows semicolon-separated statements.
    // As in std::format, we can use explicit indices like {0} and {1} to reference arguments.
    // By default, sequence copies its input range, but we don't need this here,
    // so we disable the copy by calling ref()
    mysql::results result;
    co_await conn.async_execute(
        mysql::with_params(
            "START TRANSACTION; "
            "UPDATE employee SET {0} WHERE id = {1}; "
            "SELECT first_name, last_name, salary, company_id FROM employee WHERE id = {1}; "
            "COMMIT",
            mysql::sequence(std::ref(args.updates), update_format_fn),
            args.employee_id
        ),
        result
    );

    // We ran 4 queries, so the results object will hold 4 resultsets.
    // Get the rows retrieved by the SELECT (the 3rd one).
    auto rws = result.at(2).rows();

    // If there are no rows, the given employee does not exist.
    if (rws.empty())
    {
        std::cerr << "employee_id=" << args.employee_id << " not found" << std::endl;
        exit(1);
    }

    // Print the updated employee.
    const auto employee = rws.at(0);
    std::cout << "Updated employee with id=" << args.employee_id << ":\n"
              << "  first_name: " << employee.at(0) << "\n  last_name: " << employee.at(1)
              << "\n  salary: " << employee.at(2) << "\n  company_id: " << employee.at(3) << std::endl;

    // Notify the MySQL server we want to quit, then close the underlying connection.
    co_await conn.async_close();
}

void main_impl(int argc, char** argv)
{
    // Parse the command line
    cmdline_args args = parse_cmdline_args(argc, argv);

    // Create an I/O context, required by all I/O objects
    asio::io_context ctx;

    // Launch our coroutine
    asio::co_spawn(
        ctx,
        [&] { return coro_main(args); },
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

#else

#include <iostream>

int main()
{
    std::cout << "Sorry, your compiler doesn't have the required capabilities to run this example"
              << std::endl;
}

#endif
