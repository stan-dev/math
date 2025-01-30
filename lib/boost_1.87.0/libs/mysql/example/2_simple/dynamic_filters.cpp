//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/awaitable.hpp>
#ifdef BOOST_ASIO_HAS_CO_AWAIT

//[example_dynamic_filters

/**
 * This example implements a dynamic filter using client-side SQL.
 * If you're implementing a filter with many options that can be
 * conditionally enabled, this pattern may be useful for you.
 *
 * This example uses C++20 coroutines. If you need, you can backport
 * it to C++11 by using callbacks, asio::yield_context
 * or sync functions instead of coroutines.
 *
 * This example uses the 'boost_mysql_examples' database, which you
 * can get by running db_setup.sql.
 */

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/character_set.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/format_sql.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/sequence.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/io_context.hpp>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace mysql = boost::mysql;
namespace asio = boost::asio;

// Prints an employee row to stdout
void print_employee(mysql::row_view employee)
{
    std::cout << "id: " << employee.at(0)                             // field 0: id
              << ", first_name: " << std::setw(16) << employee.at(1)  // field 1: first_name
              << ", last_name: " << std::setw(16) << employee.at(2)   // field 2: last_name
              << ", company_id: " << employee.at(3)                   // field 3: company_id
              << ", salary: " << employee.at(4) << '\n';              // field 4: salary
}

// An operator to use in a filter
enum class op_type
{
    lt,   // <
    lte,  // <=
    eq,   // =
    gt,   // >
    gte,  // >=
};

// Returns the SQL operator for the given op_type
std::string_view op_type_to_sql(op_type value)
{
    switch (value)
    {
    case op_type::lt: return "<";
    case op_type::lte: return "<=";
    case op_type::eq: return "=";
    case op_type::gte: return ">=";
    case op_type::gt: return ">";
    default: assert(false); return "=";
    }
}

// An individual filter to apply.
// For example, filter{"salary", op_type::gt, field_view(20000)} should generate a
// `salary` > 20000 condition
struct filter
{
    std::string_view field_name;    // The database column name
    op_type op;                     // The operator to apply
    mysql::field_view field_value;  // The value to check. field_view can hold any MySQL type
};

// Command line arguments
struct cmdline_args
{
    // MySQL username to use during authentication.
    std::string_view username;

    // MySQL password to use during authentication.
    std::string_view password;

    // Hostname where the MySQL server is listening.
    std::string_view server_hostname;

    // The filters to apply
    std::vector<filter> filts;

    // If order_by.has_value(), order employees using the given field
    std::optional<std::string_view> order_by;
};

// Parses the command line
static cmdline_args parse_cmdline_args(int argc, char** argv)
{
    // Available options
    constexpr std::string_view company_id_prefix = "--company-id=";
    constexpr std::string_view first_name_prefix = "--first-name=";
    constexpr std::string_view last_name_prefix = "--last-name=";
    constexpr std::string_view min_salary_prefix = "--min-salary=";
    constexpr std::string_view order_by_prefix = "--order-by=";

    // Helper function to print the usage message and exit
    auto print_usage_and_exit = [argv]() {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> [filters]\n";
        exit(1);
    };

    // Check number of arguments
    if (argc <= 4)
        print_usage_and_exit();

    // Parse the required arguments
    cmdline_args res;
    res.username = argv[1];
    res.password = argv[2];
    res.server_hostname = argv[3];

    // Parse the filters
    for (int i = 4; i < argc; ++i)
    {
        std::string_view arg = argv[i];

        // Attempt to match the argument against each prefix
        if (arg.starts_with(company_id_prefix))
        {
            auto value = arg.substr(company_id_prefix.size());
            res.filts.push_back({"company_id", op_type::eq, mysql::field_view(value)});
        }
        else if (arg.starts_with(first_name_prefix))
        {
            auto value = arg.substr(first_name_prefix.size());
            res.filts.push_back({"first_name", op_type::eq, mysql::field_view(value)});
        }
        else if (arg.starts_with(last_name_prefix))
        {
            auto value = arg.substr(last_name_prefix.size());
            res.filts.push_back({"last_name", op_type::eq, mysql::field_view(value)});
        }
        else if (arg.starts_with(min_salary_prefix))
        {
            auto value = std::stod(std::string(arg.substr(min_salary_prefix.size())));
            res.filts.push_back({"salary", op_type::gte, mysql::field_view(value)});
        }
        else if (arg.starts_with(order_by_prefix))
        {
            auto field_name = arg.substr(order_by_prefix.size());

            // For security, validate the passed field against a set of whitelisted fields
            if (field_name != "id" && field_name != "first_name" && field_name != "last_name" &&
                field_name != "salary")
            {
                std::cerr << "Order-by: invalid field " << field_name << std::endl;
                print_usage_and_exit();
            }
            res.order_by = field_name;
        }
        else
        {
            std::cerr << "Unrecognized option: " << arg << std::endl;
            print_usage_and_exit();
        }
    }

    // We should have at least one filter
    if (res.filts.empty())
    {
        std::cerr << "At least one filter should be specified" << std::endl;
        print_usage_and_exit();
    }

    return res;
}

// Composes a SELECT query to retrieve employees according to the passed filters.
// We allow an optional ORDER BY clause that must be added dynamically,
// so we can't express our query as a single format string.
// This function uses format_sql_to to build a query string incrementally.
// format_sql_to requires us to pass a format_options value, containing configuration
// options like the current character set. Use any_connection::format_opts to obtain it.
// If your use case allows you to express your query as a single format string, use with_params, instead.
std::string compose_get_employees_query(
    mysql::format_options opts,
    const std::vector<filter>& filts,
    std::optional<std::string_view> order_by
)
{
    // A format context allows composing queries incrementally.
    // This is required because we need to add the ORDER BY clause conditionally
    mysql::format_context ctx(opts);

    // Adds an individual filter to the context. Used by sequence()
    auto filter_format_fn = [](filter item, mysql::format_context_base& elm_ctx) {
        // {:i} formats a string as a SQL identifier. {:r} outputs raw SQL.
        // filter{"key", op_type::eq, field_view(42)} would get formatted as "`key` = 42"
        mysql::format_sql_to(
            elm_ctx,
            "{:i} {:r} {}",
            item.field_name,
            op_type_to_sql(item.op),
            item.field_value
        );
    };

    // Add the query with the filters to ctx.
    // sequence() will invoke filter_format_fn for each element in filts,
    // using the string " AND " as glue, to separate filters
    // By default, sequence copies its input range, but we don't need this here,
    // so we disable the copy by calling ref()
    mysql::format_sql_to(
        ctx,
        "SELECT id, first_name, last_name, company_id, salary FROM employee WHERE {}",
        mysql::sequence(std::ref(filts), filter_format_fn, " AND ")
    );

    // Add the order by
    if (order_by)
    {
        // identifier formats a string as a SQL identifier, instead of a string literal.
        // For instance, this may generate "ORDER BY `first_name`"
        mysql::format_sql_to(ctx, " ORDER BY {:i}", *order_by);
    }

    // Get our generated query. get() returns a system::result<std::string>, which
    // will contain errors if any of the args couldn't be formatted. This can happen
    // if you pass string values containing invalid UTF-8.
    // value() will throw an exception if that's the case.
    return std::move(ctx).get().value();
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
    params.password = args.password;
    params.database = "boost_mysql_examples";

    // Connect to the server
    co_await conn.async_connect(params);

    // Compose the query. format_opts() returns a system::result<format_options>,
    // containing the options required by format_context. format_opts() may return
    // an error if the connection doesn't know which character set is using -
    // use async_set_character_set if this happens.
    std::string query = compose_get_employees_query(conn.format_opts().value(), args.filts, args.order_by);

    // Execute the query as usual. Note that the query was generated
    // client-side. Appropriately using format_sql_to makes this approach secure.
    // with_params uses this same technique under the hood.
    // Casting to string_view saves a copy in async_execute
    mysql::results result;
    co_await conn.async_execute(std::string_view(query), result);

    // Print the employees
    for (mysql::row_view employee : result.rows())
    {
        print_employee(employee);
    }

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
        // You will only get this type of exceptions if you use with_diagnostics.
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
