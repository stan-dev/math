//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_EXAMPLE_ORDER_MANAGEMENT_PARSE_CMDLINE_HPP
#define BOOST_MYSQL_EXAMPLE_ORDER_MANAGEMENT_PARSE_CMDLINE_HPP

#include <boost/mysql/string_view.hpp>

#include <boost/variant2/variant.hpp>

#include <cstdint>
#include <iostream>
#include <string>

/**
 * Our command line tool implements several sub-commands. Each sub-command
 * has a set of arguments. We define a struct for each sub-command.
 */

struct get_products_args
{
    std::string search;
};

struct create_order_args
{
};

struct get_order_args
{
    std::int64_t order_id;
};

struct get_orders_args
{
};

struct add_line_item_args
{
    std::int64_t order_id;
    std::int64_t product_id;
    std::int64_t quantity;
};

struct remove_line_item_args
{
    std::int64_t line_item_id;
};

struct checkout_order_args
{
    std::int64_t order_id;
};

struct complete_order_args
{
    std::int64_t order_id;
};

// A variant type that can represent arguments for any of the sub-commands
using any_command = boost::variant2::variant<
    get_products_args,
    get_order_args,
    get_orders_args,
    create_order_args,
    add_line_item_args,
    remove_line_item_args,
    checkout_order_args,
    complete_order_args>;

// In-memory representation of the command-line arguments once parsed.
struct cmdline_args
{
    const char* username;
    const char* password;
    const char* host;
    any_command cmd;
};

// Call on error to print usage and exit
[[noreturn]] inline void usage(boost::mysql::string_view program_name)
{
    std::cerr << "Usage: " << program_name << " <username> <password> <server-hostname> <command> args...\n"
              << "Available commands:\n"
                 "    get-products <search-term>\n"
                 "    create-order\n"
                 "    get-order <order-id>\n"
                 "    get-orders\n"
                 "    add-line-item <order-id> <product-id> <quantity>\n"
                 "    remove-line-item <line-item-id>\n"
                 "    checkout-order <order-id>\n"
                 "    complete-order <order-id>"
              << std::endl;
    exit(1);
}

// Helper function to parse a sub-command
inline any_command parse_subcommand(
    boost::mysql::string_view program_name,
    boost::mysql::string_view cmd_name,
    int argc_rest,
    char** argv_rest
)
{
    if (cmd_name == "get-products")
    {
        if (argc_rest != 1)
        {
            usage(program_name);
        }
        return get_products_args{argv_rest[0]};
    }
    else if (cmd_name == "create-order")
    {
        if (argc_rest != 0)
        {
            usage(program_name);
        }
        return create_order_args{};
    }
    else if (cmd_name == "get-order")
    {
        if (argc_rest != 1)
        {
            usage(program_name);
        }
        return get_order_args{std::stoi(argv_rest[0])};
    }
    else if (cmd_name == "get-orders")
    {
        if (argc_rest != 0)
        {
            usage(program_name);
        }
        return get_orders_args{};
    }
    else if (cmd_name == "add-line-item")
    {
        if (argc_rest != 3)
        {
            usage(program_name);
        }
        return add_line_item_args{
            std::stoi(argv_rest[0]),
            std::stoi(argv_rest[1]),
            std::stoi(argv_rest[2]),
        };
    }
    else if (cmd_name == "remove-line-item")
    {
        if (argc_rest != 1)
        {
            usage(program_name);
        }
        return remove_line_item_args{
            std::stoi(argv_rest[0]),
        };
    }
    else if (cmd_name == "checkout-order")
    {
        if (argc_rest != 1)
        {
            usage(program_name);
        }
        return checkout_order_args{std::stoi(argv_rest[0])};
    }
    else if (cmd_name == "complete-order")
    {
        if (argc_rest != 1)
        {
            usage(program_name);
        }
        return complete_order_args{std::stoi(argv_rest[0])};
    }
    else
    {
        usage(program_name);
    }
}

// Parses the entire command line
inline cmdline_args parse_cmdline_args(int argc, char** argv)
{
    if (argc < 5)
    {
        usage(argv[0]);
    }
    return cmdline_args{
        argv[1],
        argv[2],
        argv[3],
        parse_subcommand(argv[0], argv[4], argc - 5, argv + 5),
    };
}

#endif
