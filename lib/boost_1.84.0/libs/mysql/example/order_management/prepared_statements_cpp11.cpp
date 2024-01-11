//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_prepared_statements_cpp11

/**
 * This example implements a very simple command-line order manager
 * for an online store, using prepared statements. You can find the table
 * definitions in example/order_management/db_setup.sql. Be sure to run this file before the example.
 * This example assumes you are connecting to a localhost MySQL server.
 *
 * The order system is intentionally very simple, and has the following tables:
 *  - products: the list of items our store sells, with price and description.
 *  - orders: the main object. Orders have a status field that can be draft, pending_payment or complete.
 *  - order_items: an order may have 0 to n line items. Each item refers to a single product.
 *
 * Orders are created empty, in a draft state. Line items can be added or removed.
 * Orders are then checked out, which transitions them to pending_payment.
 * After that, payment would happen through an external system. Once completed, an
 * order is confirmed, transitioning it to the complete status.
 * In the real world, flow would be much more complex, but this is enough for an example.
 *
 * We'll be using the untyped interface to retrieve results from MySQL.
 * This makes use of the results, rows_view, row_view and field_view classes.
 * If you prefer typing your rows statically, you may prefer using the "typed interface",
 * which uses static_results instead.
 */

#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ssl/context.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

// This header contains boilerplate code to parse the command line
// arguments into structs. Parsing the command line yields a cmdline_args,
// an alias for a boost::variant2::variant holding the command line
// arguments for any of the subcommands. We will use it via visit().
#include "parse_cmdline.hpp"

namespace mysql = boost::mysql;

namespace {

// This visitor executes a sub-command and prints the results to stdout.
struct visitor
{
    mysql::tcp_ssl_connection& conn;

    // Retrieves an order with its items. If the order exists, at least one record is returned.
    // If the order has line items, a record per item is returned. If the order has no items,
    // a single record is returned, and it will have its item_xxx fields set to NULL.
    mysql::results get_order_with_items(std::int64_t order_id) const
    {
        // Prepare a statement, since the order_id is provided by the user
        mysql::statement stmt = conn.prepare_statement(R"%(
            SELECT
                ord.id AS order_id,
                ord.status AS order_status,
                item.id AS item_id,
                item.quantity AS item_quantity,
                prod.price AS item_unit_price
            FROM orders ord
            LEFT JOIN order_items item ON ord.id = item.order_id
            LEFT JOIN products prod ON item.product_id = prod.id
            WHERE ord.id = ?
        )%");

        // Execute it
        mysql::results result;
        conn.execute(stmt.bind(order_id), result);

        return result;
    }

    // Prints an order with its line items to stdout. Each row has the format
    // described in get_order_with_items
    static void print_order_with_items(mysql::rows_view ord_items)
    {
        // Print the order
        std::cout << "Order: id=" << ord_items.at(0).at(0) << ", status=" << ord_items.at(0).at(1) << '\n';

        // Print the items (3rd to 5th fields). These will be NULL if the order
        // contains no items.
        if (ord_items.at(0).at(2).is_null())
        {
            std::cout << "No line items\n";
        }
        else
        {
            for (mysql::row_view item : ord_items)
            {
                std::cout << "  Line item: id=" << item.at(2) << ", quantity=" << item.at(3)
                          << ", unit_price=" << item.at(4).as_int64() / 100.0 << "$\n";
            }
        }
    }

    // get-products <search-term>: full text search of the products table.
    // use this command to search the store for available products
    void operator()(const get_products_args& args) const
    {
        // Our SQL contains a user-supplied paremeter (the search term),
        // so we will be using a prepared statement
        mysql::statement stmt = conn.prepare_statement(R"%(
            SELECT id, short_name, descr, price
            FROM products
            WHERE MATCH(short_name, descr) AGAINST(?)
            LIMIT 5
        )%");

        // Execute it
        mysql::results products;
        conn.execute(stmt.bind(args.search), products);

        // Print the results to stdout
        std::cout << "Your search returned the following products:\n";
        for (mysql::row_view prod : products.rows())
        {
            std::cout << "* ID: " << prod.at(0) << '\n'
                      << "  Short name: " << prod.at(1) << '\n'
                      << "  Description: " << prod.at(2) << '\n'
                      << "  Price: " << prod.at(3).as_int64() / 100.0 << "$" << std::endl;
        }
        std::cout << std::endl;
    }

    // create-order: creates a new order. Orders are always created empty. This command
    // requires no arguments
    void operator()(const create_order_args&) const
    {
        // Our SQL doesn't have parameters, so we can use a text query.
        mysql::results result;
        conn.execute("INSERT INTO orders VALUES ()", result);

        // Print the results to stdout. results::last_insert_id() returns the ID of
        // the newly inserted order.
        std::cout << "Order: id=" << result.last_insert_id() << ", status=draft" << std::endl;
    }

    // get-order <order-id>: retrieves order details
    void operator()(const get_order_args& args) const
    {
        // Retrieve the order with its items
        mysql::results result = get_order_with_items(args.order_id);

        // If we didn't find any order, issue an error
        if (result.rows().empty())
        {
            throw std::runtime_error("Can't find order with id=" + std::to_string(args.order_id));
        }

        // Print the order to stdout
        std::cout << "Retrieved order\n";
        print_order_with_items(result.rows());
    }

    // get-orders: lists all orders. Orders are listed without their line items.
    void operator()(const get_orders_args&) const
    {
        // Since this query doesn't have parameters, we don't need a prepared statement,
        // and we can use a text query instead.
        mysql::results result;
        conn.execute("SELECT id, `status` FROM orders", result);

        // Print the results to stdout
        if (result.rows().empty())
        {
            std::cout << "No orders found" << std::endl;
        }
        else
        {
            for (mysql::row_view order : result.rows())
            {
                std::cout << "Order: id=" << order.at(0) << ", status=" << order.at(1) << '\n';
            }
        }
    }

    // add-line-item <order-id> <product-id> <quantity>: adds a line item to a given order
    void operator()(const add_line_item_args& args) const
    {
        // We will need to run several statements atomically, so we start a transaction.
        mysql::results result;
        conn.execute("START TRANSACTION", result);

        // To add a line item, we require the order to be in a draft status. Get the order to check this fact.
        mysql::statement stmt = conn.prepare_statement("SELECT `status` FROM orders WHERE id = ?");
        conn.execute(stmt.bind(args.order_id), result);
        if (result.rows().empty())
        {
            // There is no such order
            throw std::runtime_error("Order with id=" + std::to_string(args.order_id) + " not found");
        }
        else if (result.rows()[0].at(0).as_string() != "draft")
        {
            // The order is no longer editable
            throw std::runtime_error("Order with id=" + std::to_string(args.order_id) + " is not editable");
        }

        // Insert the new line item. If the given product does not exist, the INSERT will fail
        // because of product_id's FOREIGN KEY constraint.
        stmt = conn.prepare_statement(
            "INSERT INTO order_items (order_id, product_id, quantity) VALUES (?, ?, ?)"
        );
        conn.execute(stmt.bind(args.order_id, args.product_id, args.quantity), result);

        // We can use results::last_insert_id to get the ID of the new line item.
        auto new_line_item_id = result.last_insert_id();

        // Retrieve the full order details
        mysql::results order_results = get_order_with_items(args.order_id);

        // We're done - commit the transaction
        conn.execute("COMMIT", result);

        // Print the results to stdout
        std::cout << "Created line item: id=" << new_line_item_id << "\n";
        print_order_with_items(order_results.rows());
    }

    // remove-line-item <line-item-id>: removes an item from an order
    void operator()(const remove_line_item_args& args) const
    {
        // We will need to run several statements atomically, so we start a transaction.
        mysql::results result;
        conn.execute("START TRANSACTION", result);

        // To remove a line item, we require the order to be in a draft status. Get the order to check it.
        auto stmt = conn.prepare_statement(R"%(
            SELECT orders.id, orders.`status`
            FROM orders
            JOIN order_items items ON (orders.id = items.order_id)
            WHERE items.id = ?
        )%");
        conn.execute(stmt.bind(args.line_item_id), result);
        if (result.rows().empty())
        {
            // The query hasn't matched any row - the supplied line item ID is not valid
            throw std::runtime_error(
                "The order item with id=" + std::to_string(args.line_item_id) + " does not exist"
            );
        }
        mysql::row_view order = result.rows()[0];
        if (order.at(1).as_string() != "draft")
        {
            // The order is no longer editable
            throw std::runtime_error("The order is not in an editable state");
        }

        // Remove the line item
        stmt = conn.prepare_statement("DELETE FROM order_items WHERE id = ?");
        conn.execute(stmt.bind(args.line_item_id), result);

        // Retrieve the full order details
        mysql::results order_results = get_order_with_items(order.at(0).as_int64());

        // We're done - commit the transaction
        conn.execute("COMMIT", result);

        // Print results to stdout
        std::cout << "Removed line item from order\n";
        print_order_with_items(order_results.rows());
    }

    // checkout-order <order-id>: marks an order as ready for checkout
    void operator()(const checkout_order_args& args) const
    {
        // We will need to run several statements atomically, so we start a transaction.
        mysql::results result;
        conn.execute("START TRANSACTION", result);

        // To checkout an order, we require it to be in a draft status. Check this fact.
        mysql::statement stmt = conn.prepare_statement("SELECT `status` FROM orders WHERE id = ?");
        conn.execute(stmt.bind(args.order_id), result);
        if (result.rows().empty())
        {
            // No order matched
            throw std::runtime_error("Order with id=" + std::to_string(args.order_id) + " not found");
        }
        else if (result.rows()[0].at(0).as_string() != "draft")
        {
            // The order is no longer editable
            throw std::runtime_error(
                "Order with id=" + std::to_string(args.order_id) + " cannot be checked out"
            );
        }

        // Update the order status
        stmt = conn.prepare_statement("UPDATE orders SET `status` = 'pending_payment' WHERE id = ?");
        conn.execute(stmt.bind(args.order_id), result);

        // Calculate the total amount to pay. SUM() returns a DECIMAL, which has a bigger
        // range than integers. DECIMAL is represented in C++ as a string. We use CAST to obtain
        // an uint64_t. If the CAST overflows, the max value for uint64_t will be returned.
        // We will be limiting our orders to USD 1bn, so overflow will be detected.
        stmt = conn.prepare_statement(R"%(
            SELECT CAST(
                IFNULL(SUM(prod.price * item.quantity), 0)
                AS UNSIGNED
            )
            FROM order_items item
            JOIN products prod ON item.product_id = prod.id
            WHERE item.order_id = ?;
        )%");
        conn.execute(stmt.bind(args.order_id), result);
        std::uint64_t total_amount = result.rows().at(0).at(0).as_uint64();

        // Verify that the total amount meets our criteria
        if (total_amount == 0)
        {
            throw std::runtime_error("The order doesn't have any line item");
        }
        else if (total_amount > 1000 * 1000 * 100)
        {
            throw std::runtime_error("Order amount of " + std::to_string(total_amount) + " exceeds limit");
        }

        // Retrieve the full order details
        mysql::results order_results = get_order_with_items(args.order_id);

        // We're done - commit the transaction
        conn.execute("COMMIT", result);

        // Print the results to stdout
        std::cout << "Checked out order. The total amount to pay is: " << total_amount / 100.0 << "$\n";
        print_order_with_items(order_results.rows());
    }

    // complete-order <order-id>: marks an order as completed
    void operator()(const complete_order_args& args) const
    {
        // We will need to run several statements atomically, so we start a transaction.
        mysql::results result;
        conn.execute("START TRANSACTION", result);

        // To complete an order, we require it to be in a pending_payment status. Check this fact.
        auto stmt = conn.prepare_statement("SELECT `status` FROM orders WHERE id = ?");
        conn.execute(stmt.bind(args.order_id), result);
        if (result.rows().empty())
        {
            // Order not found
            throw std::runtime_error("Order with id=" + std::to_string(args.order_id) + " not found");
        }
        else if (result.rows()[0].at(0).as_string() != "pending_payment")
        {
            throw std::runtime_error(
                "Order with id=" + std::to_string(args.order_id) + " is not in pending_payment status"
            );
        }

        // Update status
        stmt = conn.prepare_statement("UPDATE orders SET `status` = 'complete' WHERE id = ?");
        conn.execute(stmt.bind(args.order_id), result);

        // Retrieve the full order details
        mysql::results order_results = get_order_with_items(args.order_id);

        // We're done - commit the transaction
        conn.execute("COMMIT", result);

        // Print the results to stdout
        std::cout << "Completed order\n";
        print_order_with_items(order_results.rows());
    }
};

void main_impl(int argc, char** argv)
{
    // Parse command line arguments
    auto args = parse_cmdline_args(argc, argv);

    // I/O context and connection. We use SSL because MySQL 8+ default settings require it.
    boost::asio::io_context ctx;
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    mysql::tcp_ssl_connection conn(ctx, ssl_ctx);

    // Resolver for hostname resolution
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());

    // Connection params
    mysql::handshake_params params(
        args.username,                  // username
        args.password,                  // password
        "boost_mysql_order_management"  // database to use
    );

    // Hostname resolution
    auto endpoints = resolver.resolve(args.host, mysql::default_port_string);

    // TCP and MySQL level connect
    conn.connect(*endpoints.begin(), params);

    // Execute the command
    boost::variant2::visit(visitor{conn}, args.cmd);

    // Close the connection
    conn.close();
}

}  // namespace

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
