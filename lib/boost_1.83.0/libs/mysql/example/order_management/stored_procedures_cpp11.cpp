//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_stored_procedures

#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ssl/context.hpp>

#include <iostream>
#include <string>

/**
 * This example implements a very simple command-line order manager
 * for an online store, using stored procedures. You can find the procedure
 * definitions in example/db_setup_stored_procedures.sql. Be sure to run this file before the example.
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

    // Prints the details of an order to stdout. An order here is represented as a row
    static void print_order(mysql::row_view order)
    {
        std::cout << "Order: id=" << order.at(0) << ", status=" << order.at(1) << '\n';
    }

    // Prints the details of an order line item, again represented as a row
    static void print_line_item(mysql::row_view item)
    {
        std::cout << "  Line item: id=" << item.at(0) << ", quantity=" << item.at(1)
                  << ", unit_price=" << item.at(2).as_int64() / 100.0 << "$\n";
    }

    // Procedures that manipulate orders return two resultsets: one describing
    // the order and another with the line items the order has. Some of them
    // return only the order resultset. These functions print order details to stdout
    static void print_order_with_items(
        mysql::resultset_view order_resultset,
        mysql::resultset_view line_items_resultset
    )
    {
        // First resultset: order information. Always a single row
        print_order(order_resultset.rows().at(0));

        // Second resultset: all order line items
        mysql::rows_view line_items = line_items_resultset.rows();
        if (line_items.empty())
        {
            std::cout << "No line items\n";
        }
        else
        {
            for (mysql::row_view item : line_items)
            {
                print_line_item(item);
            }
        }
    }

    // get-products <search-term>: full text search of the products table
    void operator()(const get_products_args& args) const
    {
        // We need to pass user-supplied params to CALL, so we use a statement
        auto stmt = conn.prepare_statement("CALL get_products(?)");

        mysql::results result;
        conn.execute(stmt.bind(args.search), result);
        auto products = result.front();
        std::cout << "Your search returned the following products:\n";
        for (auto product : products.rows())
        {
            std::cout << "* ID: " << product.at(0) << '\n'
                      << "  Short name: " << product.at(1) << '\n'
                      << "  Description: " << product.at(2) << '\n'
                      << "  Price: " << product.at(3).as_int64() / 100.0 << "$" << std::endl;
        }
        std::cout << std::endl;
    }

    // create-order: creates a new order
    void operator()(const create_order_args&) const
    {
        // Since create_order doesn't have user-supplied params, we can use a text query
        mysql::results result;
        conn.execute("CALL create_order()", result);

        // Print the result to stdout. create_order() returns a resultset for
        // the newly created order, with only 1 row.
        std::cout << "Created order\n";
        print_order(result.at(0).rows().at(0));
    }

    // get-order <order-id>: retrieves order details
    void operator()(const get_order_args& args) const
    {
        // The order_id is supplied by the user, so we use a prepared statement
        auto stmt = conn.prepare_statement("CALL get_order(?)");

        // Execute the statement
        mysql::results result;
        conn.execute(stmt.bind(args.order_id), result);

        // Print the result to stdout. get_order() returns a resultset for
        // the retrieved order and another for the line items. If the order can't
        // be found, get_order() raises an error using SIGNAL, which will make
        // execute() fail with an exception.
        std::cout << "Retrieved order\n";
        print_order_with_items(result.at(0), result.at(1));
    }

    // get-orders: lists all orders
    void operator()(const get_orders_args&) const
    {
        // Since get_orders doesn't have user-supplied params, we can use a text query
        mysql::results result;
        conn.execute("CALL get_orders()", result);

        // Print results to stdout. get_orders() succeeds even if no order is found.
        // get_orders() only lists orders, not line items.
        mysql::rows_view orders = result.front().rows();
        if (orders.empty())
        {
            std::cout << "No orders found" << std::endl;
        }
        else
        {
            for (mysql::row_view order : result.front().rows())
            {
                print_order(order);
            }
        }
    }

    // add-line-item <order-id> <product-id> <quantity>: adds a line item to a given order
    void operator()(const add_line_item_args& args) const
    {
        // add_line_item has several user-supplied arguments, so we must use a statement.
        // The 4th argument is an OUT parameter. If we bind it by passing a ? marker,
        // we will get an extra resultset with just its value.
        auto stmt = conn.prepare_statement("CALL add_line_item(?, ?, ?, ?)");

        // We still have to pass a value to the 4th argument, even if it's an OUT parameter.
        // The value will be ignored, so we can pass nullptr.
        mysql::results result;
        conn.execute(stmt.bind(args.order_id, args.product_id, args.quantity, nullptr), result);

        // We can use results::out_params() to access the extra resultset containing
        // the OUT parameter
        auto new_line_item_id = result.out_params().at(0).as_int64();

        // Print the results to stdout
        std::cout << "Created line item: id=" << new_line_item_id << "\n";
        print_order_with_items(result.at(0), result.at(1));
    }

    // remove-line-item <line-item-id>: removes an item from an order
    void operator()(const remove_line_item_args& args) const
    {
        // remove_line_item has user-supplied parameters, so we use a statement
        auto stmt = conn.prepare_statement("CALL remove_line_item(?)");

        // Run the procedure
        mysql::results result;
        conn.execute(stmt.bind(args.line_item_id), result);

        // Print results to stdout
        std::cout << "Removed line item from order\n";
        print_order_with_items(result.at(0), result.at(1));
    }

    // checkout-order <order-id>: marks an order as ready for checkout
    void operator()(const checkout_order_args& args) const
    {
        // checkout_order has user-supplied parameters, so we use a statement.
        // The 2nd parameter represents the total order amount and is an OUT parameter.
        auto stmt = conn.prepare_statement("CALL checkout_order(?, ?)");

        // Execute the statement
        mysql::results result;
        conn.execute(stmt.bind(args.order_id, nullptr), result);

        // We can use results::out_params() to access the extra resultset containing
        // the OUT parameter
        auto total_amount = result.out_params().at(0).as_int64();

        // Print the results to stdout
        std::cout << "Checked out order. The total amount to pay is: " << total_amount / 100.0 << "$\n";
        print_order_with_items(result.at(0), result.at(1));
    }

    // complete-order <order-id>: marks an order as completed
    void operator()(const complete_order_args& args) const
    {
        // complete_order has user-supplied parameters, so we use a statement.
        auto stmt = conn.prepare_statement("CALL complete_order(?)");

        // Execute the statement
        mysql::results result;
        conn.execute(stmt.bind(args.order_id), result);

        // Print the results to stdout
        std::cout << "Completed order\n";
        print_order_with_items(result.at(0), result.at(1));
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
        // If a store procedure fails (e.g. because a SIGNAL statement was executed), an error
        // like this will be raised.
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