//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_stored_procedures_cpp14

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
 * We'll be using the static interface to retrieve results from MySQL.
 * This makes use of the static_results<T1, T2...> class template.
 * To use it, we need to define a set of structs/tuples describing the shape
 * of our rows. Boost.MySQL will parse the received rows into these types.
 * The static interface requires C++14 to work.
 *
 * Row types may be plain structs or std::tuple's. If we use plain structs, we need
 * to use BOOST_DESCRIBE_STRUCT on them. This adds the structs the required reflection
 * data, so Boost.MySQL knows how to parse rows into them.
 */

#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/static_results.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/describe/class.hpp>
#include <boost/optional/optional.hpp>

#include <iostream>
#include <string>
#include <tuple>

// This header contains boilerplate code to parse the command line
// arguments into structs. Parsing the command line yields a cmdline_args,
// an alias for a boost::variant2::variant holding the command line
// arguments for any of the subcommands. We will use it via visit().
#include "parse_cmdline.hpp"

// Including any of the static interface headers brings this macro into
// scope if the static interface is supported.
#ifdef BOOST_MYSQL_CXX14

namespace mysql = boost::mysql;

namespace {

// An order retrieved by our system.
struct order
{
    // The unique database ID of the object.
    std::int64_t id;

    // The order status (draft, pending_payment, complete).
    std::string status;
};
BOOST_DESCRIBE_STRUCT(order, (), (id, status))

// A line item, associated to an order and to a product.
// Our queries don't retrieve the order or product ID, so
// we don't include them in this struct.
struct order_item
{
    // The unique database ID of the object.
    std::int64_t id;

    // The number of units of this product that the user wants to buy.
    std::int64_t quantity;

    // The product's unit price, in cents of USD.
    std::int64_t unit_price;
};
BOOST_DESCRIBE_STRUCT(order_item, (), (id, quantity, unit_price))

// A product, as listed in the store product catalog.
struct product
{
    // The unique database ID of the object.
    std::int64_t id;

    // A short name for the product. Can be used as a title.
    std::string short_name;

    // The product's description. This field can be NULL in the DB,
    // so we use boost::optional<T> for it. If you're using C++17 or higher,
    // you can use std::optional instead.
    boost::optional<std::string> descr;

    // The product's unit price, in cents of USD.
    std::int64_t price;
};
BOOST_DESCRIBE_STRUCT(product, (), (id, short_name, descr, price))

// An empty row type. This can be used to describe empty resultsets,
// like the ones returned by INSERT or CALL.
using empty = std::tuple<>;

// This visitor executes a sub-command and prints the results to stdout.
struct visitor
{
    mysql::tcp_ssl_connection& conn;

    // Prints the details of an order to stdout
    static void print_order(const order& ord)
    {
        std::cout << "Order: id=" << ord.id << ", status=" << ord.status << '\n';
    }

    // Prints the details of an order line item
    static void print_line_item(const order_item& item)
    {
        std::cout << "  Line item: id=" << item.id << ", quantity=" << item.quantity
                  << ", unit_price=" << item.unit_price / 100.0 << "$\n";
    }

    // Prints an order with its line items to stdout
    static void print_order_with_items(const order& ord, boost::span<const order_item> items)
    {
        print_order(ord);

        if (items.empty())
        {
            std::cout << "No line items\n";
        }
        else
        {
            for (const auto& item : items)
            {
                print_line_item(item);
            }
        }
    }

    // get-products <search-term>: full text search of the products table.
    // use this command to search the store for available products
    void operator()(const get_products_args& args) const
    {
        // We need to pass user-supplied params to CALL, so we use a statement
        auto stmt = conn.prepare_statement("CALL get_products(?)");

        // get_products returns two resultsets:
        //   1. A collection of products
        //   2. An empty resultset describing the effects of the CALL statement
        mysql::static_results<product, empty> products;
        conn.execute(stmt.bind(args.search), products);

        // Print the results to stdout. By default, rows() returns the rows for the 1st resultset.
        std::cout << "Your search returned the following products:\n";
        for (const product& prod : products.rows())
        {
            std::cout << "* ID: " << prod.id << '\n'
                      << "  Short name: " << prod.short_name << '\n'
                      << "  Description: " << (prod.descr ? *prod.descr : "") << '\n'
                      << "  Price: " << prod.price / 100.0 << "$" << std::endl;
        }
        std::cout << std::endl;
    }

    // create-order: creates a new order. Orders are always created empty. This command
    // requires no arguments
    void operator()(const create_order_args&) const
    {
        // Since create_order doesn't have user-supplied params, we can use a text query.
        // create_order returns two resultsets:
        //   1. The created order. This is always a single row.
        //   2. An empty resultset describing the effects of the CALL statement
        mysql::static_results<order, empty> result;
        conn.execute("CALL create_order()", result);

        // Print the result to stdout. create_order() returns a resultset for
        // the newly created order, with only 1 row.
        std::cout << "Created order\n";
        print_order(result.rows()[0]);
    }

    // get-order <order-id>: retrieves order details
    void operator()(const get_order_args& args) const
    {
        // The order_id is supplied by the user, so we use a prepared statement
        auto stmt = conn.prepare_statement("CALL get_order(?)");

        // get_order returns three resultsets:
        //   1. The retrieved order. This is always a single row.
        //   2. A collection of line items for this order.
        //   3. An empty resultset describing the effects of the CALL statement
        // If the order can't be found, get_order() raises an error using SIGNAL, which will make
        // execute() fail with an exception.
        mysql::static_results<order, order_item, empty> result;
        conn.execute(stmt.bind(args.order_id), result);

        // Print the result to stdout.
        // rows<i>() can be used to access rows for the i-th resultset.
        // rows() means rows<0>().
        std::cout << "Retrieved order\n";
        print_order_with_items(result.rows<0>()[0], result.rows<1>());
    }

    // get-orders: lists all orders
    void operator()(const get_orders_args&) const
    {
        // Since get_orders doesn't have user-supplied params, we can use a text query
        // get_orders returns two resultsets:
        //   1. A collection of orders.
        //   2. An empty resultset describing the effects of the CALL statement
        mysql::static_results<order, empty> result;
        conn.execute("CALL get_orders()", result);

        // Print results to stdout. get_orders() succeeds even if no order is found.
        // get_orders() only lists orders, not line items.
        if (result.rows().empty())
        {
            std::cout << "No orders found" << std::endl;
        }
        else
        {
            for (const order& ord : result.rows())
            {
                print_order(ord);
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

        // add_line_item returns four resultsets:
        //   1. The affected order. Always a single row.
        //   2. A collection of line items for the affected order.
        //   3. An OUT params resultset, containing the ID of the newly created line item. Single row.
        //      MySQL always marks OUT params as nullable.
        //   4. An empty resultset describing the effects of the CALL statement
        using out_params_t = std::tuple<boost::optional<std::int64_t>>;
        mysql::static_results<order, order_item, out_params_t, empty> result;

        // We still have to pass a value to the 4th argument, even if it's an OUT parameter.
        // The value will be ignored, so we can pass nullptr.
        conn.execute(stmt.bind(args.order_id, args.product_id, args.quantity, nullptr), result);

        // We can access the OUT param as we access any other resultset
        auto new_line_item_id = std::get<0>(result.rows<2>()[0]).value();

        // Print the results to stdout
        std::cout << "Created line item: id=" << new_line_item_id << "\n";
        print_order_with_items(result.rows<0>()[0], result.rows<1>());
    }

    // remove-line-item <line-item-id>: removes an item from an order
    void operator()(const remove_line_item_args& args) const
    {
        // remove_line_item has user-supplied parameters, so we use a statement
        auto stmt = conn.prepare_statement("CALL remove_line_item(?)");

        // remove_line_item returns three resultsets:
        //   1. The affected order. Always a single row.
        //   2. A collection of line items for the affected order.
        //   3. An empty resultset describing the effects of the CALL statement
        mysql::static_results<order, order_item, empty> result;
        conn.execute(stmt.bind(args.line_item_id), result);

        // Print results to stdout
        std::cout << "Removed line item from order\n";
        print_order_with_items(result.rows<0>()[0], result.rows<1>());
    }

    // checkout-order <order-id>: marks an order as ready for checkout
    void operator()(const checkout_order_args& args) const
    {
        // checkout_order has user-supplied parameters, so we use a statement.
        // The 2nd parameter represents the total order amount and is an OUT parameter.
        auto stmt = conn.prepare_statement("CALL checkout_order(?, ?)");

        // checkout_order returns four resultsets:
        //   1. The affected order. Always a single row.
        //   2. A collection of line items for the affected order.
        //   3. An OUT params resultset, containing the total amount to pay, in USD cents. Single row.
        //      MySQL always marks OUT params as nullable.
        //   4. An empty resultset describing the effects of the CALL statement
        using out_params_t = std::tuple<boost::optional<std::int64_t>>;
        mysql::static_results<order, order_item, out_params_t, empty> result;
        conn.execute(stmt.bind(args.order_id, nullptr), result);

        // We can access the OUT param as we access any other resultset
        auto total_amount = std::get<0>(result.rows<2>()[0]).value_or(0);

        // Print the results to stdout
        std::cout << "Checked out order. The total amount to pay is: " << total_amount / 100.0 << "$\n";
        print_order_with_items(result.rows<0>()[0], result.rows<1>());
    }

    // complete-order <order-id>: marks an order as completed
    void operator()(const complete_order_args& args) const
    {
        // complete_order has user-supplied parameters, so we use a statement.
        auto stmt = conn.prepare_statement("CALL complete_order(?)");

        // complete_order returns three resultsets:
        //   1. The affected order. Always a single row.
        //   2. A collection of line items for the affected order.
        //   3. An empty resultset describing the effects of the CALL statement
        mysql::static_results<order, order_item, empty> result;
        conn.execute(stmt.bind(args.order_id), result);

        // Print the results to stdout
        std::cout << "Completed order\n";
        print_order_with_items(result.rows<0>()[0], result.rows<1>());
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

#else

int main()
{
    std::cout << "Sorry, your compiler doesn't have the required capabilities to run this example"
              << std::endl;
}

#endif

//]
