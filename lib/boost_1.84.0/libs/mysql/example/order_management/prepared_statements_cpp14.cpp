//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_prepared_statements_cpp14

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
 * We'll be using the static interface to retrieve results from MySQL.
 * This makes use of the static_results<RowType> class template.
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
#include <stdexcept>
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

// An order with its line items. This record type is returned by JOINs
// from the orders and order_items tables. We use this type to retrieve both
// an order and its line items in a single operation.
// If the order contains no line items, the item_xxx fields are NULL.
struct order_with_items
{
    // The ID of the order
    std::int64_t order_id;

    // The status of the order
    std::string order_status;

    // The ID of the line item, or NULL if the order doesn't have any
    boost::optional<std::int64_t> item_id;

    // The number of units of this product that the user wants to buy,
    // or NULL if the order doesn't have line items
    boost::optional<std::int64_t> item_quantity;

    // The product's unit price, in cents of USD, or NULL if the order
    // doesn't have line items
    boost::optional<std::int64_t> item_unit_price;

    bool has_item() const
    {
        return item_id.has_value() && item_quantity.has_value() && item_unit_price.has_value();
    }
};
BOOST_DESCRIBE_STRUCT(
    order_with_items,
    (),
    (order_id, order_status, item_id, item_quantity, item_unit_price)
);

// An empty row type. This can be used to describe empty resultsets,
// like the ones returned by INSERT or CALL.
using empty = std::tuple<>;

// This visitor executes a sub-command and prints the results to stdout.
struct visitor
{
    mysql::tcp_ssl_connection& conn;

    static void print_order(const order& ord)
    {
        std::cout << "Order: id=" << ord.id << ", status=" << ord.status << '\n';
    }

    // Retrieves an order with its items. If the order exists, at least one record is returned.
    // If the order has line items, a record per item is returned. If the order has no items,
    // a single record is returned, and it will have its item_xxx fields set to NULL.
    mysql::static_results<order_with_items> get_order_with_items(std::int64_t order_id) const
    {
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

        mysql::static_results<order_with_items> result;
        conn.execute(stmt.bind(order_id), result);
        return result;
    }

    // Prints an order with its line items to stdout
    static void print_order_with_items(boost::span<const order_with_items> ord_items)
    {
        assert(!ord_items.empty());

        // Print the order
        std::cout << "Order: id=" << ord_items[0].order_id << ", status=" << ord_items[0].order_status
                  << '\n';

        // Print the items
        if (!ord_items[0].has_item())
        {
            std::cout << "No line items\n";
        }
        else
        {
            for (const auto& item : ord_items)
            {
                std::cout << "  Line item: id=" << *item.item_id << ", quantity=" << *item.item_quantity
                          << ", unit_price=" << *item.item_unit_price / 100.0 << "$\n";
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

        // The product struct describes the shape of the rows that
        // we expect the server to send.
        mysql::static_results<product> products;
        conn.execute(stmt.bind(args.search), products);

        // Print the results to stdout
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
        // Since this is an INSERT, we don't expect any row to be returned.
        // empty is an alias for std::tuple<>, which tells static_results to expect
        // an empty resultset.
        mysql::static_results<empty> result;
        conn.execute("INSERT INTO orders VALUES ()", result);

        // We can use static_results::last_insert_id() to retrieve the ID of the newly
        // created object. last_insert_id() returns always a uint64_t. Our schema uses
        // plain INTs for the id field, so this cast is safe.
        order ord{static_cast<std::int64_t>(result.last_insert_id()), "draft"};
        print_order(ord);
    }

    // get-order <order-id>: retrieves order details
    void operator()(const get_order_args& args) const
    {
        // Retrieve the order with its items
        mysql::static_results<order_with_items> result = get_order_with_items(args.order_id);

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
        mysql::static_results<order> result;
        conn.execute("SELECT id, `status` FROM orders", result);

        // Print the results to stdout
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
        // We will need to run several statements atomically, so we start a transaction.
        mysql::static_results<empty> empty_results;
        conn.execute("START TRANSACTION", empty_results);

        // To add a line item, we require the order to be in a draft status. Get the order to check this fact.
        mysql::statement stmt = conn.prepare_statement("SELECT id, `status` FROM orders WHERE id = ?");
        mysql::static_results<order> orders;
        conn.execute(stmt.bind(args.order_id), orders);
        if (orders.rows().empty())
        {
            // There is no such order
            throw std::runtime_error("Order with id=" + std::to_string(args.order_id) + " not found");
        }
        else if (orders.rows()[0].status != "draft")
        {
            // The order is no longer editable
            throw std::runtime_error("Order with id=" + std::to_string(args.order_id) + " is not editable");
        }

        // Insert the new line item. If the given product does not exist, the INSERT will fail
        // because of product_id's FOREIGN KEY constraint.
        stmt = conn.prepare_statement(
            "INSERT INTO order_items (order_id, product_id, quantity) VALUES (?, ?, ?)"
        );
        conn.execute(stmt.bind(args.order_id, args.product_id, args.quantity), empty_results);

        // We can use static_results::last_insert_id to get the ID of the new line item.
        auto new_line_item_id = empty_results.last_insert_id();

        // Retrieve the full order details
        mysql::static_results<order_with_items> order_results = get_order_with_items(args.order_id);

        // We're done - commit the transaction
        conn.execute("COMMIT", empty_results);

        // Print the results to stdout
        std::cout << "Created line item: id=" << new_line_item_id << "\n";
        print_order_with_items(order_results.rows());
    }

    // remove-line-item <line-item-id>: removes an item from an order
    void operator()(const remove_line_item_args& args) const
    {
        // We will need to run several statements atomically, so we start a transaction.
        mysql::static_results<empty> empty_results;
        conn.execute("START TRANSACTION", empty_results);

        // To remove a line item, we require the order to be in a draft status. Get the order to check this
        // fact.
        mysql::static_results<order> orders;
        auto stmt = conn.prepare_statement(R"%(
            SELECT orders.id, orders.`status`
            FROM orders
            JOIN order_items items ON (orders.id = items.order_id)
            WHERE items.id = ?
        )%");
        conn.execute(stmt.bind(args.line_item_id), orders);
        if (orders.rows().empty())
        {
            // The query hasn't matched any row - the supplied line item ID is not valid
            throw std::runtime_error(
                "The order item with id=" + std::to_string(args.line_item_id) + " does not exist"
            );
        }
        const order& ord = orders.rows()[0];
        if (ord.status != "draft")
        {
            // The order is no longer editable
            throw std::runtime_error("The order is not in an editable state");
        }

        // Remove the line item
        stmt = conn.prepare_statement("DELETE FROM order_items WHERE id = ?");
        conn.execute(stmt.bind(args.line_item_id), empty_results);

        // Retrieve the full order details
        mysql::static_results<order_with_items> order_results = get_order_with_items(ord.id);

        // We're done - commit the transaction
        conn.execute("COMMIT", empty_results);

        // Print results to stdout
        std::cout << "Removed line item from order\n";
        print_order_with_items(order_results.rows());
    }

    // checkout-order <order-id>: marks an order as ready for checkout
    void operator()(const checkout_order_args& args) const
    {
        // We will need to run several statements atomically, so we start a transaction.
        mysql::static_results<empty> empty_results;
        conn.execute("START TRANSACTION", empty_results);

        // To checkout an order, we require it to be in a draft status. Check this fact.
        mysql::statement stmt = conn.prepare_statement("SELECT id, `status` FROM orders WHERE id = ?");
        mysql::static_results<order> orders;
        conn.execute(stmt.bind(args.order_id), orders);
        if (orders.rows().empty())
        {
            // No order matched
            throw std::runtime_error("Order with id=" + std::to_string(args.order_id) + " not found");
        }
        else if (orders.rows()[0].status != "draft")
        {
            // The order is no longer editable
            throw std::runtime_error(
                "Order with id=" + std::to_string(args.order_id) + " cannot be checked out"
            );
        }

        // Update the order status
        stmt = conn.prepare_statement("UPDATE orders SET `status` = 'pending_payment' WHERE id = ?");
        conn.execute(stmt.bind(args.order_id), empty_results);

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
        mysql::static_results<std::tuple<std::uint64_t>> amount_results;
        conn.execute(stmt.bind(args.order_id), amount_results);
        std::uint64_t total_amount = std::get<0>(amount_results.rows()[0]);

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
        mysql::static_results<order_with_items> order_results = get_order_with_items(args.order_id);

        // We're done - commit the transaction
        conn.execute("COMMIT", empty_results);

        // Print the results to stdout
        std::cout << "Checked out order. The total amount to pay is: " << total_amount / 100.0 << "$\n";
        print_order_with_items(order_results.rows());
    }

    // complete-order <order-id>: marks an order as completed
    void operator()(const complete_order_args& args) const
    {
        // We will need to run several statements atomically, so we start a transaction.
        mysql::static_results<empty> empty_results;
        conn.execute("START TRANSACTION", empty_results);

        // To complete an order, we require it to be in a pending_payment status. Check this fact.
        auto stmt = conn.prepare_statement("SELECT id, `status` FROM orders WHERE id = ?");
        mysql::static_results<order> orders;
        conn.execute(stmt.bind(args.order_id), orders);
        if (orders.rows().empty())
        {
            // Order not found
            throw std::runtime_error("Order with id=" + std::to_string(args.order_id) + " not found");
        }
        else if (orders.rows()[0].status != "pending_payment")
        {
            throw std::runtime_error(
                "Order with id=" + std::to_string(args.order_id) + " is not in pending_payment status"
            );
        }

        // Update status
        stmt = conn.prepare_statement("UPDATE orders SET `status` = 'complete' WHERE id = ?");
        conn.execute(stmt.bind(args.order_id), empty_results);

        // Retrieve the full order details
        mysql::static_results<order_with_items> order_results = get_order_with_items(args.order_id);

        // We're done - commit the transaction
        conn.execute("COMMIT", empty_results);

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

#else

int main()
{
    std::cout << "Sorry, your compiler doesn't have the required capabilities to run this example"
              << std::endl;
}

#endif

//]
