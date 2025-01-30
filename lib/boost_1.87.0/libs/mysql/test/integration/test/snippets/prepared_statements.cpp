//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/results.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/test/unit_test.hpp>

#include "test_integration/snippets/get_connection.hpp"

#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
#include <optional>
#endif

using namespace boost::mysql;

namespace {

int get_int_value_from_user() { return 42; }

//[prepared_statements_execute
// description, price and show_in_store are not trusted, since they may
// have been read from a file or an HTTP endpoint
void insert_product(
    tcp_connection& conn,
    const statement& stmt,
    string_view description,
    int price,
    bool show_in_store
)
{
    results result;
    conn.execute(stmt.bind(description, price, show_in_store), result);
}
//]

#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
//[prepared_statements_execute_null
// description, price and show_in_store are not trusted, since they may
// have been read from a file or an HTTP endpoint
void insert_product(
    tcp_connection& conn,
    const statement& stmt,
    std::optional<string_view> description,
    int price,
    bool show_in_store
)
{
    // If description has a value, a string will be sent to the server; otherwise, a NULL will
    results result;
    conn.execute(stmt.bind(description, price, show_in_store), result);
}
//]
#endif

//[prepared_statements_execute_iterator_range
void exec_statement(tcp_connection& conn, const statement& stmt, const std::vector<field>& params)
{
    results result;
    conn.execute(stmt.bind(params.begin(), params.end()), result);
}
//]

BOOST_AUTO_TEST_CASE(section_prepared_statements)
{
    auto& conn = test::get_connection();
    {
        //[prepared_statements_prepare
        // Table setup
        const char* table_definition = R"%(
            CREATE TEMPORARY TABLE products (
                id INT PRIMARY KEY AUTO_INCREMENT,
                description VARCHAR(256),
                price INT NOT NULL,
                show_in_store TINYINT
            )
        )%";
        results result;
        conn.execute(table_definition, result);

        // Prepare a statement to insert into this table
        statement stmt = conn.prepare_statement(
            "INSERT INTO products (description, price, show_in_store) VALUES (?, ?, ?)"
        );
        //]

        // Run the functions to verify that everything works
        insert_product(conn, stmt, string_view("This is a product"), 2000, true);
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
        insert_product(conn, stmt, std::optional<string_view>(), 2000, true);
#endif
        exec_statement(conn, stmt, {field_view("abc"), field_view(2000), field_view(1)});
        conn.execute("DROP TABLE products", result);
    }
    {
        //[prepared_statements_casting_table
        const char* table_definition = "CREATE TEMPORARY TABLE my_table(my_field TINYINT)";
        //]

        results r;
        conn.execute(table_definition, r);

        //[prepared_statements_casting_execute
        int value = get_int_value_from_user();
        auto stmt = conn.prepare_statement("INSERT INTO my_table VALUES (?)");

        results result;
        conn.execute(stmt.bind(value), result);
        //]
    }
}

}  // namespace
