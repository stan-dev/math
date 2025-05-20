//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/resultset_view.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/static_results.hpp>

#include <boost/core/ignore_unused.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "test_common/ci_server.hpp"
#include "test_integration/snippets/describe.hpp"
#include "test_integration/snippets/get_connection.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

namespace {

std::int64_t get_employee_id() { return 42; }

BOOST_AUTO_TEST_CASE(section_multi_resultset)
{
    auto& conn = get_connection();

    {
        //[multi_resultset_call_dynamic

        // We're using the dynamic interface. results can stored multiple resultsets
        results result;

        // The procedure parameter, employee_id, will likely be obtained from an untrusted source,
        // so we will use a prepared statement
        statement get_employee_stmt = conn.prepare_statement("CALL get_employees(?)");

        // Obtain the parameters required to call the statement, e.g. from a file or HTTP message
        std::int64_t employee_id = get_employee_id();

        // Call the statement
        conn.execute(get_employee_stmt.bind(employee_id), result);

        // results can be used as a random-access collection of resultsets.
        // result.at(0).rows() returns the matched companies, if any
        rows_view matched_company = result.at(0).rows();

        // We can do the same to access the matched employees
        rows_view matched_employees = result.at(1).rows();

        // Use matched_company and matched_employees as required
        //]

        boost::ignore_unused(matched_company);
        boost::ignore_unused(matched_employees);
    }
#ifdef BOOST_MYSQL_CXX14
    {
        //[multi_resultset_call_static
        // We must list all the resultset types the operation returns as template arguments
        static_results<company, employee, empty> result;
        conn.execute("CALL get_employees('HGS')", result);

        // We can use rows<0>() to access the rows for the first resultset
        if (result.rows<0>().empty())
        {
            std::cout << "Company not found" << std::endl;
        }
        else
        {
            const company& comp = result.rows<0>()[0];
            std::cout << "Company name: " << comp.name << ", tax_id: " << comp.tax_id << std::endl;
        }

        // rows<1>() will return the rows for the second resultset
        for (const employee& emp : result.rows<1>())
        {
            std::cout << "Employee " << emp.first_name << " " << emp.last_name << std::endl;
        }
        //]
    }
#endif
    {
        //[multi_resultset_out_params
        // To retrieve output parameters, you must use prepared statements. Text queries don't support this
        // We specify placeholders for both IN and OUT parameters
        statement stmt = conn.prepare_statement("CALL create_employee(?, ?, ?, ?)");

        // When executing the statement, we provide an actual value for the IN parameters,
        // and a dummy value for the OUT parameter. This value will be ignored, but it's required by the
        // protocol
        results result;
        conn.execute(stmt.bind("HGS", "John", "Doe", nullptr), result);

        // Retrieve output parameters. This row_view has an element per
        // OUT or INOUT parameter that used a ? placeholder
        row_view output_params = result.out_params();
        std::int64_t new_employee_id = output_params.at(0).as_int64();
        //]

        boost::ignore_unused(new_employee_id);
    }
}

BOOST_AUTO_TEST_CASE(section_multi_resultset_multi_queries)
{
    boost::asio::io_context ctx;
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());
    tcp_connection conn(ctx.get_executor());

    auto endpoint = *resolver.resolve(get_hostname(), default_port_string).begin();

    //[multi_resultset_multi_queries
    // The username and password to use
    boost::mysql::handshake_params params(
        mysql_username,         // username, as a string
        mysql_password,         // password, as a string
        "boost_mysql_examples"  // database
    );

    // Allows running multiple semicolon-separated in a single call.
    // We must set this before calling connect
    params.set_multi_queries(true);

    // Connect to the server specifying that we want support for multi-queries
    conn.connect(endpoint, params);

    // We can now use the multi-query feature.
    // This will result in three resultsets, one per query.
    results result;
    conn.execute(
        R"(
            CREATE TEMPORARY TABLE posts (
                id INT PRIMARY KEY AUTO_INCREMENT,
                title VARCHAR (256),
                body TEXT
            );
            INSERT INTO posts (title, body) VALUES ('Breaking news', 'Something happened!');
            SELECT COUNT(*) FROM posts;
        )",
        result
    );
    //]

    //[multi_resultset_results_as_collection
    // result is actually a random-access collection of resultsets.
    // The INSERT is the 2nd query, so we can access its resultset like this:
    boost::mysql::resultset_view insert_result = result.at(1);

    // A resultset has metadata, rows, and additional data, like the last insert ID:
    std::int64_t post_id = insert_result.last_insert_id();

    // The SELECT result is the third one, so we can access it like this:
    boost::mysql::resultset_view select_result = result.at(2);

    // select_result is a view that points into result.
    // We can take ownership of it using the resultset class:
    boost::mysql::resultset owning_select_result(select_result);  // valid even after result is destroyed

    // We can access rows of resultset objects as usual:
    std::int64_t num_posts = owning_select_result.rows().at(0).at(0).as_int64();
    //]

    boost::ignore_unused(post_id);
    boost::ignore_unused(num_posts);
}

}  // namespace