//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/static_execution_state.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <iostream>
#include <string>

#include "test_integration/snippets/describe.hpp"
#include "test_integration/snippets/get_connection.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

namespace {

static std::string get_company_id() { return "HGS"; }

BOOST_AUTO_TEST_CASE(section_multi_function)
{
    auto& conn = get_connection();

    {
        //[multi_function_setup
        const char* table_definition = R"%(
            CREATE TEMPORARY TABLE posts (
                id INT PRIMARY KEY AUTO_INCREMENT,
                title VARCHAR (256) NOT NULL,
                body TEXT NOT NULL
            )
        )%";
        //]

        results result;
        conn.execute(table_definition, result);
        conn.execute(
            R"%(
                INSERT INTO posts (title, body) VALUES
                    ('Post 1', 'A very long post body'),
                    ('Post 2', 'An even longer post body')
            )%",
            result
        );

        //[multi_function_dynamic_start
        // st will hold information about the operation being executed.
        // It must be passed to any successive operations for this execution
        execution_state st;

        // Sends the query and reads response and meta, but not the rows
        conn.start_execution("SELECT title, body FROM posts", st);
        //]

        //[multi_function_dynamic_read
        // st.complete() returns true once the OK packet is received
        while (!st.complete())
        {
            // row_batch will be valid until conn performs the next network operation
            rows_view row_batch = conn.read_some_rows(st);

            for (row_view post : row_batch)
            {
                // Process post as required
                std::cout << "Title:" << post.at(0) << std::endl;
            }
        }
        //]
    }
#ifdef BOOST_MYSQL_CXX14
    {
        //[multi_function_static_start
        // st will hold information about the operation being executed.
        // It must be passed to any successive operations for this execution
        static_execution_state<post> st;

        // Sends the query and reads response and meta, but not the rows.
        // If there is any schema mismatch between the declared row type and
        // what the server returned, start_execution will detect it and fail
        conn.start_execution("SELECT id, title, body FROM posts", st);
        //]

        //[multi_function_static_read
        // storage will be filled with the read rows. You can use any other contiguous range.
        std::array<post, 20> posts;

        // st.complete() returns true once the OK packet is received
        while (!st.complete())
        {
            std::size_t read_rows = conn.read_some_rows(st, boost::span<post>(posts));
            for (const post& p : boost::span<post>(posts.data(), read_rows))
            {
                // Process post as required
                std::cout << "Title " << p.title << std::endl;
            }
        }
        //]

        results result;
        conn.execute("DROP TABLE posts", result);
    }
#endif
    {
        //[multi_function_stored_procedure_dynamic
        // Get the company ID to retrieve, possibly from the user
        std::string company_id = get_company_id();

        // Call the procedure
        execution_state st;
        statement stmt = conn.prepare_statement("CALL get_employees(?)");
        conn.start_execution(stmt.bind(company_id), st);

        // The above code will generate 3 resultsets
        // Read the 1st one, which contains the matched companies
        while (st.should_read_rows())
        {
            rows_view company_batch = conn.read_some_rows(st);

            // Use the retrieved companies as required
            for (row_view company : company_batch)
            {
                std::cout << "Company: " << company.at(1).as_string() << "\n";
            }
        }

        // Move on to the 2nd one, containing the employees for these companies
        conn.read_resultset_head(st);
        while (st.should_read_rows())
        {
            rows_view employee_batch = conn.read_some_rows(st);

            // Use the retrieved employees as required
            for (row_view employee : employee_batch)
            {
                std::cout << "Employee " << employee.at(0).as_string() << " " << employee.at(1).as_string()
                          << "\n";
            }
        }

        // The last one is an empty resultset containing information about the
        // CALL statement itself. We're not interested in this
        conn.read_resultset_head(st);
        assert(st.complete());
        //]
    }
#ifdef BOOST_MYSQL_CXX14
    {
        //[multi_function_stored_procedure_static
        // Get the company ID to retrieve, possibly from the user
        std::string company_id = get_company_id();

        // Our procedure generates three resultsets. We must pass each row type
        // to static_execution_state as template parameters
        using empty = std::tuple<>;
        static_execution_state<company, employee, empty> st;

        // Call the procedure
        statement stmt = conn.prepare_statement("CALL get_employees(?)");
        conn.start_execution(stmt.bind(company_id), st);

        // Read the 1st one, which contains the matched companies
        std::array<company, 5> companies;
        while (st.should_read_rows())
        {
            std::size_t read_rows = conn.read_some_rows(st, boost::span<company>(companies));

            // Use the retrieved companies as required
            for (const company& c : boost::span<company>(companies.data(), read_rows))
            {
                std::cout << "Company: " << c.name << "\n";
            }
        }

        // Move on to the 2nd one, containing the employees for these companies
        conn.read_resultset_head(st);
        std::array<employee, 20> employees;
        while (st.should_read_rows())
        {
            std::size_t read_rows = conn.read_some_rows(st, boost::span<employee>(employees));

            // Use the retrieved companies as required
            for (const employee& emp : boost::span<employee>(employees.data(), read_rows))
            {
                std::cout << "Employee " << emp.first_name << " " << emp.last_name << "\n";
            }
        }

        // The last one is an empty resultset containing information about the
        // CALL statement itself. We're not interested in this
        conn.read_resultset_head(st);
        assert(st.complete());
        //]
    }
#endif
}

}  // namespace