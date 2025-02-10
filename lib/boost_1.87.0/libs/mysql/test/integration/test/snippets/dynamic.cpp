//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/field.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

#include "test_integration/snippets/get_connection.hpp"

using namespace boost::mysql;

namespace {

BOOST_AUTO_TEST_CASE(section_dynamic)
{
    auto& conn = test::get_connection();

    {
        //[dynamic_views
        // Populate a results object
        results result;
        conn.execute("SELECT 'Hello world'", result);

        // results::rows() returns a rows_view. The underlying memory is owned by the results object
        rows_view all_rows = result.rows();

        // Indexing a rows_view yields a row_view. The underlying memory is owned by the results object
        row_view first_row = all_rows.at(0);

        // Indexing a row_view yields a field_view. The underlying memory is owned by the results object
        field_view first_field = first_row.at(0);  // Contains the string "Hello world"

        //]
        BOOST_TEST(first_field.as_string() == "Hello world");

        //[dynamic_taking_ownership
        // You may use all_rows_owning after result has gone out of scope
        rows all_rows_owning{all_rows};

        // You may use first_row_owning after result has gone out of scope
        row first_row_owning{first_row};

        // You may use first_field_owning after result has gone out of scope
        field first_field_owning{first_field};
        //]
    }
    {
        //[dynamic_using_fields
        results result;
        conn.execute("SELECT 'abc', 42", result);

        // Obtain a field's underlying value using the is_xxx and get_xxx accessors
        field_view f = result.rows().at(0).at(0);  // f points to the string "abc"
        if (f.is_string())
        {
            // we know it's a string, unchecked access
            string_view s = f.get_string();
            std::cout << s << std::endl;  // Use the string as required
        }
        else
        {
            // Oops, something went wrong - schema mismatch?
        }

        // Alternative: use the as_xxx accessor
        f = result.rows().at(0).at(1);
        std::int64_t value = f.as_int64();  // Checked access. Throws if f doesn't contain an int
        std::cout << value << std::endl;    // Use the int as required

        //]
    }
    {
        //[dynamic_handling_nulls
        results result;

        // Create some test data
        conn.execute(
            R"%(
                CREATE TEMPORARY TABLE products (
                    id VARCHAR(50) PRIMARY KEY,
                    description VARCHAR(256)
                )
            )%",
            result
        );
        conn.execute("INSERT INTO products VALUES ('PTT', 'Potatoes'), ('CAR', NULL)", result);

        // Retrieve the data. Note that some fields are NULL
        conn.execute("SELECT id, description FROM products", result);

        for (row_view r : result.rows())
        {
            field_view description_fv = r.at(1);
            if (description_fv.is_null())
            {
                // Handle the NULL value
                // Note: description_fv.is_string() will return false here; NULL is represented as a separate
                // type
                std::cout << "No description for product_id " << r.at(0) << std::endl;
            }
            else
            {
                // Handle the non-NULL case. Get the underlying value and use it as you want
                // If there is any schema mismatch (and description was not defined as VARCHAR), this will
                // throw
                string_view description = description_fv.as_string();

                // Use description as required
                std::cout << "product_id " << r.at(0) << ": " << description << std::endl;
            }
        }
        //]

        conn.execute("DROP TABLE products", result);
    }
    {
        //[dynamic_field_accessor_references
        field f("my_string");            // constructs a field that owns the string "my_string"
        std::string& s = f.as_string();  // s points into f's storage
        s.push_back('2');                // f now holds "my_string2"

        //]

        BOOST_TEST(s == "my_string2");
    }
    {
        //[dynamic_field_assignment
        field f("my_string");  // constructs a field that owns the string "my_string"
        f = 42;                // destroys "my_string" and stores the value 42 as an int64

        //]

        BOOST_TEST(f.as_int64() == 42);
    }
}

}  // namespace
