//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/pfr.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/static_results.hpp>

#include <boost/config.hpp>
#include <boost/test/unit_test.hpp>

#include <tuple>

#include "test_integration/snippets/describe.hpp"
#include "test_integration/snippets/get_connection.hpp"

#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
#include <optional>
#endif

using namespace boost::mysql;

// PFR types can't be placed in anonymous namespaces
namespace snippets_static {

//[describe_statistics
struct statistics
{
    std::string company;
    double average;
    double max_value;
};
BOOST_DESCRIBE_STRUCT(statistics, (), (company, average, max_value))
//]

#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
//[describe_post_v2
struct post_v2
{
    int id;
    std::string title;
    std::optional<std::string> body;  // body may be NULL
};
BOOST_DESCRIBE_STRUCT(post_v2, (), (id, title, body))
//]
#endif

//[describe_post_pfr
// post_v3 doesn't contain any metadata - we're not using BOOST_DESCRIBE_STRUCT here
struct post_v3
{
    int id;
    std::string title;
    std::string body;
};
//]

//
// Comparison table - we want all type definitions to be similar
//
namespace descr_type {
//[static_comparison_describe_struct
// Definition should be at namespace scope
struct post
{
    int id;
    std::string title;
    std::string body;
};
BOOST_DESCRIBE_STRUCT(post, (), (id, title, body))
//]
}  // namespace descr_type

namespace pfr_type {
//[static_comparison_pfr_struct
// Definition should be at namespace scope
struct post
{
    int id;
    std::string title;
    std::string body;
};
//]
}  // namespace pfr_type

BOOST_AUTO_TEST_CASE(section_static)
{
#ifdef BOOST_MYSQL_CXX14

    auto& conn = test::get_connection();

    {
        //[static_setup
        constexpr const char* table_definition = R"%(
            CREATE TEMPORARY TABLE posts (
                id INT PRIMARY KEY AUTO_INCREMENT,
                title VARCHAR (256) NOT NULL,
                body TEXT NOT NULL
            )
        )%";
        //]

        results r;
        conn.execute(table_definition, r);
    }

    {
        using test::post;

        //[static_query
        static_results<post> result;
        conn.execute("SELECT id, title, body FROM posts", result);

        for (const post& p : result.rows())
        {
            // Process the post as required
            std::cout << "Title: " << p.title << "\n" << p.body << "\n";
        }
        //]
    }
    {
        //[static_field_order
        // Summing 0e0 is MySQL way to cast a DECIMAL field to DOUBLE
        constexpr const char* sql = R"%(
            SELECT
                IFNULL(AVG(salary), 0.0) + 0e0 AS average,
                IFNULL(MAX(salary), 0.0) + 0e0 AS max_value,
                company_id AS company
            FROM employee
            GROUP BY company_id
        )%";

        static_results<statistics> result;
        conn.execute(sql, result);
        //]
    }
#if BOOST_PFR_CORE_NAME_ENABLED
    {
        //[static_pfr_by_name
        // pfr_by_name is a marker type. It tells static_results to use
        // Boost.PFR for reflection, instead of Boost.Describe.
        static_results<pfr_by_name<post_v3>> result;

        // As with Boost.Describe, query fields are matched to struct
        // members by name. This means that the fields in the query
        // may appear in any order.
        conn.execute("SELECT body, id, title FROM posts", result);

        // Note that result.rows() is a span of post_v3 objects,
        // rather than pfr_by_name<post_v3> objects. post_v3
        // is the underlying row type for pfr_by_name<post_v3>
        for (const post_v3& p : result.rows())
        {
            // Process the post as required
            std::cout << "Title: " << p.title << "\n" << p.body << "\n";
        }
        //]
    }
#endif
#if BOOST_PFR_USE_CPP17
    {
        //[static_pfr_by_position
        // pfr_by_position is another marker type.
        // Fields in post_v3 must appear in the same order as in the query,
        // as matching will be done by position.
        static_results<pfr_by_position<post_v3>> result;
        conn.execute("SELECT id, title, body FROM posts", result);

        // The underlying row type is post_v3
        for (const post_v3& p : result.rows())
        {
            // Process the post as required
            std::cout << "Title: " << p.title << "\n" << p.body << "\n";
        }
        //]
    }
#endif
    {
        //[static_tuples
        static_results<std::tuple<std::int64_t>> result;
        conn.execute("SELECT COUNT(*) FROM employee", result);
        std::cout << "Number of employees: " << std::get<0>(result.rows()[0]) << "\n";
        //]
    }
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
    {
        //[static_nulls_table
        constexpr const char* table_definition = R"%(
            CREATE TEMPORARY TABLE posts_v2 (
                id INT PRIMARY KEY AUTO_INCREMENT,
                title VARCHAR (256) NOT NULL,
                body TEXT
            )
        )%";
        //]

        // Verify that post_v2's definition is correct
        results r;
        conn.execute(table_definition, r);
        static_results<post_v2> result;
        conn.execute("SELECT * FROM posts_v2", result);
        conn.execute("DROP TABLE posts_v2", r);
    }
#endif  // BOOST_NO_CXX17_HDR_OPTIONAL
    {
        using namespace descr_type;
        //[static_comparison_describe
        // Usage
        static_results<post> result;
        conn.execute("SELECT title, body, id FROM posts", result);
        //]
    }
#if BOOST_PFR_CORE_NAME_ENABLED
    {
        using namespace pfr_type;
        //[static_comparison_pfr_by_name
        // Usage
        static_results<pfr_by_name<post>> result;
        conn.execute("SELECT title, body, id FROM posts", result);
        //]
    }
#endif
#if BOOST_PFR_USE_CPP17
    {
        using namespace pfr_type;
        //[static_comparison_pfr_by_position
        // Usage
        static_results<pfr_by_position<post>> result;
        conn.execute("SELECT id, title, body FROM posts", result);
        //]
    }
#endif
    {
        //[static_comparison_tuples
        using tuple_t = std::tuple<int, std::string, std::string>;
        static_results<tuple_t> result;
        conn.execute("SELECT id, title, body FROM posts", result);
        //]
    }

    {
        results r;
        conn.execute("DROP TABLE posts", r);
    }

#endif  // BOOST_MYSQL_CXX14
}

}  // namespace snippets_static