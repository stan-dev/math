//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SNIPPETS_DESCRIBE_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SNIPPETS_DESCRIBE_HPP

#include <boost/describe/class.hpp>
#include <boost/optional/optional.hpp>

#include <cstdint>
#include <string>
#include <tuple>

namespace boost {
namespace mysql {
namespace test {

//[describe_post
// We can use a plain struct with ints and strings to describe our rows.
// This must be placed at the namespace level
struct post
{
    int id;
    std::string title;
    std::string body;
};

// We use BOOST_DESCRIBE_STRUCT to add reflection capabilities to post.
// We must list all the fields that should be populated by Boost.MySQL
BOOST_DESCRIBE_STRUCT(post, (), (id, title, body))
//]

//[describe_stored_procedures
// Describes the first resultset
struct company
{
    std::string id;
    std::string name;
    std::string tax_id;
};
BOOST_DESCRIBE_STRUCT(company, (), (id, name, tax_id))

// Describes the second resultset
struct employee
{
    std::string first_name;
    std::string last_name;
    boost::optional<std::uint64_t> salary;
};
BOOST_DESCRIBE_STRUCT(employee, (), (first_name, last_name, salary))

// The last resultset will always be empty.
// We can use an empty tuple to represent it.
using empty = std::tuple<>;
//]

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
