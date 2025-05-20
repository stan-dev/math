//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_TEST_EXECUTION_PROCESSOR_STATIC_EXECUTION_PROCESSOR_HELPERS_HPP
#define BOOST_MYSQL_TEST_UNIT_TEST_EXECUTION_PROCESSOR_STATIC_EXECUTION_PROCESSOR_HELPERS_HPP

#include <boost/mysql/metadata.hpp>
#include <boost/mysql/metadata_collection_view.hpp>

#include <boost/describe/class.hpp>
#include <boost/describe/operators.hpp>

// Throughout execution_processor tests we use a set of common values for rows.
// Definitions here to reduce duplication

namespace boost {
namespace mysql {
namespace test {

// Row types
struct row1
{
    std::string fvarchar;
    std::int16_t ftiny;
};
BOOST_DESCRIBE_STRUCT(row1, (), (fvarchar, ftiny))

using row1_tuple = std::tuple<std::int16_t, std::string>;

struct row2
{
    std::int64_t fbigint;
};
BOOST_DESCRIBE_STRUCT(row2, (), (fbigint))

using row2_tuple = std::tuple<std::int64_t>;

// This is doing field reordering by name
struct row3
{
    double fdouble;
    std::int8_t ftiny;
    float ffloat;
};
BOOST_DESCRIBE_STRUCT(row3, (), (fdouble, ftiny, ffloat))

using row3_tuple = std::tuple<float, double, std::int8_t>;

struct empty
{
};
BOOST_DESCRIBE_STRUCT(empty, (), ())

// For tests verifying that field selection works
struct row3_selection
{
    std::int8_t ftiny;
    float ffloat;
};
BOOST_DESCRIBE_STRUCT(row3_selection, (), (ftiny, ffloat))
using row3_selection_tuple = std::tuple<float, double>;

#ifdef BOOST_DESCRIBE_CXX14
using boost::describe::operators::operator==;
using boost::describe::operators::operator<<;
#endif

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
