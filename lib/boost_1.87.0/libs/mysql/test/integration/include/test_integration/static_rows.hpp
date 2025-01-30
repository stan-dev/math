//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_STATIC_ROWS_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_STATIC_ROWS_HPP

#include <boost/describe/class.hpp>
#include <boost/describe/operators.hpp>
#include <boost/optional/optional.hpp>
#include <boost/optional/optional_io.hpp>

#include <tuple>

namespace boost {
namespace mysql {
namespace test {

struct row_multifield
{
    boost::optional<float> field_nullable;
    std::int32_t field_int;
    std::string field_varchar;
};
BOOST_DESCRIBE_STRUCT(row_multifield, (), (field_nullable, field_int, field_varchar))

struct row_multifield_bad
{
    std::string field_varchar;
    float field_nullable;
    std::string field_int;
    int field_missing;
};
BOOST_DESCRIBE_STRUCT(row_multifield_bad, (), (field_varchar, field_nullable, field_int, field_missing))

struct row_2fields
{
    boost::optional<int> id;
    boost::optional<std::string> field_varchar;
};
BOOST_DESCRIBE_STRUCT(row_2fields, (), (id, field_varchar))

using empty = std::tuple<>;

#ifdef BOOST_DESCRIBE_CXX14
using boost::describe::operators::operator==;
using boost::describe::operators::operator<<;
#endif

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
