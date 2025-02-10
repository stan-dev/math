//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/test/unit_test.hpp>

#include "test_integration/snippets/get_connection.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

namespace {

BOOST_AUTO_TEST_CASE(section_metadata)
{
    auto& conn = get_connection();

    //[metadata
    // By default, a connection has metadata_mode::minimal
    results result;
    conn.execute("SELECT 1 AS my_field", result);
    string_view colname = result.meta()[0].column_name();

    // colname will be empty because conn.meta_mode() == metadata_mode::minimal
    BOOST_TEST(colname == "");

    // If you are using metadata names, set the connection's metadata_mode
    conn.set_meta_mode(metadata_mode::full);
    conn.execute("SELECT 1 AS my_field", result);
    colname = result.meta()[0].column_name();
    BOOST_TEST(colname == "my_field");
    //]

    conn.set_meta_mode(metadata_mode::minimal);
}

}  // namespace
