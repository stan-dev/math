//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/resultset_view.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/check_meta.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_resultset_view)

BOOST_AUTO_TEST_CASE(null_view)
{
    resultset_view v;
    BOOST_TEST(!v.has_value());
}

BOOST_AUTO_TEST_CASE(valid_view)
{
    results result;
    exec_access(get_iface(result))
        .meta({meta_builder().type(column_type::tinyint).build_coldef()})
        .row(42)
        .ok(ok_builder().affected_rows(4).last_insert_id(5).warnings(6).info("2nd").out_params(true).build());

    auto v = result.at(0);
    BOOST_TEST_REQUIRE(v.has_value());
    BOOST_TEST(v.rows() == makerows(1, 42));
    check_meta(v.meta(), {column_type::tinyint});
    BOOST_TEST(v.affected_rows() == 4u);
    BOOST_TEST(v.last_insert_id() == 5u);
    BOOST_TEST(v.warning_count() == 6u);
    BOOST_TEST(v.info() == "2nd");
    BOOST_TEST(v.is_out_params());
}

BOOST_AUTO_TEST_SUITE_END()
