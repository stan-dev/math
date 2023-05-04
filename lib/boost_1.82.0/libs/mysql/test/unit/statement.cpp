//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/statement.hpp>

#include <boost/test/unit_test.hpp>

#include "create_statement.hpp"
#include "test_common.hpp"

using namespace boost::mysql::test;
using boost::mysql::statement;

namespace {

BOOST_AUTO_TEST_SUITE(test_statement_)

BOOST_AUTO_TEST_CASE(default_ctor)
{
    statement stmt;
    BOOST_TEST(!stmt.valid());
}

BOOST_AUTO_TEST_CASE(member_fns)
{
    auto stmt = create_statement(3, 1);

    BOOST_TEST(stmt.valid());
    BOOST_TEST(stmt.num_params() == 3u);
    BOOST_TEST(stmt.id() == 1u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
