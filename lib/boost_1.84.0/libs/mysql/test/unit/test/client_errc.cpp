//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/stringize.hpp"

using namespace boost::mysql;

namespace {

BOOST_AUTO_TEST_SUITE(test_client_errc)

BOOST_AUTO_TEST_SUITE(error_to_string_)

BOOST_AUTO_TEST_CASE(regular)
{
    BOOST_TEST(error_code(client_errc::sequence_number_mismatch).message() == "Mismatched sequence numbers");
}

BOOST_AUTO_TEST_CASE(unknown_error)
{
    BOOST_TEST(error_code(static_cast<client_errc>(0xfffefdfc)).message() == "<unknown MySQL client error>");
}

BOOST_AUTO_TEST_CASE(coverage)
{
    // Check that no value causes problems.
    // Ensure that all branches of the switch/case are covered
    for (int i = 0; i < 20; ++i)
    {
        BOOST_CHECK_NO_THROW(error_code(static_cast<client_errc>(i)).message());
    }
}
BOOST_AUTO_TEST_SUITE_END()  // error_to_string_

BOOST_AUTO_TEST_CASE(error_code_from_errc)
{
    error_code code(client_errc::protocol_value_error);
    BOOST_TEST(code.value() == static_cast<int>(client_errc::protocol_value_error));
    BOOST_TEST(&code.category() == &boost::mysql::get_client_category());  // categories are not printable
}

BOOST_AUTO_TEST_SUITE_END()  // test_client_errc

}  // namespace
