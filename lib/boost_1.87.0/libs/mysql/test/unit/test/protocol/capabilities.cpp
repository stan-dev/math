//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/impl/internal/protocol/capabilities.hpp>

#include <boost/test/unit_test.hpp>

using namespace boost::mysql::detail;

BOOST_AUTO_TEST_SUITE(test_capabilities)

constexpr capabilities rhs{CLIENT_CONNECT_WITH_DB | CLIENT_SSL | CLIENT_COMPRESS};

BOOST_AUTO_TEST_CASE(has_bit_set)
{
    capabilities caps(CLIENT_COMPRESS);
    BOOST_TEST(caps.has(CLIENT_COMPRESS));
}

BOOST_AUTO_TEST_CASE(has_bit_not_set)
{
    capabilities caps(CLIENT_COMPRESS);
    BOOST_TEST(!caps.has(CLIENT_SSL));
}

BOOST_AUTO_TEST_CASE(has_multiple_bits_set)
{
    capabilities caps(CLIENT_CONNECT_WITH_DB | CLIENT_SSL | CLIENT_COMPRESS);
    for (int i = 0; i < 32; ++i)
    {
        std::uint32_t cap_bit = 1 << i;
        bool is_set = cap_bit == CLIENT_CONNECT_WITH_DB || cap_bit == CLIENT_SSL ||
                      cap_bit == CLIENT_COMPRESS;
        BOOST_TEST(caps.has(cap_bit) == is_set);
    }
}

BOOST_AUTO_TEST_CASE(has_all_has_none)
{
    capabilities lhs(0);
    BOOST_TEST(!lhs.has_all(rhs));
}

BOOST_AUTO_TEST_CASE(has_all_has_some_but_not_all)
{
    capabilities lhs(CLIENT_CONNECT_WITH_DB | CLIENT_COMPRESS);
    BOOST_TEST(!lhs.has_all(rhs));
}

BOOST_AUTO_TEST_CASE(has_all_has_some_but_not_all_plus_unrelated)
{
    capabilities lhs(CLIENT_CONNECT_WITH_DB | CLIENT_COMPRESS | CLIENT_TRANSACTIONS);
    BOOST_TEST(!lhs.has_all(rhs));
}

BOOST_AUTO_TEST_CASE(has_all_has_only_the_requested_ones)
{
    capabilities lhs(rhs);
    BOOST_TEST(lhs.has_all(rhs));
}

BOOST_AUTO_TEST_CASE(has_all_has_the_requested_ones_and_others)
{
    capabilities lhs = rhs | capabilities(CLIENT_TRANSACTIONS);
    BOOST_TEST(lhs.has_all(rhs));
}

BOOST_AUTO_TEST_SUITE_END()  // test_capabilities
