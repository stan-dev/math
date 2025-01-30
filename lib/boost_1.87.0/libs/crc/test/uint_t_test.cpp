// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/crc.hpp>
#include <boost/integer.hpp>
#include <boost/core/lightweight_test_trait.hpp>

int main()
{
    BOOST_TEST_TRAIT_SAME( boost::uint_t< 8>::fast, boost::crc_detail::uint_t< 8>::fast );
    BOOST_TEST_TRAIT_SAME( boost::uint_t<16>::fast, boost::crc_detail::uint_t<16>::fast );
    BOOST_TEST_TRAIT_SAME( boost::uint_t<32>::fast, boost::crc_detail::uint_t<32>::fast );
    BOOST_TEST_TRAIT_SAME( boost::uint_t<64>::fast, boost::crc_detail::uint_t<64>::fast );

    return boost::report_errors();
}
