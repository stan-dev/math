//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/impl/internal/protocol/static_buffer.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"

using namespace boost::mysql::test;
using boost::mysql::detail::static_buffer;

BOOST_AUTO_TEST_SUITE(test_static_buffer)

struct fixture
{
    static constexpr std::size_t max_size_value = 32;
    using string_type = static_buffer<max_size_value>;

    static std::vector<std::uint8_t> original_midsize() { return {0x01, 0x02, 0x03}; }
    static std::vector<std::uint8_t> original_maxsize()
    {
        return std::vector<std::uint8_t>(max_size_value, 0xde);
    }

    std::vector<std::uint8_t> midsize = original_midsize();
    std::vector<std::uint8_t> maxsize = original_maxsize();

    void wipe_midsize() { midsize = {0xa, 0xb, 0xc}; }
    void wipe_maxsize() { maxsize = std::vector<std::uint8_t>(max_size_value, 0xaa); }
};

// Default ctor.
BOOST_FIXTURE_TEST_CASE(default_constructor, fixture)
{
    string_type v;
    BOOST_TEST(v.to_span().size() == 0u);
}

// clear
BOOST_FIXTURE_TEST_CASE(clear_empty, fixture)
{
    string_type v;
    v.clear();
    BOOST_TEST(v.to_span().size() == 0u);
}

BOOST_FIXTURE_TEST_CASE(clear_not_empty, fixture)
{
    string_type v;
    const std::array<std::uint8_t, 5> data{
        {0, 1, 2, 3, 4}
    };
    v.append(data.data(), data.size());
    BOOST_TEST(v.to_span().size() == 5u);
    v.clear();
    BOOST_TEST(v.to_span().size() == 0u);
}

// append
BOOST_FIXTURE_TEST_CASE(append_from_empty_to_empty, fixture)
{
    string_type v;
    v.append(midsize.data(), 0);
    wipe_midsize();
    BOOST_TEST(v.to_span().size() == 0u);
}

BOOST_FIXTURE_TEST_CASE(append_from_empty_to_midsize, fixture)
{
    string_type v;
    v.append(midsize.data(), midsize.size());
    wipe_midsize();
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(v.to_span(), original_midsize());
}

BOOST_FIXTURE_TEST_CASE(append_from_empty_to_maxsize, fixture)
{
    string_type v;
    v.append(maxsize.data(), maxsize.size());
    wipe_maxsize();
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(v.to_span(), original_maxsize());
}

BOOST_FIXTURE_TEST_CASE(append_from_midsize_to_midsize, fixture)
{
    // Initial
    string_type v;
    std::vector<std::uint8_t> initial{2, 2, 2};
    v.append(initial.data(), initial.size());

    // Append more data
    v.append(midsize.data(), midsize.size());
    wipe_midsize();  // Check it was actually copied

    // Verify
    auto expected = concat(initial, original_midsize());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(v.to_span(), expected);
}

BOOST_FIXTURE_TEST_CASE(append_from_midsize_to_maxsize, fixture)
{
    // Initial
    string_type v;
    v.append(midsize.data(), midsize.size());

    // Append
    std::vector<std::uint8_t> newbuff(max_size_value - midsize.size(), 1);
    v.append(newbuff.data(), newbuff.size());
    wipe_midsize();  // Verify that we actually copied the data

    // Verify
    auto expected = concat(original_midsize(), newbuff);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(v.to_span(), expected);
}

BOOST_AUTO_TEST_SUITE_END()
