//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_TEST_PROTOCOL_SERIALIZATION_TEST_HPP
#define BOOST_MYSQL_TEST_UNIT_TEST_PROTOCOL_SERIALIZATION_TEST_HPP

#include <boost/mysql/impl/internal/protocol/impl/deserialization_context.hpp>
#include <boost/mysql/impl/internal/protocol/impl/serialization_context.hpp>

#include <boost/asio/buffer.hpp>
#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstring>
#include <initializer_list>

#include "test_common/assert_buffer_equals.hpp"
#include "test_unit/serialize_to_vector.hpp"

namespace boost {
namespace mysql {
namespace test {

// A special buffer for deserialization tests. Allocates the exact size of the serialized message (contrary to
// std::vector), making it easier for sanitizers to detect overruns
class deserialization_buffer
{
    std::size_t size_;
    std::unique_ptr<std::uint8_t[]> data_;

public:
    explicit deserialization_buffer(std::size_t size) : size_(size), data_(new std::uint8_t[size]) {}
    deserialization_buffer(std::size_t size, std::uint8_t value) : deserialization_buffer(size)
    {
        std::memset(data(), value, size);
    }
    deserialization_buffer(span<const std::uint8_t> data) : size_(data.size())
    {
        if (!data.empty())
        {
            data_.reset(new std::uint8_t[data.size()]);
            std::memcpy(data_.get(), data.data(), data.size());
        }
    }
    deserialization_buffer(const std::vector<std::uint8_t>& data)
        : deserialization_buffer(span<const std::uint8_t>(data))
    {
    }
    deserialization_buffer(std::initializer_list<std::uint8_t> data)
        : size_(data.size()), data_(new std::uint8_t[data.size()])
    {
        std::copy(data.begin(), data.end(), data_.get());
    }
    span<const std::uint8_t> to_span() const noexcept { return {data_.get(), size_}; }
    operator span<const std::uint8_t>() const noexcept { return to_span(); }
    std::uint8_t* data() noexcept { return data_.get(); }
    const std::uint8_t* data() const noexcept { return data_.get(); }
    std::size_t size() const noexcept { return size_; }
};

template <class T>
void do_serialize_test(T value, span<const std::uint8_t> expected)
{
    // Serialize
    auto actual = serialize_to_vector([&](detail::serialization_context& ctx) { value.serialize(ctx); });

    // Check
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(expected, actual);
}

template <class T>
void do_deserialize_test(T value, span<const std::uint8_t> serialized)
{
    deserialization_buffer buffer(serialized);
    auto expected_first = buffer.data() + buffer.size();
    detail::deserialization_context ctx(buffer);
    T actual{};
    detail::deserialize_errc err = actual.deserialize(ctx);

    // No error
    BOOST_TEST(err == detail::deserialize_errc::ok);

    // Iterator advanced
    BOOST_TEST(ctx.first() == expected_first);

    // Actual value
    BOOST_TEST(actual == value);
}

template <class T>
void do_deserialize_extra_space_test(T value, span<const std::uint8_t> serialized)
{
    // Create a buffer with extra data
    deserialization_buffer buffer(serialized.size() + 1);
    std::memcpy(buffer.data(), serialized.data(), serialized.size());
    buffer.data()[serialized.size()] = 0xff;

    // Deserialize
    detail::deserialization_context ctx(buffer);
    T actual{};
    detail::deserialize_errc err = actual.deserialize(ctx);

    // No error
    BOOST_TEST(err == detail::deserialize_errc::ok);

    // Iterator advanced
    BOOST_TEST(ctx.first() == buffer.data() + serialized.size());

    // Actual value
    BOOST_TEST(actual == value);
}

template <class T>
void do_deserialize_not_enough_space_test(span<const std::uint8_t> serialized)
{
    // Create a new buffer with one less byte
    deserialization_buffer buffer(serialized.subspan(0, serialized.size() - 1));
    detail::deserialization_context ctx(buffer.to_span());

    T value{};
    detail::deserialize_errc err = value.deserialize(ctx);
    BOOST_TEST(err == detail::deserialize_errc::incomplete_message);
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
