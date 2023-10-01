//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include "test_unit/test_stream.hpp"

#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/any_stream.hpp>

#include <boost/asio/buffer.hpp>
#include <boost/asio/compose.hpp>
#include <boost/asio/coroutine.hpp>
#include <boost/asio/error.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/post.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <set>
#include <vector>

#include "test_common/buffer_concat.hpp"
#include "test_common/tracker_executor.hpp"

using namespace boost::mysql::test;
using boost::mysql::error_code;

static boost::asio::io_context ctx;

std::size_t boost::mysql::test::test_stream::get_size_to_read(std::size_t buffer_size) const
{
    auto it = read_break_offsets_.upper_bound(num_bytes_read_);
    std::size_t max_bytes_by_break = it == read_break_offsets_.end() ? std::size_t(-1)
                                                                     : *it - num_bytes_read_;
    return (std::min)({num_unread_bytes(), buffer_size, max_bytes_by_break});
}

std::size_t boost::mysql::test::test_stream::do_read(asio::mutable_buffer buff, error_code& ec)
{
    // Fail count
    error_code err = fail_count_.maybe_fail();
    if (err)
    {
        ec = err;
        return 0;
    }

    // If the user requested some bytes but we don't have any,
    // fail. In the real world, the stream would block until more
    // bytes are received, but this is a test, and this condition
    // indicates an error.
    if (num_unread_bytes() == 0 && buff.size() != 0)
    {
        ec = boost::asio::error::eof;
        return 0;
    }

    // Actually read
    std::size_t bytes_to_transfer = get_size_to_read(buff.size());
    if (bytes_to_transfer)
    {
        std::memcpy(buff.data(), bytes_to_read_.data() + num_bytes_read_, bytes_to_transfer);
        num_bytes_read_ += bytes_to_transfer;
    }

    // Clear errors
    ec = error_code();

    return bytes_to_transfer;
}

std::size_t boost::mysql::test::test_stream::do_write(asio::const_buffer buff, error_code& ec)
{
    // Fail count
    error_code err = fail_count_.maybe_fail();
    if (err)
    {
        ec = err;
        return 0;
    }

    // Actually write
    std::size_t num_bytes_to_transfer = (std::min)(buff.size(), write_break_size_);
    concat(bytes_written_, buff.data(), num_bytes_to_transfer);

    // Clear errors
    ec = error_code();

    return num_bytes_to_transfer;
}

struct boost::mysql::test::test_stream::read_op : boost::asio::coroutine
{
    test_stream& stream_;
    asio::mutable_buffer buff_;

    read_op(test_stream& stream, asio::mutable_buffer buff) noexcept : stream_(stream), buff_(buff){};

    template <class Self>
    void operator()(Self& self)
    {
        BOOST_ASIO_CORO_REENTER(*this)
        {
            BOOST_ASIO_CORO_YIELD boost::asio::post(stream_.get_executor(), std::move(self));
            {
                error_code err;
                std::size_t bytes_read = stream_.do_read(buff_, err);
                self.complete(err, bytes_read);
            }
        }
    }
};

struct boost::mysql::test::test_stream::write_op : boost::asio::coroutine
{
    test_stream& stream_;
    asio::const_buffer buff_;

    write_op(test_stream& stream, asio::const_buffer buff) noexcept : stream_(stream), buff_(buff){};

    template <class Self>
    void operator()(Self& self)
    {
        BOOST_ASIO_CORO_REENTER(*this)
        {
            BOOST_ASIO_CORO_YIELD boost::asio::post(stream_.get_executor(), std::move(self));
            {
                error_code err;
                std::size_t bytes_written = stream_.do_write(buff_, err);
                self.complete(err, bytes_written);
            }
        }
    }
};

boost::mysql::test::test_stream::executor_type boost::mysql::test::test_stream::get_executor()
{
    return create_tracker_executor(ctx.get_executor(), &executor_info_);
}

// Reading
std::size_t boost::mysql::test::test_stream::read_some(asio::mutable_buffer buff, error_code& ec)
{
    return do_read(buff, ec);
}
void boost::mysql::test::test_stream::async_read_some(
    asio::mutable_buffer buff,
    asio::any_completion_handler<void(error_code, std::size_t)> handler
)
{
    boost::asio::async_compose<
        asio::any_completion_handler<void(error_code, std::size_t)>,
        void(error_code, std::size_t)>(read_op(*this, buff), handler, get_executor());
}

// Writing
std::size_t boost::mysql::test::test_stream::write_some(boost::asio::const_buffer buff, error_code& ec)
{
    return do_write(buff, ec);
}

void boost::mysql::test::test_stream::async_write_some(
    boost::asio::const_buffer buff,
    asio::any_completion_handler<void(error_code, std::size_t)> handler
)
{
    boost::asio::async_compose<
        asio::any_completion_handler<void(error_code, std::size_t)>,
        void(error_code, std::size_t)>(write_op(*this, buff), handler, get_executor());
}

test_stream& boost::mysql::test::test_stream::add_bytes(span<const std::uint8_t> bytes)
{
    concat(bytes_to_read_, bytes.data(), bytes.size());
    return *this;
}

test_stream& boost::mysql::test::test_stream::add_break(std::size_t byte_num)
{
    BOOST_ASSERT(byte_num <= bytes_to_read_.size());
    read_break_offsets_.insert(byte_num);
    return *this;
}

template class boost::mysql::detail::any_stream_impl<boost::mysql::test::test_stream>;
