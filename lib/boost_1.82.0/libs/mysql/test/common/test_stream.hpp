//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_TEST_STREAM_HPP
#define BOOST_MYSQL_TEST_COMMON_TEST_STREAM_HPP

#include <boost/mysql/error_code.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/system_executor.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <set>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

// Inspired by Beast's fail count
class fail_count
{
    std::size_t fail_after_;
    std::size_t num_calls_{0};
    error_code err_;

public:
    static constexpr std::size_t never_fail = std::size_t(-1);
    explicit fail_count(
        std::size_t fail_after = never_fail,
        error_code err = make_error_code(std::errc::io_error)
    ) noexcept
        : fail_after_(fail_after), err_(err)
    {
    }
    error_code maybe_fail() noexcept { return ++num_calls_ >= fail_after_ ? err_ : error_code(); }
};

class test_stream
{
public:
    // Executors
    using executor_type = boost::asio::any_io_executor;
    executor_type get_executor() noexcept { return executor_; }

    using lowest_layer_type = test_stream;
    lowest_layer_type& lowest_layer() { return *this; }

    struct read_behavior
    {
        std::vector<std::uint8_t> bytes_to_read;
        std::set<std::size_t> break_offsets;

        read_behavior() = default;
        read_behavior(std::vector<std::uint8_t> bytes, std::set<std::size_t> offsets = {})
            : bytes_to_read(std::move(bytes)), break_offsets(std::move(offsets))
        {
            assert(break_offsets.empty() || *break_offsets.rbegin() < bytes_to_read.size());
        }
    };

    // Constructors
    inline test_stream(
        fail_count fc = fail_count(),
        executor_type ex = boost::asio::system_executor()
    );

    inline test_stream(
        std::vector<std::uint8_t> bytes_to_read,
        fail_count fc = fail_count(),
        executor_type ex = boost::asio::system_executor()
    );

    inline test_stream(
        read_behavior b,
        fail_count fc = fail_count(),
        executor_type ex = boost::asio::system_executor()
    );

    // Setting test behavior
    inline void add_message(const std::vector<std::uint8_t>& bytes, bool separate_reads = true);
    inline void set_read_behavior(read_behavior b);
    void set_write_break_size(std::size_t size) noexcept { write_break_size_ = size; }
    void set_fail_count(const fail_count& fc) noexcept { fail_count_ = fc; }
    void set_executor(executor_type ex) { executor_ = ex; }

    // Getting test results
    std::size_t num_bytes_read() const noexcept { return num_bytes_read_; }
    std::size_t num_unread_bytes() const noexcept
    {
        return bytes_to_read_.size() - num_bytes_read_;
    }
    const std::vector<std::uint8_t>& bytes_written() const noexcept { return bytes_written_; }

    // Stream operations
    template <class MutableBufferSequence>
    std::size_t read_some(const MutableBufferSequence& buffers, error_code& ec)
    {
        return do_read(buffers, ec);
    }

    template <class ConstBufferSequence>
    std::size_t write_some(const ConstBufferSequence& buffers, error_code& ec)
    {
        return do_write(buffers, ec);
    }

    template <
        class MutableBufferSequence,
        BOOST_ASIO_COMPLETION_TOKEN_FOR(void(::boost::mysql::error_code, std::size_t))
            CompletionToken>
    BOOST_ASIO_INITFN_AUTO_RESULT_TYPE(CompletionToken, void(error_code, std::size_t))
    async_read_some(const MutableBufferSequence& buffers, CompletionToken&& token);

    template <
        class ConstBufferSequence,
        BOOST_ASIO_COMPLETION_TOKEN_FOR(void(::boost::mysql::error_code, std::size_t))
            CompletionToken>
    BOOST_ASIO_INITFN_AUTO_RESULT_TYPE(CompletionToken, void(error_code, std::size_t))
    async_write_some(ConstBufferSequence const& buffers, CompletionToken&& token);

private:
    std::vector<std::uint8_t> bytes_to_read_;
    std::set<std::size_t> read_break_offsets_;
    std::size_t num_bytes_read_{0};
    std::vector<std::uint8_t> bytes_written_;
    fail_count fail_count_;
    std::size_t write_break_size_{1024};  // max number of bytes to be written in each write_some
    executor_type executor_;

    template <class MutableBufferSequence>
    struct read_op;
    template <class ConstBufferSequence>
    struct write_op;

    inline std::size_t get_size_to_read(std::size_t buffer_size) const;

    template <class MutableBufferSequence>
    std::size_t do_read(const MutableBufferSequence& buffers, error_code& ec);

    template <class ConstBufferSequence>
    std::size_t do_write(const ConstBufferSequence& buffers, error_code& ec);
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#include "impl/test_stream.hpp"

#endif