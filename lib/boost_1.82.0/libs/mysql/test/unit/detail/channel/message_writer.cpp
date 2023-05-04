//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>

#include <boost/mysql/detail/channel/message_writer.hpp>

#include <boost/asio/bind_executor.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/test/unit_test.hpp>

#include <vector>

#include "assert_buffer_equals.hpp"
#include "buffer_concat.hpp"
#include "create_message.hpp"
#include "test_stream.hpp"

using boost::asio::buffer;
using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::detail::message_writer;
using boost::mysql::test::concat_copy;
using boost::mysql::test::create_message;
using boost::mysql::test::fail_count;
using boost::mysql::test::test_stream;

namespace {

// Machinery to be able to cover the sync and async
// functions with the same test code
class writer_fns
{
public:
    virtual ~writer_fns() {}
    virtual void write(
        message_writer& writer,
        test_stream& stream,
        boost::asio::const_buffer buffer,
        std::uint8_t& seqnum,
        error_code& ec
    ) = 0;
    virtual const char* name() const noexcept = 0;
};

class sync_writer_fns : public writer_fns
{
public:
    void write(
        message_writer& writer,
        test_stream& stream,
        boost::asio::const_buffer buffer,
        std::uint8_t& seqnum,
        error_code& ec
    ) final override
    {
        writer.write(stream, buffer, seqnum, ec);
    }
    const char* name() const noexcept final override { return "sync"; };
};

class async_writer_fns : public writer_fns
{
public:
    void write(
        message_writer& writer,
        test_stream& stream,
        boost::asio::const_buffer buffer,
        std::uint8_t& seqnum,
        error_code& err
    ) final override
    {
        boost::asio::io_context ctx;
        writer.async_write(
            stream,
            buffer,
            seqnum,
            boost::asio::bind_executor(ctx.get_executor(), [&](error_code ec) { err = ec; })
        );
        ctx.run();
    }
    const char* name() const noexcept final override { return "async"; };
};

sync_writer_fns sync;
async_writer_fns async;
writer_fns* all_reader_fns[] = {&sync, &async};

BOOST_AUTO_TEST_SUITE(test_message_writer)

BOOST_AUTO_TEST_CASE(success)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_writer writer(8);
            test_stream stream;
            std::vector<std::uint8_t> msg{0x01, 0x02, 0x03};
            std::uint8_t seqnum = 4;
            auto expected_msg = create_message(seqnum, msg);
            error_code err(client_errc::wrong_num_params);

            fns->write(writer, stream, buffer(msg), seqnum, err);

            BOOST_TEST(err == error_code());
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffer(stream.bytes_written()), buffer(expected_msg));
            BOOST_TEST(seqnum == 5);
        }
    }
}

BOOST_AUTO_TEST_CASE(empty_message)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_writer writer(8);
            test_stream stream;
            std::vector<std::uint8_t> msg{};
            std::uint8_t seqnum = 4;
            auto expected_msg = create_message(seqnum, msg);
            error_code err(client_errc::wrong_num_params);

            fns->write(writer, stream, buffer(msg), seqnum, err);

            BOOST_TEST(err == error_code());
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffer(stream.bytes_written()), buffer(expected_msg));
            BOOST_TEST(seqnum == 5);
        }
    }
}

BOOST_AUTO_TEST_CASE(short_writes)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_writer writer(8);
            test_stream stream;
            stream.set_write_break_size(2);  // 2 byte reads
            std::vector<std::uint8_t> msg{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07};
            std::uint8_t seqnum = 4;
            auto expected_msg = create_message(seqnum, msg);
            error_code err(client_errc::wrong_num_params);

            fns->write(writer, stream, buffer(msg), seqnum, err);

            BOOST_TEST(err == error_code());
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffer(stream.bytes_written()), buffer(expected_msg));
            BOOST_TEST(seqnum == 5);
        }
    }
}

BOOST_AUTO_TEST_CASE(multi_frame)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_writer writer(8);
            test_stream stream;
            std::vector<std::uint8_t> msg_frame_1{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
            std::vector<std::uint8_t> msg_frame_2{0x11, 0x12};
            auto msg = concat_copy(msg_frame_1, msg_frame_2);
            std::uint8_t seqnum = 4;
            auto expected_msg = create_message(seqnum, msg_frame_1, seqnum + 1, msg_frame_2);
            error_code err(client_errc::wrong_num_params);

            fns->write(writer, stream, buffer(msg), seqnum, err);

            BOOST_TEST(err == error_code());
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffer(stream.bytes_written()), buffer(expected_msg));
            BOOST_TEST(seqnum == 6);
        }
    }
}

BOOST_AUTO_TEST_CASE(error)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_writer writer(8);
            test_stream stream(fail_count(0, error_code(client_errc::server_unsupported)));
            std::vector<std::uint8_t> msg{0x01, 0x02, 0x03};
            std::uint8_t seqnum = 4;
            auto expected_msg = create_message(seqnum, msg);
            error_code err(client_errc::wrong_num_params);

            fns->write(writer, stream, buffer(msg), seqnum, err);

            BOOST_TEST(err == error_code(client_errc::server_unsupported));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace