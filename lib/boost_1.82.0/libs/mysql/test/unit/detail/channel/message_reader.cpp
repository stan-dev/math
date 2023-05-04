//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/channel/message_reader.hpp>

#include <boost/asio/bind_executor.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/system_executor.hpp>
#include <boost/test/unit_test.hpp>

#include <system_error>

#include "assert_buffer_equals.hpp"
#include "buffer_concat.hpp"
#include "create_message.hpp"
#include "test_stream.hpp"

using boost::asio::buffer;
using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::detail::message_reader;
using boost::mysql::test::create_message;
using boost::mysql::test::fail_count;
using boost::mysql::test::test_stream;

namespace {

// Machinery to be able to cover the sync and async
// functions with the same test code
class reader_fns
{
public:
    virtual ~reader_fns() {}
    virtual void read_some(
        message_reader&,
        test_stream&,
        error_code&,
        bool keep_messages = false
    ) = 0;
    virtual boost::asio::const_buffer read_one(
        message_reader&,
        test_stream&,
        std::uint8_t& seqnum,
        error_code& ec,
        bool keep_messages = false
    ) = 0;
    virtual const char* name() const noexcept = 0;
};

class sync_reader_fns : public reader_fns
{
public:
    void read_some(
        message_reader& reader,
        test_stream& stream,
        error_code& err,
        bool keep_messages = false
    ) final override
    {
        reader.read_some(stream, err, keep_messages);
    }
    boost::asio::const_buffer read_one(
        message_reader& reader,
        test_stream& stream,
        std::uint8_t& seqnum,
        error_code& ec,
        bool keep_messages = false
    ) final override
    {
        return reader.read_one(stream, seqnum, ec, keep_messages);
    }
    const char* name() const noexcept final override { return "sync"; };
};

class async_reader_fns : public reader_fns
{
public:
    void read_some(
        message_reader& reader,
        test_stream& stream,
        error_code& err,
        bool keep_messages = false
    ) final override
    {
        boost::asio::io_context ctx;
        reader.async_read_some(
            stream,
            boost::asio::bind_executor(ctx.get_executor(), [&](error_code ec) { err = ec; }),
            keep_messages
        );
        ctx.run();
    }
    boost::asio::const_buffer read_one(
        message_reader& reader,
        test_stream& stream,
        std::uint8_t& seqnum,
        error_code& err,
        bool keep_messages = false
    ) final override
    {
        boost::asio::io_context ctx;
        boost::asio::const_buffer res;
        reader.async_read_one(
            stream,
            seqnum,
            boost::asio::bind_executor(
                ctx.get_executor(),
                [&](error_code ec, boost::asio::const_buffer b) {
                    err = ec;
                    res = b;
                }
            ),
            keep_messages
        );
        ctx.run();
        return res;
    }
    const char* name() const noexcept final override { return "async"; };
};

sync_reader_fns sync;
async_reader_fns async;
reader_fns* all_reader_fns[] = {&sync, &async};

BOOST_AUTO_TEST_SUITE(test_message_reader)

BOOST_AUTO_TEST_SUITE(read_some)

BOOST_AUTO_TEST_CASE(message_fits_in_buffer)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
            test_stream stream(create_message(seqnum, msg_body));
            error_code err(client_errc::server_unsupported);

            // Doesn't have a message initially
            BOOST_TEST(!reader.has_message());

            // Read succesfully
            fns->read_some(reader, stream, err);
            BOOST_TEST(err == error_code());
            BOOST_REQUIRE(reader.has_message());
            BOOST_TEST(stream.num_unread_bytes() == 0u);

            // Get next message and validate it
            auto msg = reader.get_next_message(seqnum, err);
            BOOST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg_body));

            // There isn't another message
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(fragmented_message_fits_in_buffer)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
            test_stream stream(test_stream::read_behavior(
                create_message(seqnum, msg_body),
                {3, 5}  // break the message at bytes 3 and 5
            ));
            error_code err(client_errc::server_unsupported);

            // Doesn't have a message initially
            BOOST_TEST(!reader.has_message());

            // Read succesfully
            fns->read_some(reader, stream, err);
            BOOST_TEST(err == error_code());
            BOOST_REQUIRE(reader.has_message());
            BOOST_TEST(stream.num_unread_bytes() == 0u);

            // Get next message and validate it
            auto msg = reader.get_next_message(seqnum, err);
            BOOST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg_body));

            // There isn't another message
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(message_doesnt_fit_in_buffer)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(0);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
            test_stream stream(create_message(seqnum, msg_body));
            error_code err(client_errc::server_unsupported);

            // Doesn't have a message initially
            BOOST_TEST(!reader.has_message());

            // Read succesfully
            fns->read_some(reader, stream, err);
            BOOST_TEST(err == error_code());
            BOOST_REQUIRE(reader.has_message());
            BOOST_TEST(reader.buffer().size() >= msg_body.size());
            BOOST_TEST(stream.num_unread_bytes() == 0u);

            // Get next message and validate it
            auto msg = reader.get_next_message(seqnum, err);
            BOOST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg_body));

            // There isn't another message
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(two_messages)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            std::uint8_t seqnum1 = 2;
            std::uint8_t seqnum2 = 5;
            std::vector<std::uint8_t> msg1_body{0x01, 0x02, 0x03};
            std::vector<std::uint8_t> msg2_body{0x05, 0x06, 0x07, 0x08};
            test_stream stream(create_message(seqnum1, msg1_body, seqnum2, msg2_body));
            error_code err(client_errc::server_unsupported);

            // Doesn't have a message initially
            BOOST_TEST(!reader.has_message());

            // Read succesfully
            fns->read_some(reader, stream, err);
            BOOST_TEST(err == error_code());
            BOOST_REQUIRE(reader.has_message());
            BOOST_TEST(stream.num_unread_bytes() == 0u);

            // Get next message and validate it
            auto msg = reader.get_next_message(seqnum1, err);
            BOOST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum1 == 3);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg1_body));

            // Reading again does nothing
            fns->read_some(reader, stream, err);
            BOOST_TEST(err == error_code());

            // Get the 2nd message and validate it
            msg = reader.get_next_message(seqnum2, err);
            BOOST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum2 == 6u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg2_body));

            // There isn't another message
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(previous_message_keep_messages_false)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            std::uint8_t seqnum1 = 2;
            std::uint8_t seqnum2 = 5;
            std::vector<std::uint8_t> msg1_body{0x01, 0x02, 0x03};
            std::vector<std::uint8_t> msg2_body{0x05, 0x06, 0x07};
            test_stream stream;
            stream.add_message(create_message(seqnum1, msg1_body));
            stream.add_message(create_message(seqnum2, msg2_body));
            error_code err(client_errc::server_unsupported);

            // Read and get 1st message
            fns->read_some(reader, stream, err, false);
            BOOST_TEST(err == error_code());
            BOOST_REQUIRE(reader.has_message());
            auto msg1 = reader.get_next_message(seqnum1, err);
            BOOST_REQUIRE(err == error_code());

            // Read and get 2nd message
            fns->read_some(reader, stream, err, false);
            BOOST_TEST(err == error_code());
            BOOST_REQUIRE(reader.has_message());
            auto msg2 = reader.get_next_message(seqnum2, err);
            BOOST_REQUIRE(err == error_code());
            BOOST_TEST(stream.num_unread_bytes() == 0u);
            BOOST_TEST(seqnum2 == 6);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg2, buffer(msg2_body));
            BOOST_TEST(!reader.has_message());

            // The 2nd message is located where the 1st message was
            BOOST_TEST(msg1.data() == msg2.data());
        }
    }
}

BOOST_AUTO_TEST_CASE(previous_message_keep_messages_true)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            std::uint8_t seqnum1 = 2;
            std::uint8_t seqnum2 = 5;
            std::vector<std::uint8_t> msg1_body{0x01, 0x02, 0x03};
            std::vector<std::uint8_t> msg2_body{0x05, 0x06, 0x07};
            test_stream stream;
            stream.add_message(create_message(seqnum1, msg1_body));
            stream.add_message(create_message(seqnum2, msg2_body));
            error_code err(client_errc::server_unsupported);

            // Read and get 1st message
            fns->read_some(reader, stream, err, true);
            BOOST_TEST(err == error_code());
            BOOST_REQUIRE(reader.has_message());
            auto msg1 = reader.get_next_message(seqnum1, err);
            BOOST_REQUIRE(err == error_code());

            // Read and get 2nd message
            fns->read_some(reader, stream, err, true);
            BOOST_TEST(err == error_code());
            BOOST_REQUIRE(reader.has_message());
            auto msg2 = reader.get_next_message(seqnum2, err);
            BOOST_REQUIRE(err == error_code());
            BOOST_TEST(stream.num_unread_bytes() == 0u);
            BOOST_TEST(seqnum2 == 6u);
            BOOST_TEST(!reader.has_message());

            // Both messages are valid
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg1, buffer(msg1_body));
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg2, buffer(msg2_body));
        }
    }
}

BOOST_AUTO_TEST_CASE(error)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            test_stream stream(fail_count(0, client_errc::wrong_num_params));
            error_code err(client_errc::server_unsupported);

            // Read with error
            fns->read_some(reader, stream, err);
            BOOST_TEST(!reader.has_message());
            BOOST_TEST(err == error_code(client_errc::wrong_num_params));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

// Cases specific to get_next_message
BOOST_AUTO_TEST_SUITE(get_next_message)

BOOST_AUTO_TEST_CASE(multiframe_message)
{
    message_reader reader(512, 8);
    std::uint8_t seqnum = 2;
    test_stream stream(
        create_message(seqnum, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08}, 3, {0x09, 0x0a})
    );
    std::vector<std::uint8_t>
        expected_msg{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a};
    error_code err(client_errc::server_unsupported);

    // Read succesfully
    reader.read_some(stream, err);
    BOOST_TEST(err == error_code());
    BOOST_REQUIRE(reader.has_message());
    BOOST_TEST(stream.num_unread_bytes() == 0u);

    // Get next message and validate it
    auto msg = reader.get_next_message(seqnum, err);
    BOOST_REQUIRE(err == error_code());
    BOOST_TEST(seqnum == 4u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(expected_msg));

    // There isn't another message
    BOOST_TEST(!reader.has_message());
}

BOOST_AUTO_TEST_CASE(seqnum_overflow)
{
    message_reader reader(512);
    std::uint8_t seqnum = 0xff;
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
    test_stream stream(create_message(seqnum, msg_body));
    error_code err(client_errc::server_unsupported);

    // Read succesfully
    reader.read_some(stream, err);
    BOOST_TEST(err == error_code());
    BOOST_REQUIRE(reader.has_message());
    BOOST_TEST(stream.num_unread_bytes() == 0u);

    // Get next message and validate it
    auto msg = reader.get_next_message(seqnum, err);
    BOOST_REQUIRE(err == error_code());
    BOOST_TEST(seqnum == 0u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg_body));
}

BOOST_AUTO_TEST_CASE(error_passed_seqnum_mismatch)
{
    message_reader reader(512);
    test_stream stream(create_message(2, {0x01, 0x02, 0x03}));
    error_code err(client_errc::server_unsupported);

    // Read succesfully
    reader.read_some(stream, err);
    BOOST_TEST(err == error_code());
    BOOST_REQUIRE(reader.has_message());
    BOOST_TEST(stream.num_unread_bytes() == 0u);

    // Passed-in seqnum is invalid
    std::uint8_t bad_seqnum = 0;
    reader.get_next_message(bad_seqnum, err);
    BOOST_TEST(err == make_error_code(client_errc::sequence_number_mismatch));
    BOOST_TEST(bad_seqnum == 0u);
}

BOOST_AUTO_TEST_CASE(error_intermediate_frame_seqnum_mismatch)
{
    message_reader reader(512, 8);  // frames are broken each 8 bytes
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
    std::uint8_t seqnum = 2;
    test_stream stream(create_message(
        seqnum,
        {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08},
        4,
        {0x11, 0x12, 0x13, 0x14}  // the right seqnum would be 3
    ));
    error_code err(client_errc::server_unsupported);

    // Read succesfully
    reader.read_some(stream, err);
    BOOST_TEST(err == error_code());
    BOOST_REQUIRE(reader.has_message());
    BOOST_TEST(stream.num_unread_bytes() == 0u);

    // The read frame has a mismatched seqnum
    reader.get_next_message(seqnum, err);
    BOOST_TEST(err == make_error_code(client_errc::sequence_number_mismatch));
    BOOST_TEST(seqnum == 2u);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(read_one)

BOOST_AUTO_TEST_CASE(success)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
            test_stream stream(create_message(seqnum, msg_body));
            error_code err(client_errc::server_unsupported);

            // Read succesfully
            auto msg = fns->read_one(reader, stream, seqnum, err);
            BOOST_TEST(err == error_code());
            BOOST_TEST(seqnum == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg_body));

            BOOST_TEST(stream.num_unread_bytes() == 0u);
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(cached_message)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            std::uint8_t seqnum1 = 2;
            std::uint8_t seqnum2 = 8;
            std::vector<std::uint8_t> msg1_body{0x01, 0x02, 0x03};
            std::vector<std::uint8_t> msg2_body{0x04, 0x05};
            test_stream stream(create_message(seqnum1, msg1_body, seqnum2, msg2_body));
            error_code err(client_errc::server_unsupported);

            // Read succesfully
            auto msg = fns->read_one(reader, stream, seqnum1, err);
            BOOST_TEST(err == error_code());
            BOOST_TEST(seqnum1 == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg1_body));
            BOOST_TEST(reader.has_message());

            // Read again
            msg = fns->read_one(reader, stream, seqnum2, err);
            BOOST_TEST(err == error_code());
            BOOST_TEST(seqnum2 == 9u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, buffer(msg2_body));
            BOOST_TEST(stream.num_unread_bytes() == 0u);
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(error_in_read)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
            test_stream stream(
                create_message(seqnum, msg_body),
                fail_count(0, error_code(client_errc::wrong_num_params))
            );
            error_code err(client_errc::server_unsupported);

            fns->read_one(reader, stream, seqnum, err);
            BOOST_TEST(err == error_code(client_errc::wrong_num_params));
            BOOST_TEST(seqnum == 2u);
        }
    }
}

BOOST_AUTO_TEST_CASE(seqnum_mismatch)
{
    for (auto* fns : all_reader_fns)
    {
        BOOST_TEST_CONTEXT(fns->name())
        {
            message_reader reader(512);
            test_stream stream(create_message(2, {0x01, 0x02, 0x03}));
            error_code err(client_errc::server_unsupported);

            std::uint8_t bad_seqnum = 42;
            fns->read_one(reader, stream, bad_seqnum, err);
            BOOST_TEST(err == error_code(client_errc::sequence_number_mismatch));
            BOOST_TEST(bad_seqnum == 42u);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace