//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/any_stream.hpp>
#include <boost/mysql/detail/any_stream_impl.hpp>

#include <boost/mysql/impl/internal/channel/message_reader.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::span;
using boost::mysql::client_errc;
using boost::mysql::error_code;

BOOST_AUTO_TEST_SUITE(test_message_reader)

using netfun_maker_some = netfun_maker_fn<void, any_stream&, message_reader&>;
using netfun_maker_one = netfun_maker_fn<
    span<const std::uint8_t>,
    any_stream&,
    message_reader&,
    std::uint8_t&>;

struct
{
    netfun_maker_some::signature read_some;
    netfun_maker_one::signature read_one;
    const char* name;
} all_fns[] = {
    {netfun_maker_some::sync_errc_noerrinfo(&read_some_messages),
     netfun_maker_one::sync_errc_noerrinfo(&read_one_message),
     "sync" },
    {netfun_maker_some::async_noerrinfo(&async_read_some_messages),
     netfun_maker_one::async_noerrinfo(&async_read_one_message),
     "async"},
};

struct fixture
{
    test_any_stream stream;

    test_stream& inner_stream() noexcept { return cast<test_stream>(stream); }
};

BOOST_AUTO_TEST_SUITE(read_some)

BOOST_AUTO_TEST_CASE(message_fits_in_buffer)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
            fix.inner_stream().add_bytes(create_frame(seqnum, msg_body));
            error_code err = client_errc::unknown_auth_plugin;

            // Doesn't have a message initially
            BOOST_TEST(!reader.has_message());

            // Read succesfully
            fn.read_some(fix.stream, reader).validate_no_error();
            BOOST_TEST_REQUIRE(reader.has_message());
            BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);

            // Get next message and validate it
            auto msg = reader.get_next_message(seqnum, err);
            BOOST_TEST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, msg_body);

            // There isn't another message
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(fragmented_message_fits_in_buffer)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
            fix.inner_stream()
                .add_bytes(create_frame(seqnum, msg_body))
                .add_break(3)
                .add_break(5);  // break the message at bytes 3 and 5
            error_code err(client_errc::server_unsupported);

            // Doesn't have a message initially
            BOOST_TEST(!reader.has_message());

            // Read succesfully
            fn.read_some(fix.stream, reader).validate_no_error();
            BOOST_TEST_REQUIRE(reader.has_message());
            BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);

            // Get next message and validate it
            auto msg = reader.get_next_message(seqnum, err);
            BOOST_TEST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, msg_body);

            // There isn't another message
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(message_doesnt_fit_in_buffer)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(0);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
            fix.inner_stream().add_bytes(create_frame(seqnum, msg_body));
            error_code err(client_errc::server_unsupported);

            // Doesn't have a message initially
            BOOST_TEST(!reader.has_message());

            // Read succesfully
            fn.read_some(fix.stream, reader).validate_no_error();
            BOOST_TEST_REQUIRE(reader.has_message());
            BOOST_TEST(reader.buffer().size() >= msg_body.size());
            BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);

            // Get next message and validate it
            auto msg = reader.get_next_message(seqnum, err);
            BOOST_TEST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, msg_body);

            // There isn't another message
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(two_messages)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            std::uint8_t seqnum1 = 2;
            std::uint8_t seqnum2 = 5;
            std::vector<std::uint8_t> msg1_body{0x01, 0x02, 0x03};
            std::vector<std::uint8_t> msg2_body{0x05, 0x06, 0x07, 0x08};
            fix.inner_stream()
                .add_bytes(create_frame(seqnum1, msg1_body))
                .add_bytes(create_frame(seqnum2, msg2_body));
            error_code err(client_errc::server_unsupported);

            // Doesn't have a message initially
            BOOST_TEST(!reader.has_message());

            // Read succesfully
            fn.read_some(fix.stream, reader).validate_no_error();
            BOOST_TEST_REQUIRE(reader.has_message());
            BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);

            // Get next message and validate it
            auto msg1 = reader.get_next_message(seqnum1, err);
            BOOST_TEST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum1 == 3);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg1, msg1_body);

            // Get the 2nd message and validate it
            auto msg2 = reader.get_next_message(seqnum2, err);
            BOOST_TEST_REQUIRE(err == error_code());
            BOOST_TEST(seqnum2 == 6u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg2, msg2_body);

            // There isn't another message
            BOOST_TEST(!reader.has_message());

            // Reading again does nothing
            fn.read_some(fix.stream, reader);
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(previous_message)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            std::uint8_t seqnum1 = 2;
            std::uint8_t seqnum2 = 5;
            std::vector<std::uint8_t> msg1_body{0x01, 0x02, 0x03};
            std::vector<std::uint8_t> msg2_body{0x05, 0x06, 0x07};
            fix.inner_stream()
                .add_bytes(create_frame(seqnum1, msg1_body))
                .add_break()
                .add_bytes(create_frame(seqnum2, msg2_body));
            error_code err(client_errc::server_unsupported);

            // Read and get 1st message
            fn.read_some(fix.stream, reader).validate_no_error();
            BOOST_TEST_REQUIRE(reader.has_message());
            auto msg1 = reader.get_next_message(seqnum1, err);
            BOOST_TEST_REQUIRE(err == error_code());

            // Read and get 2nd message
            fn.read_some(fix.stream, reader).validate_no_error();
            BOOST_TEST_REQUIRE(reader.has_message());
            auto msg2 = reader.get_next_message(seqnum2, err);
            BOOST_TEST_REQUIRE(err == error_code());
            BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);
            BOOST_TEST(seqnum2 == 6);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg2, msg2_body);
            BOOST_TEST(!reader.has_message());

            // The 2nd message is located where the 1st message was
            BOOST_TEST(msg1.data() == msg2.data());
        }
    }
}

BOOST_AUTO_TEST_CASE(error)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            fix.inner_stream().set_fail_count(fail_count(0, client_errc::wrong_num_params));

            // Read with error
            fn.read_some(fix.stream, reader).validate_error_exact(client_errc::wrong_num_params);
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

// Cases specific to get_next_message
BOOST_AUTO_TEST_SUITE(get_next_message)

BOOST_AUTO_TEST_CASE(multiframe_message)
{
    fixture fix;
    message_reader reader(512, 8);
    std::uint8_t seqnum = 2;
    fix.inner_stream()
        .add_bytes(create_frame(seqnum, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08}))
        .add_bytes(create_frame(3, {0x09, 0x0a}));
    std::vector<std::uint8_t> expected_msg{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a};
    error_code err(client_errc::server_unsupported);

    // Read succesfully
    reader.read_some(fix.stream, err);
    BOOST_TEST(err == error_code());
    BOOST_TEST_REQUIRE(reader.has_message());
    BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);

    // Get next message and validate it
    auto msg = reader.get_next_message(seqnum, err);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(seqnum == 4u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, expected_msg);

    // There isn't another message
    BOOST_TEST(!reader.has_message());
}

BOOST_AUTO_TEST_CASE(seqnum_overflow)
{
    fixture fix;
    message_reader reader(512);
    std::uint8_t seqnum = 0xff;
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
    fix.inner_stream().add_bytes(create_frame(seqnum, msg_body));
    error_code err(client_errc::server_unsupported);

    // Read succesfully
    reader.read_some(fix.stream, err);
    BOOST_TEST(err == error_code());
    BOOST_TEST_REQUIRE(reader.has_message());
    BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);

    // Get next message and validate it
    auto msg = reader.get_next_message(seqnum, err);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(seqnum == 0u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, msg_body);
}

BOOST_AUTO_TEST_CASE(error_passed_seqnum_mismatch)
{
    fixture fix;
    message_reader reader(512);
    fix.inner_stream().add_bytes(create_frame(2, {0x01, 0x02, 0x03}));
    error_code err(client_errc::server_unsupported);

    // Read succesfully
    reader.read_some(fix.stream, err);
    BOOST_TEST(err == error_code());
    BOOST_TEST_REQUIRE(reader.has_message());
    BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);

    // Passed-in seqnum is invalid
    std::uint8_t bad_seqnum = 0;
    reader.get_next_message(bad_seqnum, err);
    BOOST_TEST(err == make_error_code(client_errc::sequence_number_mismatch));
    BOOST_TEST(bad_seqnum == 0u);
}

BOOST_AUTO_TEST_CASE(error_intermediate_frame_seqnum_mismatch)
{
    fixture fix;
    message_reader reader(512, 8);  // frames are broken each 8 bytes
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
    std::uint8_t seqnum = 2;
    fix.inner_stream()
        .add_bytes(create_frame(seqnum, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08}))
        .add_bytes(create_frame(4, {0x11, 0x12, 0x13, 0x14}));  // the right seqnum would be 3
    error_code err(client_errc::server_unsupported);

    // Read succesfully
    reader.read_some(fix.stream, err);
    BOOST_TEST(err == error_code());
    BOOST_TEST_REQUIRE(reader.has_message());
    BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);

    // The read frame has a mismatched seqnum
    reader.get_next_message(seqnum, err);
    BOOST_TEST(err == client_errc::sequence_number_mismatch);
    BOOST_TEST(seqnum == 2u);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(read_one)

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
            fix.inner_stream().add_bytes(create_frame(seqnum, msg_body));

            // Read succesfully
            auto msg = fn.read_one(fix.stream, reader, seqnum).get();
            BOOST_TEST(seqnum == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, msg_body);

            BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(cached_message)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            std::uint8_t seqnum1 = 2;
            std::uint8_t seqnum2 = 8;
            std::vector<std::uint8_t> msg1_body{0x01, 0x02, 0x03};
            std::vector<std::uint8_t> msg2_body{0x04, 0x05};
            fix.inner_stream()
                .add_bytes(create_frame(seqnum1, msg1_body))
                .add_bytes(create_frame(seqnum2, msg2_body));

            // Read succesfully
            auto msg = fn.read_one(fix.stream, reader, seqnum1).get();
            BOOST_TEST(seqnum1 == 3u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, msg1_body);
            BOOST_TEST(reader.has_message());

            // Read again
            msg = fn.read_one(fix.stream, reader, seqnum2).get();
            BOOST_TEST(seqnum2 == 9u);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, msg2_body);
            BOOST_TEST(fix.inner_stream().num_unread_bytes() == 0u);
            BOOST_TEST(!reader.has_message());
        }
    }
}

BOOST_AUTO_TEST_CASE(error_in_read)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            std::uint8_t seqnum = 2;
            std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
            fix.inner_stream()
                .add_bytes(create_frame(seqnum, msg_body))
                .set_fail_count(fail_count(0, client_errc::wrong_num_params));

            fn.read_one(fix.stream, reader, seqnum).validate_error_exact(client_errc::wrong_num_params);
            BOOST_TEST(seqnum == 2u);
        }
    }
}

BOOST_AUTO_TEST_CASE(seqnum_mismatch)
{
    for (auto fn : all_fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            fixture fix;
            message_reader reader(512);
            fix.inner_stream().add_bytes(create_frame(2, {0x01, 0x02, 0x03}));

            std::uint8_t bad_seqnum = 42;
            fn.read_one(fix.stream, reader, bad_seqnum)
                .validate_error_exact(client_errc::sequence_number_mismatch);
            BOOST_TEST(bad_seqnum == 42u);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
