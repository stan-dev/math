//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/impl/internal/sansio/message_reader.hpp>
#include <boost/mysql/impl/internal/sansio/read_buffer.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_frame.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::span;
using boost::mysql::client_errc;
using boost::mysql::error_code;

using u8vec = std::vector<std::uint8_t>;

BOOST_AUTO_TEST_SUITE(test_message_reader)

class reader_fixture
{
public:
    message_reader reader;
    std::uint8_t seqnum{42};

    reader_fixture(
        std::vector<std::uint8_t> contents,
        std::size_t buffsize = 512,
        std::size_t max_size = static_cast<std::size_t>(-1)
    )
        : reader(buffsize, max_size, 64),  // max frame size is 64
          contents_(std::move(contents)),
          buffer_first_(reader.internal_buffer().first())
    {
    }

    void set_contents(u8vec value)
    {
        contents_ = std::move(value);
        bytes_written_ = 0u;
    }

    // Reads bytes until reader.done() or all bytes in contents have been read.
    // Resizes the buffer as required
    void read_until_completion()
    {
        while (!reader.done() && remaining_bytes())
        {
            auto ec = reader.prepare_buffer();
            BOOST_TEST(ec == error_code());
            std::size_t bytes_to_copy = (std::min)(reader.buffer().size(), contents_.size() - bytes_written_);
            read_bytes(bytes_to_copy);
        }
        BOOST_TEST(reader.done());
        BOOST_TEST(remaining_bytes() == 0u);
    }

    // Simulates a read of num_bytes against the read buffer, then processes the result.
    // Doesn't resize the buffer
    void read_bytes(std::size_t num_bytes)
    {
        // Simulate a write against the buffer
        if (num_bytes)
        {
            assert(num_bytes <= reader.buffer().size());
            std::memcpy(reader.buffer().data(), contents_.data() + bytes_written_, num_bytes);
            bytes_written_ += num_bytes;
        }

        // Trigger the op
        reader.resume(num_bytes);
    }

    span<const std::uint8_t> check_message(const std::vector<std::uint8_t>& expected)
    {
        BOOST_TEST_REQUIRE(reader.done());
        BOOST_TEST_REQUIRE(reader.error() == error_code());
        auto msg = reader.message();
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, expected);
        return msg;
    }

    void record_buffer_first() noexcept { buffer_first_ = reader.internal_buffer().first(); }
    void check_buffer_stability() { BOOST_TEST(reader.internal_buffer().first() == buffer_first_); }
    std::size_t buffsize() const noexcept { return reader.internal_buffer().size(); }

private:
    std::vector<std::uint8_t> contents_;
    std::size_t bytes_written_{0};
    const std::uint8_t* buffer_first_;

    std::size_t remaining_bytes() const noexcept { return contents_.size() - bytes_written_; }
};

// Parsing algorithm. Without buffer relocations or short reads
BOOST_AUTO_TEST_CASE(parsing_algorithm_success)
{
    struct
    {
        const char* name;
        u8vec input;
        u8vec expected_msg;
        std::uint8_t expected_seqnum;
    } test_cases[] = {
        {"empty_message", create_empty_frame(42), {}, 43},
        {"one_frame", create_frame(42, {0x01, 0x02, 0x03}), {0x01, 0x02, 0x03}, 43},
        {"one_frame_max_size",
         buffer_builder().add(create_frame(42, u8vec(64, 0x04))).add(create_empty_frame(43)).build(),
         u8vec(64, 0x04),
         44},
        {"two_frames",
         buffer_builder().add(create_frame(42, u8vec(64, 0x04))).add(create_frame(43, {0x05, 0x06})).build(),
         concat(u8vec(64, 0x04), {0x05, 0x06}),
         44},
        {"two_frames_max_size",
         buffer_builder()
             .add(create_frame(42, u8vec(64, 0x04)))
             .add(create_frame(43, u8vec(64, 0x05)))
             .add(create_empty_frame(44))
             .build(),
         concat(u8vec(64, 0x04), u8vec(64, 0x05)),
         45},
        {"three_frames",
         buffer_builder()
             .add(create_frame(42, u8vec(64, 0x04)))
             .add(create_frame(43, u8vec(64, 0x05)))
             .add(create_frame(44, {0x0a}))
             .build(),
         buffer_builder().add(u8vec(64, 0x04)).add(u8vec(64, 0x05)).add({0x0a}).build(),
         45},
    };

    for (auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            reader_fixture fix(tc.input);
            fix.reader.prepare_read(fix.seqnum);
            BOOST_TEST(!fix.reader.done());

            // Receive the message
            fix.read_bytes(tc.input.size());

            // Check
            fix.check_message(tc.expected_msg);
            BOOST_TEST(fix.seqnum == tc.expected_seqnum);

            // Buffer didn't reallocate
            fix.check_buffer_stability();
        }
    }
}

BOOST_AUTO_TEST_CASE(seqnum_overflow)
{
    // message to be parsed
    reader_fixture fix(
        buffer_builder()
            .add(create_frame(255, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(0, {0x05, 0x06, 0x07}))
            .build(),
        64 + 16
    );
    fix.seqnum = 255;

    // Setup
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());

    // all in one
    fix.read_bytes(64 + 4 * 2 + 3);
    fix.check_message(concat(std::vector<std::uint8_t>(64, 0x04), {0x05, 0x06, 0x07}));
    BOOST_TEST(fix.seqnum == 1u);

    // Buffer didn't reallocate
    fix.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(seqnum_mismatch)
{
    struct
    {
        const char* name;
        u8vec input;
    } test_cases[] = {
        {"1st_frame", create_frame(1, {0x01, 0x02})},
        {"2nd_frame", concat(create_frame(42, u8vec(64, 0x04)), create_frame(44, {0x01}))},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            reader_fixture fix(tc.input);
            fix.reader.prepare_read(fix.seqnum);
            fix.read_bytes(tc.input.size());
            BOOST_TEST(fix.reader.done());
            BOOST_TEST(fix.reader.error() == error_code(client_errc::sequence_number_mismatch));
        }
    }
}

// Long read: we received two messages at once.
// We don't consume the next message while parsing the first one
// We don't get rid of the first message while there's space for the second one
BOOST_AUTO_TEST_CASE(long_read)
{
    // message to be parsed
    std::vector<std::uint8_t> first_msg_body{0x01, 0x02, 0x03};
    std::vector<std::uint8_t> second_msg_body{0x04, 0x05, 0x06, 0x07};
    reader_fixture fix(concat(create_frame(42, first_msg_body), create_frame(2, second_msg_body)));

    // Prepare first read
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());

    // The read yields two messages at once
    fix.read_bytes(15);
    fix.check_message(first_msg_body);
    BOOST_TEST(fix.seqnum == 43u);

    // We can read the 2nd message, too
    std::uint8_t seqnum2 = 2u;
    fix.reader.prepare_read(seqnum2);
    fix.check_message(second_msg_body);
    BOOST_TEST(fix.seqnum == 43u);  // old seqnum not updated
    BOOST_TEST(seqnum2 == 3u);      // new seqnum updated

    // Buffer shouldn't reallocate
    fix.check_buffer_stability();
}

// Short reads
BOOST_AUTO_TEST_CASE(short_reads_multiple)
{
    // message to be parsed
    reader_fixture fix(create_frame(42, {0x01, 0x02, 0x03}));
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());

    // 1 byte in the header received
    fix.read_bytes(1);
    BOOST_TEST(!fix.reader.done());

    // Another 2 bytes received
    fix.read_bytes(2);
    BOOST_TEST(!fix.reader.done());

    // Header fully received
    fix.read_bytes(1);
    BOOST_TEST(!fix.reader.done());

    // 1 byte in body received
    fix.read_bytes(1);
    BOOST_TEST(!fix.reader.done());

    // body fully received
    fix.read_bytes(2);
    fix.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(fix.seqnum == 43u);

    // Buffer shouldn't reallocate
    fix.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(short_reads_header_size)
{
    // message to be parsed
    reader_fixture fix(create_frame(42, {0x01, 0x02, 0x03}));
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());

    // Full header received
    fix.read_bytes(4);
    BOOST_TEST(!fix.reader.done());

    // Full body received
    fix.read_bytes(3);
    fix.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(fix.seqnum == 43u);

    // Buffer didn't reallocate
    fix.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(short_reads_two_frames)
{
    // message to be parsed
    reader_fixture fix(
        buffer_builder()
            .add(create_frame(42, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(43, {0x05, 0x06, 0x07}))
            .build(),
        64 + 16
    );
    auto expected_message = concat(std::vector<std::uint8_t>(64, 0x04), {0x05, 0x06, 0x07});

    // Setup
    fix.reader.prepare_read(fix.seqnum);

    // part of header 1
    fix.read_bytes(3);
    BOOST_TEST(!fix.reader.done());

    // header 1 full
    fix.read_bytes(1);
    BOOST_TEST(!fix.reader.done());

    // part of body 1
    fix.read_bytes(64 - 8);
    BOOST_TEST(!fix.reader.done());

    // rest of body 1
    fix.read_bytes(8);
    BOOST_TEST(!fix.reader.done());

    // part of header 2
    fix.read_bytes(1);
    BOOST_TEST(!fix.reader.done());

    // another part of header 2
    fix.read_bytes(2);
    BOOST_TEST(!fix.reader.done());

    // rest of header 2 and part of body 1
    fix.read_bytes(2);
    BOOST_TEST(!fix.reader.done());

    // another part of body 2
    fix.read_bytes(1);
    BOOST_TEST(!fix.reader.done());

    // remaining of body 2
    fix.read_bytes(1);
    fix.check_message(expected_message);
    BOOST_TEST(fix.seqnum == 44u);

    // Buffer shouldn't reallocate
    fix.check_buffer_stability();
}

// Buffer resizing
BOOST_AUTO_TEST_CASE(buffer_resizing_not_enough_space)
{
    // Setup
    reader_fixture fix(create_frame(42, u8vec(50, 0x04)), 0);
    BOOST_TEST(fix.buffsize() == 0u);

    // Prepare read. The buffer hasn't resized.
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());
    BOOST_TEST(fix.buffsize() == 0u);

    // Resize the buffer
    auto ec = fix.reader.prepare_buffer();
    BOOST_TEST(ec == error_code());
    fix.record_buffer_first();
    BOOST_TEST(fix.buffsize() == 4u);

    // Read the header. The buffer didn't reallocate
    fix.read_bytes(4);
    BOOST_TEST(!fix.reader.done());
    fix.check_buffer_stability();

    // Resize the buffer again
    ec = fix.reader.prepare_buffer();
    BOOST_TEST(ec == error_code());
    fix.record_buffer_first();
    BOOST_TEST(fix.buffsize() == 50u);

    // Finish reading
    fix.read_bytes(50);
    fix.check_message(u8vec(50, 0x04));
    BOOST_TEST(fix.seqnum == 43u);
}

BOOST_AUTO_TEST_CASE(buffer_resizing_old_messages_removed)
{
    // prepare_buffer removes old messages
    // so the buffer doesn't grow indefinitely

    // Setup
    reader_fixture fix(create_frame(42, u8vec(60, 0x04)), 0);

    // Parse an entire message, to make space in the buffer
    fix.reader.prepare_read(fix.seqnum);
    fix.read_until_completion();
    fix.check_message(u8vec(60, 0x04));

    // Record size, as this should not increase
    BOOST_TEST(fix.buffsize() == 60u);

    // Parse new messages
    for (std::uint8_t i = 0u; i < 100u; ++i)
    {
        // Setup
        u8vec msg_body(50, i);
        fix.seqnum = i;
        fix.set_contents(create_frame(i, msg_body));

        // Prepare read
        fix.reader.prepare_read(fix.seqnum);

        // Read the message into the buffer and trigger the op until completion.
        // This will call prepare_buffer() internally
        fix.read_until_completion();

        // Check results
        BOOST_TEST_REQUIRE(fix.reader.error() == error_code());
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.reader.message(), msg_body);
    }

    // Buffer size should be the same
    BOOST_TEST(fix.buffsize() == 60u);
}

BOOST_AUTO_TEST_CASE(buffer_resizing_size_eq_max_size)
{
    // Reading a frame of exactly max_size works
    // Setup
    reader_fixture fix(create_frame(42, u8vec(32, 0x04)), 0, 32);
    BOOST_TEST(fix.buffsize() == 0u);

    // Prepare read. The buffer hasn't resized.
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());
    BOOST_TEST(fix.buffsize() == 0u);

    // Execute the read successfully
    fix.read_until_completion();
    fix.check_message(u8vec(32, 0x04));
    BOOST_TEST(fix.seqnum == 43u);
}

BOOST_AUTO_TEST_CASE(buffer_resizing_max_size_exceeded)
{
    // Setup
    reader_fixture fix(create_frame(42, u8vec(50, 0x04)), 16, 32);
    BOOST_TEST(fix.buffsize() == 16u);
    fix.record_buffer_first();

    // Prepare read. The buffer hasn't resized.
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());
    BOOST_TEST(fix.buffsize() == 16u);
    fix.check_buffer_stability();

    // We have enough size for the header
    auto ec = fix.reader.prepare_buffer();
    BOOST_TEST(ec == error_code());
    BOOST_TEST(fix.buffsize() == 16u);
    fix.record_buffer_first();

    // Read the header. The buffer didn't reallocate
    fix.read_bytes(4);
    BOOST_TEST(!fix.reader.done());
    fix.check_buffer_stability();

    // Resizing the buffer here would require exceeding max size and fails
    ec = fix.reader.prepare_buffer();
    BOOST_TEST(ec == client_errc::max_buffer_size_exceeded);
}

BOOST_AUTO_TEST_CASE(buffer_resizing_max_size_exceeded_subsequent_frames)
{
    // Setup
    reader_fixture fix(create_frame(42, u8vec(90, 0x04)), 80, 80);
    BOOST_TEST(fix.buffsize() == 80u);
    fix.record_buffer_first();

    // Prepare read
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());

    // Read header, body 1, header 2, part of body 2
    auto ec = fix.reader.prepare_buffer();
    BOOST_TEST(ec == error_code());
    fix.read_bytes(80);
    BOOST_TEST(!fix.reader.done());

    // Resizing the buffer here would require exceeding max size and fails
    ec = fix.reader.prepare_buffer();
    BOOST_TEST(ec == client_errc::max_buffer_size_exceeded);
}

// Keep parsing state
BOOST_AUTO_TEST_CASE(keep_state_continuation)
{
    // Setup. We use a multiframe message to verify that we update the sequence number reference correctly
    u8vec msg1_body{0x01, 0x02, 0x03};
    u8vec msg2_body(65, 0x04);
    reader_fixture fix(
        buffer_builder()
            .add(create_frame(42, msg1_body))
            .add(create_frame(43, u8vec(64, 0x04)))
            .add(create_frame(44, {0x04}))
            .build(),
        16
    );

    // Read the first message and part of the second
    fix.reader.prepare_read(fix.seqnum);
    fix.read_bytes(16);
    auto msg = fix.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(fix.seqnum == 43u);

    // Prepare the second read. We don't have enough bytes or buffer space
    fix.reader.prepare_read(fix.seqnum);
    BOOST_TEST(!fix.reader.done());
    fix.check_buffer_stability();                      // Didn't reallocate
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, msg1_body);  // Old message still valid
    BOOST_TEST(fix.seqnum == 44u);                     // Updated to the last received seqnum

    // Prepare a read as a continuation. This will not throw away the parsing state
    std::uint8_t new_seqnum{fix.seqnum};
    fix.reader.prepare_read(new_seqnum, true);
    fix.read_until_completion();
    fix.check_message(msg2_body);
    BOOST_TEST(fix.seqnum == 44u);  // Old seqnum not updated
    BOOST_TEST(new_seqnum == 45u);  // New seqnum updated
}

BOOST_AUTO_TEST_CASE(keep_state_done)
{
    // Passing keep_state=true won't have effect if the operation is already done
    reader_fixture fix(
        buffer_builder().add(create_frame(42, {0x01, 0x02, 0x03})).add(create_frame(43, {0x04, 0x05})).build()
    );

    // Read the first message
    fix.reader.prepare_read(fix.seqnum);
    fix.read_bytes(7);
    fix.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(fix.seqnum == 43u);

    // Prepare a read as a continuation. As the operation is done, this will reset parsing state
    fix.reader.prepare_read(fix.seqnum, true);
    BOOST_TEST(!fix.reader.done());
    fix.read_bytes(6);
    fix.check_message({0x04, 0x05});
    BOOST_TEST(fix.seqnum == 44u);
}

BOOST_AUTO_TEST_CASE(keep_state_initial)
{
    // Passing keep_state=true with a reader that hasn't been used works
    reader_fixture fix(create_frame(42, {0x01, 0x02, 0x03}));
    fix.reader.prepare_read(fix.seqnum, true);
    fix.read_bytes(7);
    fix.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(fix.seqnum == 43u);
}

// Resetting
BOOST_AUTO_TEST_CASE(reset_done)
{
    // Read a message until completion
    reader_fixture fix(create_frame(42, {0x01, 0x02, 0x03}));
    fix.reader.prepare_read(fix.seqnum);
    fix.read_until_completion();

    // Reset
    fix.reader.reset();

    // A new message can be read now
    fix.set_contents(create_frame(20, {0x09, 0x0a}));
    fix.seqnum = 20;
    fix.reader.prepare_read(fix.seqnum);
    fix.read_until_completion();
    fix.check_message({0x09, 0x0a});
    fix.check_buffer_stability();  // No reallocation happened
    BOOST_TEST(fix.seqnum == 21u);
}

BOOST_AUTO_TEST_CASE(reset_message_half_read)
{
    // Read part of a message
    reader_fixture fix(create_frame(42, {0x01, 0x02, 0x03}));
    fix.reader.prepare_read(fix.seqnum);
    fix.read_bytes(3);

    // Reset
    fix.reader.reset();

    // A new message can be read now
    fix.set_contents(create_frame(20, {0x09, 0x0a}));
    fix.seqnum = 20;
    fix.reader.prepare_read(fix.seqnum);
    fix.read_until_completion();
    fix.check_message({0x09, 0x0a});
    fix.check_buffer_stability();  // No reallocation happened
    BOOST_TEST(fix.seqnum == 21u);
}

BOOST_AUTO_TEST_CASE(reset_keep_state_true)
{
    // Read part of a message
    reader_fixture fix(create_frame(42, {0x01, 0x02, 0x03}));
    fix.reader.prepare_read(fix.seqnum, true);
    fix.read_bytes(3);

    // Reset
    fix.reader.reset();

    // A new message can be read now
    fix.set_contents(create_frame(20, {0x09, 0x0a}));
    fix.seqnum = 20;
    fix.reader.prepare_read(fix.seqnum);
    fix.read_until_completion();
    fix.check_message({0x09, 0x0a});
    fix.check_buffer_stability();  // No reallocation happened
    BOOST_TEST(fix.seqnum == 21u);
}

BOOST_AUTO_TEST_SUITE_END()
