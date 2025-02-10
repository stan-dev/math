//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/impl/internal/sansio/read_buffer.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/printing.hpp"

using namespace boost::mysql::detail;
using boost::mysql::client_errc;
using boost::mysql::error_code;

BOOST_AUTO_TEST_SUITE(test_read_buffer)

// Records the buffer first pointer and size to verify the buffer
// didn't do any re-allocation
class stability_checker
{
    read_buffer& b_;
    const std::uint8_t* first_;
    std::size_t total_size_;

public:
    stability_checker(read_buffer& b) noexcept : b_(b), first_(b.first()), total_size_(b.size()) {}
    void check_stability()
    {
        BOOST_TEST(b_.first() == first_);
        BOOST_TEST(b_.size() == total_size_);
    }
    void check_reallocation()
    {
        BOOST_TEST(b_.first() != first_);
        BOOST_TEST(b_.size() > total_size_);
    }
};

static void check_buffer(
    read_buffer& buff,
    const std::uint8_t* reserved_first,
    const std::uint8_t* current_message_first,
    const std::uint8_t* pending_first,
    const std::uint8_t* free_first,
    std::size_t reserved_size,
    std::size_t current_message_size,
    std::size_t pending_size,
    std::size_t free_size
)
{
    BOOST_TEST(buff.reserved_first() == reserved_first);
    BOOST_TEST(buff.current_message_first() == current_message_first);
    BOOST_TEST(buff.pending_first() == pending_first);
    BOOST_TEST(buff.free_first() == free_first);

    BOOST_TEST(buff.reserved_area().data() == reserved_first);
    BOOST_TEST(buff.current_message().data() == current_message_first);
    BOOST_TEST(buff.pending_area().data() == pending_first);
    BOOST_TEST(buff.free_area().data() == free_first);

    BOOST_TEST(buff.reserved_size() == reserved_size);
    BOOST_TEST(buff.current_message_size() == current_message_size);
    BOOST_TEST(buff.pending_size() == pending_size);
    BOOST_TEST(buff.free_size() == free_size);

    BOOST_TEST(buff.reserved_area().size() == reserved_size);
    BOOST_TEST(buff.current_message().size() == current_message_size);
    BOOST_TEST(buff.pending_area().size() == pending_size);
    BOOST_TEST(buff.free_area().size() == free_size);

    BOOST_TEST(buff.size() == reserved_size + current_message_size + pending_size + free_size);
}

static void check_buffer(
    read_buffer& buff,
    const std::vector<std::uint8_t>& reserved,
    const std::vector<std::uint8_t>& current_message,
    const std::vector<std::uint8_t>& pending,
    std::size_t free_size
)
{
    std::size_t current_message_offset = reserved.size();
    std::size_t pending_offset = current_message_offset + current_message.size();
    std::size_t free_offset = pending_offset + pending.size();

    BOOST_TEST(buff.first() != nullptr);

    check_buffer(
        buff,
        buff.first(),
        buff.first() + current_message_offset,
        buff.first() + pending_offset,
        buff.first() + free_offset,
        reserved.size(),
        current_message.size(),
        pending.size(),
        free_size
    );

    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff.reserved_area(), reserved);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff.current_message(), current_message);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff.pending_area(), pending);
}

static void check_empty_buffer(read_buffer& buff)
{
    check_buffer(buff, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, 0);
}

static void copy_to_free_area(read_buffer& buff, const std::vector<std::uint8_t>& bytes)
{
    std::copy(bytes.begin(), bytes.end(), buff.free_first());
}

BOOST_AUTO_TEST_SUITE(init_ctor)

BOOST_AUTO_TEST_CASE(some_initial_size)
{
    read_buffer buff(531);
    check_buffer(buff, {}, {}, {}, 531);
}

BOOST_AUTO_TEST_CASE(zero_initial_size)
{
    read_buffer buff(0);

    check_empty_buffer(buff);

    // Calling all other functions with 0 values on this buffer doesn't cause UB
    buff.move_to_pending(0);
    buff.move_to_current_message(0);
    buff.move_to_reserved(0);
    buff.remove_reserved();
    auto ec = buff.grow_to_fit(0);
    BOOST_TEST(ec == error_code());
    check_empty_buffer(buff);
}

BOOST_AUTO_TEST_CASE(initial_size_eq_max_size)
{
    read_buffer buff(16, 16);
    check_buffer(buff, {}, {}, {}, 16);
    BOOST_TEST(buff.max_size() == 16u);

    // Using the buffer works normally
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04});
    buff.move_to_pending(4);
    buff.move_to_current_message(3);
    buff.move_to_reserved(1);
    check_buffer(buff, {0x01}, {0x02, 0x03}, {0x04}, 12);

    // Growing works
    auto ec = buff.grow_to_fit(12);
    BOOST_TEST(ec == error_code());
    ec = buff.grow_to_fit(13);
    BOOST_TEST(ec == client_errc::max_buffer_size_exceeded);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_to_pending)

BOOST_AUTO_TEST_CASE(some_bytes)
{
    read_buffer buff(512);
    stability_checker checker(buff);
    std::vector<std::uint8_t> contents{0x01, 0x02, 0x03, 0x04};
    copy_to_free_area(buff, contents);
    buff.move_to_pending(4);

    check_buffer(buff, {}, {}, contents, 508);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(all_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    std::vector<std::uint8_t> contents(8, 0x01);
    copy_to_free_area(buff, contents);
    buff.move_to_pending(buff.size());

    check_buffer(buff, {}, {}, contents, 0);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(zero_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    buff.move_to_pending(0);

    check_buffer(buff, {}, {}, {}, 8);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(several_calls)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    std::vector<std::uint8_t> contents{0x01, 0x02, 0x03, 0x04};
    copy_to_free_area(buff, contents);
    buff.move_to_pending(2);
    buff.move_to_pending(2);

    check_buffer(buff, {}, {}, contents, 4);
    checker.check_stability();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_to_current_message)

BOOST_AUTO_TEST_CASE(some_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(2);

    check_buffer(buff, {}, {0x01, 0x02}, {0x03, 0x04, 0x05, 0x06}, 2);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(all_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(6);

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {}, 2);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(zero_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(0);

    check_buffer(buff, {}, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, 2);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(several_calls)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(2);
    buff.move_to_current_message(3);

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05}, {0x06}, 2);
    checker.check_stability();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_to_reserved)

BOOST_AUTO_TEST_CASE(some_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(5);
    buff.move_to_reserved(3);

    check_buffer(buff, {0x01, 0x02, 0x03}, {0x04, 0x05}, {0x06}, 2);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(all_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(5);
    buff.move_to_reserved(5);

    check_buffer(buff, {0x01, 0x02, 0x03, 0x04, 0x05}, {}, {0x06}, 2);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(zero_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(5);
    buff.move_to_reserved(0);

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05}, {0x06}, 2);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(several_calls)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(5);
    buff.move_to_reserved(1);
    buff.move_to_reserved(2);

    check_buffer(buff, {0x01, 0x02, 0x03}, {0x04, 0x05}, {0x06}, 2);
    checker.check_stability();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(remove_current_message_last)

BOOST_AUTO_TEST_CASE(some_bytes)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    buff.move_to_reserved(1);
    buff.remove_current_message_last(2);

    check_buffer(buff, {0x01}, {0x02, 0x03, 0x04}, {0x07, 0x08}, 10);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(all_bytes)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    buff.move_to_reserved(1);
    buff.remove_current_message_last(5);

    check_buffer(buff, {0x01}, {}, {0x07, 0x08}, 13);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(without_pending)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(8);
    buff.move_to_reserved(1);
    buff.remove_current_message_last(4);

    check_buffer(buff, {0x01}, {0x02, 0x03, 0x04}, {}, 12);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(without_reserved)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    buff.remove_current_message_last(4);

    check_buffer(buff, {}, {0x01, 0x02}, {0x07, 0x08}, 12);
    checker.check_stability();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(remove_reserved)

BOOST_AUTO_TEST_CASE(with_other_areas)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    buff.move_to_reserved(2);
    buff.remove_reserved();

    check_buffer(buff, {}, {0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 10);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(without_other_areas)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(8);
    buff.move_to_reserved(8);
    buff.remove_reserved();

    check_buffer(buff, {}, {}, {}, 16);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(zero_bytes)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    buff.remove_reserved();

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 8);
    checker.check_stability();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(grow_to_fit)

BOOST_AUTO_TEST_CASE(not_enough_space)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    auto ec = buff.grow_to_fit(100);

    BOOST_TEST(ec == error_code());
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 100);
    checker.check_reallocation();
}

BOOST_AUTO_TEST_CASE(one_missing_byte)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);

    auto ec = buff.grow_to_fit(9);

    BOOST_TEST(ec == error_code());
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 9);
    checker.check_reallocation();
}

BOOST_AUTO_TEST_CASE(enough_space)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    auto ec = buff.grow_to_fit(8);
    BOOST_TEST(ec == error_code());

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 8);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(zero_bytes)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    auto ec = buff.grow_to_fit(0);
    BOOST_TEST(ec == error_code());

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 8);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(from_size_0)
{
    // Regression check: growing from size 0 works
    read_buffer buff(0, 1024);
    check_empty_buffer(buff);

    auto ec = buff.grow_to_fit(16);
    BOOST_TEST(ec == error_code());
    check_buffer(buff, {}, {}, {}, 16);
}

BOOST_AUTO_TEST_CASE(lt_max_size)
{
    read_buffer buff(8, 16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);

    // Grow past the current size, but not reaching max size
    auto ec = buff.grow_to_fit(7);
    BOOST_TEST(ec == error_code());

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 7);
    checker.check_reallocation();
}

BOOST_AUTO_TEST_CASE(eq_max_size)
{
    read_buffer buff(8, 16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);

    // Grow past the current size, reaching max_size
    auto ec = buff.grow_to_fit(8);
    BOOST_TEST(ec == error_code());

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 8);
    checker.check_reallocation();
}

BOOST_AUTO_TEST_CASE(gt_max_size)
{
    read_buffer buff(8, 16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);

    // Try to grow past max size
    auto ec = buff.grow_to_fit(10);
    BOOST_TEST(ec == client_errc::max_buffer_size_exceeded);
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 0);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(several_grows)
{
    read_buffer buff(8, 16);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);

    // Grow with reallocation
    auto ec = buff.grow_to_fit(4);
    BOOST_TEST(ec == error_code());
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08}, 4);

    // Place some more bytes in the buffer
    copy_to_free_area(buff, {0x09, 0x0a});
    buff.move_to_pending(2);
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08, 0x09, 0x0a}, 2);

    // Grow without reallocation
    ec = buff.grow_to_fit(2);
    BOOST_TEST(ec == error_code());
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08, 0x09, 0x0a}, 2);
    copy_to_free_area(buff, {0x0b, 0x0c});
    buff.move_to_pending(2);
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c}, 0);

    // Fail when attempting to grow past max size
    ec = buff.grow_to_fit(5);
    BOOST_TEST(ec == client_errc::max_buffer_size_exceeded);
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c}, 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(reset)

BOOST_AUTO_TEST_CASE(zero_size_buffer)
{
    read_buffer buff(0, 1024);
    buff.reset();
    BOOST_TEST(buff.size() == 0u);
    BOOST_TEST(buff.max_size() == 1024u);
}

BOOST_AUTO_TEST_CASE(free_buffer)
{
    read_buffer buff(16, 1024);
    stability_checker checker(buff);
    buff.reset();
    check_buffer(buff, {}, {}, {}, 16);
    checker.check_stability();
    BOOST_TEST(buff.max_size() == 1024u);
}

BOOST_AUTO_TEST_CASE(pending_bytes)
{
    read_buffer buff(16, 512);
    stability_checker checker(buff);
    buff.move_to_pending(4);
    buff.reset();
    check_buffer(buff, {}, {}, {}, 16);
    checker.check_stability();
    BOOST_TEST(buff.max_size() == 512u);
}

BOOST_AUTO_TEST_CASE(current_message_bytes)
{
    read_buffer buff(16, 512);
    stability_checker checker(buff);
    buff.move_to_pending(4);
    buff.move_to_current_message(4);
    buff.reset();
    check_buffer(buff, {}, {}, {}, 16);
    checker.check_stability();
    BOOST_TEST(buff.max_size() == 512u);
}

BOOST_AUTO_TEST_CASE(reserved_bytes)
{
    read_buffer buff(16, 512u);
    stability_checker checker(buff);
    buff.move_to_pending(4);
    buff.move_to_current_message(4);
    buff.move_to_reserved(4);
    buff.reset();
    check_buffer(buff, {}, {}, {}, 16);
    checker.check_stability();
    BOOST_TEST(buff.max_size() == 512u);
}

BOOST_AUTO_TEST_CASE(bytes_in_all_areas)
{
    read_buffer buff(16, 1024);
    stability_checker checker(buff);
    buff.move_to_pending(10);
    buff.move_to_current_message(8);
    buff.move_to_reserved(2);
    buff.reset();
    check_buffer(buff, {}, {}, {}, 16);
    checker.check_stability();
    BOOST_TEST(buff.max_size() == 1024u);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
