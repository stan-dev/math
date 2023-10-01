//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/impl/internal/channel/read_buffer.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"

using namespace boost::mysql::detail;

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
}

static void check_buffer(
    read_buffer& buff,
    const std::vector<std::uint8_t>& reserved,
    const std::vector<std::uint8_t>& current_message,
    const std::vector<std::uint8_t>& pending
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
        buff.size() - free_offset
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

    BOOST_TEST(buff.free_size() == buff.size());
    BOOST_TEST(buff.size() >= 531u);
    check_buffer(buff, {}, {}, {});
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
    buff.grow_to_fit(0);
    check_empty_buffer(buff);
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

    check_buffer(buff, {}, {}, contents);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(all_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    std::vector<std::uint8_t> contents(buff.size(), 0x01);
    copy_to_free_area(buff, contents);
    buff.move_to_pending(buff.size());

    check_buffer(buff, {}, {}, contents);
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(zero_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    buff.move_to_pending(0);

    check_buffer(buff, {}, {}, {});
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

    check_buffer(buff, {}, {}, contents);
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

    check_buffer(buff, {}, {0x01, 0x02}, {0x03, 0x04, 0x05, 0x06});
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(all_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(6);

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {});
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(zero_bytes)
{
    read_buffer buff(8);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
    buff.move_to_pending(6);
    buff.move_to_current_message(0);

    check_buffer(buff, {}, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06});
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

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05}, {0x06});
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

    check_buffer(buff, {0x01, 0x02, 0x03}, {0x04, 0x05}, {0x06});
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

    check_buffer(buff, {0x01, 0x02, 0x03, 0x04, 0x05}, {}, {0x06});
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

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05}, {0x06});
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

    check_buffer(buff, {0x01, 0x02, 0x03}, {0x04, 0x05}, {0x06});
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

    check_buffer(buff, {0x01}, {0x02, 0x03, 0x04}, {0x07, 0x08});
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

    check_buffer(buff, {0x01}, {}, {0x07, 0x08});
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

    check_buffer(buff, {0x01}, {0x02, 0x03, 0x04}, {});
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

    check_buffer(buff, {}, {0x01, 0x02}, {0x07, 0x08});
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

    check_buffer(buff, {}, {0x03, 0x04, 0x05, 0x06}, {0x07, 0x08});
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

    check_buffer(buff, {}, {}, {});
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

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08});
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
    buff.grow_to_fit(100);

    BOOST_TEST(buff.free_size() >= 100u);
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08});
    checker.check_reallocation();
}

BOOST_AUTO_TEST_CASE(one_missing_byte)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);

    std::size_t required_size = buff.size() - 8 + 1;
    buff.grow_to_fit(required_size);

    BOOST_TEST(buff.free_size() >= required_size);
    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08});
    checker.check_reallocation();
}

BOOST_AUTO_TEST_CASE(enough_space)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    buff.grow_to_fit(8);

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08});
    checker.check_stability();
}

BOOST_AUTO_TEST_CASE(zero_bytes)
{
    read_buffer buff(16);
    stability_checker checker(buff);
    copy_to_free_area(buff, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08});
    buff.move_to_pending(8);
    buff.move_to_current_message(6);
    buff.grow_to_fit(0);

    check_buffer(buff, {}, {0x01, 0x02, 0x03, 0x04, 0x05, 0x06}, {0x07, 0x08});
    checker.check_stability();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
