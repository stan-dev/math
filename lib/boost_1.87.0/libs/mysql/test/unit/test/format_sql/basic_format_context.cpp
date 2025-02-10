//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/format_sql.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/printing.hpp"

using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_basic_format_context)

static format_options opts{utf8mb4_charset, true};
static format_options ascii_opts{ascii_charset, true};

// An OutputString with the bare minimum members
struct string_archetype
{
    std::string impl;

    // No default constructor
    string_archetype(int) {}

    // No copy operations
    string_archetype(const string_archetype&) = delete;
    string_archetype& operator=(const string_archetype&) = delete;

    // Move operations allowed
    string_archetype(string_archetype&&) = default;
    string_archetype& operator=(string_archetype&&) = default;

    // Adequate member functions
    void append(const char* data, std::size_t size) { impl.append(data, size); }
    void clear() { impl.clear(); }
};

// An OutputString that is also default constructible
struct string_archetype_defctor : string_archetype
{
    string_archetype_defctor() : string_archetype(42) {}
};

using archetype_context = basic_format_context<string_archetype>;

BOOST_AUTO_TEST_CASE(ctor_without_storage)
{
    // Construct
    basic_format_context<string_archetype_defctor> ctx(opts);

    // Check error state and options
    BOOST_TEST(ctx.error_state() == error_code());
    BOOST_TEST(ctx.format_opts().charset.name == "utf8mb4");

    // Can be used to append and result can be retrieved
    ctx.append_raw("SELECT ").append_value(42);
    BOOST_TEST(std::move(ctx).get().value().impl == "SELECT 42");
}

BOOST_AUTO_TEST_CASE(ctor_with_storage)
{
    // The string to move from. Heap helps asan detect memory errors
    std::unique_ptr<string_archetype> s(new string_archetype(100));
    s->impl = "abcd";

    // Construct
    archetype_context ctx(opts, std::move(*s));
    s.reset();

    // Check error state and options
    BOOST_TEST(ctx.error_state() == error_code());
    BOOST_TEST(ctx.format_opts().charset.name == "utf8mb4");

    // Can be used to append and result can be retrieved
    ctx.append_raw("SELECT ").append_value(42);
    BOOST_TEST(std::move(ctx).get().value().impl == "SELECT 42");
}

BOOST_AUTO_TEST_CASE(move_constructor)
{
    // Move source. Heap helps asan
    std::unique_ptr<archetype_context> source{new archetype_context(opts, string_archetype(42))};
    source->append_raw("SELECT ");

    // Construct
    archetype_context ctx(std::move(*source));
    source.reset();

    // Check error state and options
    BOOST_TEST(ctx.error_state() == error_code());
    BOOST_TEST(ctx.format_opts().charset.name == "utf8mb4");

    // Can be used to append and result can be retrieved
    ctx.append_value(42);
    BOOST_TEST(std::move(ctx).get().value().impl == "SELECT 42");
}

BOOST_AUTO_TEST_CASE(move_constructor_error)
{
    // Move source. Heap helps asan
    std::unique_ptr<archetype_context> source{new archetype_context(opts, string_archetype(42))};
    source->add_error(client_errc::extra_bytes);

    // Construct
    archetype_context ctx(std::move(*source));
    source.reset();

    // Error state is propagated
    BOOST_TEST(ctx.error_state() == client_errc::extra_bytes);
    BOOST_TEST(std::move(ctx).get().error() == client_errc::extra_bytes);
}

BOOST_AUTO_TEST_CASE(move_assign)
{
    // Move source. Heap helps asan
    std::unique_ptr<archetype_context> source{new archetype_context(opts, string_archetype(42))};
    source->append_raw("SELECT ");
    archetype_context ctx(ascii_opts, string_archetype(42));
    ctx.append_raw("abc").add_error(client_errc::wrong_num_params);

    // Assign
    ctx = std::move(*source);
    source.reset();

    // Check error state and options
    BOOST_TEST(ctx.error_state() == error_code());
    BOOST_TEST(ctx.format_opts().charset.name == "utf8mb4");

    // Can be used to append and result can be retrieved
    ctx.append_value(42);
    BOOST_TEST(std::move(ctx).get().value().impl == "SELECT 42");
}

BOOST_AUTO_TEST_CASE(move_assign_error)
{
    // Move source. Heap helps asan
    std::unique_ptr<archetype_context> source{new archetype_context(opts, string_archetype(42))};
    source->add_error(client_errc::extra_bytes);
    archetype_context ctx(ascii_opts, string_archetype(42));

    // Assign
    ctx = std::move(*source);
    source.reset();

    // Error state is propagated
    BOOST_TEST(ctx.error_state() == client_errc::extra_bytes);
    BOOST_TEST(std::move(ctx).get().error() == client_errc::extra_bytes);
}

// Spotcheck: format_context move operations work
BOOST_AUTO_TEST_CASE(string_format_context)
{
    // Constructor from storage
    std::string storage = "abcde";
    format_context ctx(opts, std::move(storage));

    // Check
    BOOST_TEST(ctx.error_state() == error_code());
    BOOST_TEST(ctx.format_opts().charset.name == "utf8mb4");

    // Move construction
    ctx.append_raw("SELECT ");
    format_context ctx2(std::move(ctx));
    ctx2.append_value(42);
    BOOST_TEST(std::move(ctx2).get().value() == "SELECT 42");

    // Move assignment
    format_context ctx3(ascii_opts);
    ctx3.append_raw("def");
    ctx = std::move(ctx3);
    BOOST_TEST(std::move(ctx).get().value() == "def");
}

// Spotcheck: append_raw works with empty strings
BOOST_AUTO_TEST_CASE(append_raw_empty_string)
{
    // Context
    format_context ctx(opts);

    // With the empty context
    ctx.append_raw("");

    // With contents already present
    ctx.append_raw("SELECT ").append_value(42).append_raw("");

    BOOST_TEST(std::move(ctx).get().value() == "SELECT 42");
}

BOOST_AUTO_TEST_SUITE_END()