//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/impl/internal/protocol/impl/protocol_types.hpp>

#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <memory>

#include "operators.hpp"
#include "serialization_test.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_common/create_basic.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::string_view;

BOOST_AUTO_TEST_SUITE(test_protocol_types)

template <std::size_t N>
const char* get_string_N()
{
    static std::string res(N, 'a');
    return res.c_str();
}

template <std::size_t N>
string_fixed<N> makesfixed(const char (&value)[N + 1])
{
    static_assert(N >= 1, "Expected a C-array literal");
    string_fixed<N> res;
    std::memcpy(res.value.data(), value, N);
    return res;
}

class test_case
{
    const char* name_;

public:
    test_case(const char* name) noexcept : name_(name) {}
    virtual ~test_case() {}
    virtual void serialize_test() = 0;
    virtual void deserialize_test() = 0;
    virtual void deserialize_space_test() = 0;
    virtual void deserialize_not_enough_space_test() = 0;
    const char* name() const noexcept { return name_; }
};

template <class T>
class test_case_impl final : public test_case
{
    T value_;
    std::vector<std::uint8_t> serialized_;
    bool do_space_;

public:
    test_case_impl(const char* name, T value, std::vector<std::uint8_t> serialized, bool do_space = true)
        : test_case(name), value_(value), serialized_(std::move(serialized)), do_space_(do_space)
    {
    }
    void serialize_test() override final { do_serialize_test(value_, serialized_); }
    void deserialize_test() override final { do_deserialize_test(value_, serialized_); }
    void deserialize_space_test() override final
    {
        if (do_space_)
            do_deserialize_extra_space_test(value_, serialized_);
    }
    void deserialize_not_enough_space_test() override final
    {
        if (do_space_)
            do_deserialize_not_enough_space_test<T>(serialized_);
    }
};

template <class T>
std::shared_ptr<test_case> make_test(
    const char* name,
    T value,
    std::vector<std::uint8_t> serialized,
    bool do_space = true
)
{
    return std::make_shared<test_case_impl<T>>(name, value, std::move(serialized), do_space);
}

std::vector<std::shared_ptr<test_case>> make_all_cases()
{
    return {
        // basic integers
        make_test("int1", int1{0xff}, {0xff}),
        make_test("int2", int2{0xfeff}, {0xff, 0xfe}),
        make_test("int3", int3{0xfdfeff}, {0xff, 0xfe, 0xfd}),
        make_test("int4", int4{0xfcfdfeff}, {0xff, 0xfe, 0xfd, 0xfc}),
        make_test("int8", int8{0xf8f9fafbfcfdfeff}, {0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8}),
        make_test("sint1_positive", int_holder<std::int8_t>{0x01}, {0x01}),
        make_test("sint1_negative", int_holder<std::int8_t>{-1}, {0xff}),
        make_test("sint2_positive", int_holder<std::int16_t>{0x0201}, {0x01, 0x02}),
        make_test("sint2_negative", int_holder<std::int16_t>{-0x101}, {0xff, 0xfe}),
        make_test("sint4_positive", int_holder<std::int32_t>{0x04030201}, {0x01, 0x02, 0x03, 0x04}),
        make_test("sint4_negative", int_holder<std::int32_t>{-0x3020101}, {0xff, 0xfe, 0xfd, 0xfc}),
        make_test(
            "sint8_positive",
            sint8{0x0807060504030201},
            {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08}
        ),
        make_test(
            "sint8_negative",
            sint8{-0x0706050403020101},
            {0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8}
        ),

        // int_lenenc
        make_test("int_lenenc_1_byte_regular", int_lenenc{1}, {0x01}),
        make_test("int_lenenc_1_byte_max", int_lenenc{250}, {0xfa}),
        make_test("int_lenenc_2_bytes_regular", int_lenenc{0xfeb7}, {0xfc, 0xb7, 0xfe}),
        make_test("int_lenenc_2_bytes_max", int_lenenc{0xffff}, {0xfc, 0xff, 0xff}),
        make_test("int_lenenc_3_bytes_regular", int_lenenc{0xa0feff}, {0xfd, 0xff, 0xfe, 0xa0}),
        make_test("int_lenenc_3_bytes_max", int_lenenc{0xffffff}, {0xfd, 0xff, 0xff, 0xff}),
        make_test(
            "int_lenenc_8_bytes_regular",
            int_lenenc{0xf8f9fafbfcfdfeff},
            {0xfe, 0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8}
        ),
        make_test(
            "int_lenenc_8_bytes_max",
            int_lenenc{0xffffffffffffffff},
            {0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}
        ),

        // string_fixed
        make_test("4c_regular_characters", makesfixed<4>("abde"), {0x61, 0x62, 0x64, 0x65}),
        make_test("3c_null_characters", makesfixed<3>("\0\1a"), {0x00, 0x01, 0x61}),
        make_test("3c_utf8_characters", string_fixed<3>{{{'\xc3', '\xb1', 'a'}}}, {0xc3, 0xb1, 0x61}),
        make_test("1c_regular_characters", makesfixed<1>("a"), {0x61}),

        // string_null
        make_test("regular_characters", string_null{"abc"}, {0x61, 0x62, 0x63, 0x00}),
        make_test("utf8_characters", string_null{"\xc3\xb1"}, {0xc3, 0xb1, 0x00}),
        make_test("empty", string_null{""}, {0x00}),

        // string_lenenc
        make_test("empty", string_lenenc{""}, {0x00}),
        make_test("1_byte_size_regular_characters", string_lenenc{"abc"}, {0x03, 0x61, 0x62, 0x63}),
        make_test("1_byte_size_null_characters", string_lenenc{makesv("a\0b")}, {0x03, 0x61, 0x00, 0x62}),
        make_test(
            "1_byte_size_max",
            string_lenenc{get_string_N<250>()},
            concat({250}, std::vector<std::uint8_t>(250, 0x61))
        ),
        make_test(
            "2_byte_size_min",
            string_lenenc{get_string_N<251>()},
            concat({0xfc, 251, 0}, std::vector<std::uint8_t>(251, 0x61))
        ),
        make_test(
            "2_byte_size_max",
            string_lenenc{get_string_N<0xffff>()},
            concat({0xfc, 0xff, 0xff}, std::vector<std::uint8_t>(0xffff, 0x61))
        ),
        make_test(
            "3_byte_size_min",
            string_lenenc{get_string_N<0x10000>()},
            concat({0xfd, 0x00, 0x00, 0x01}, std::vector<std::uint8_t>(0x10000, 0x61))
        ),
    };
}

std::vector<std::shared_ptr<test_case>> all_cases()
{
    static auto res = make_all_cases();
    return res;
}

BOOST_AUTO_TEST_CASE(serialize)
{
    for (const auto& sample : all_cases())
    {
        BOOST_TEST_CONTEXT(sample->name()) { sample->serialize_test(); }
    }
}

BOOST_AUTO_TEST_CASE(deserialize)
{
    for (const auto& sample : all_cases())
    {
        BOOST_TEST_CONTEXT(sample->name()) { sample->deserialize_test(); }
    }
}

BOOST_AUTO_TEST_CASE(deserialize_extra_space)
{
    for (const auto& sample : all_cases())
    {
        BOOST_TEST_CONTEXT(sample->name()) { sample->deserialize_space_test(); }
    }
}
BOOST_AUTO_TEST_CASE(deserialize_not_enough_space)
{
    for (const auto& sample : all_cases())
    {
        BOOST_TEST_CONTEXT(sample->name()) { sample->deserialize_not_enough_space_test(); }
    }
}

// string_eof can be serialized/deserialized, but is space sensitive, so extra space/not enough space tests
// shouldn't be performed
BOOST_AUTO_TEST_CASE(string_eof_)
{
    struct
    {
        const char* name;
        string_eof value;
        std::vector<std::uint8_t> serialized;
    } test_cases[] = {
        {"regular_characters", string_eof{"abc"},          {0x61, 0x62, 0x63}},
        {"null_characters",    string_eof{makesv("a\0b")}, {0x61, 0x00, 0x62}},
        {"empty",              string_eof{""},             {}                },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name << ", serialization") { do_serialize_test(tc.value, tc.serialized); }
        BOOST_TEST_CONTEXT(tc.name << ", deserialization") { do_deserialize_test(tc.value, tc.serialized); }
    }
}

BOOST_AUTO_TEST_SUITE_END()
