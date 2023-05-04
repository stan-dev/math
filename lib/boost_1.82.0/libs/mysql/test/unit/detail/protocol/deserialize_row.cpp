//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Tests for both deserialize_binary_row() and deserialize_text_row()

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata.hpp>

#include <boost/mysql/detail/auxiliar/string_view_offset.hpp>
#include <boost/mysql/detail/protocol/capabilities.hpp>
#include <boost/mysql/detail/protocol/common_messages.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>
#include <boost/mysql/detail/protocol/db_flavor.hpp>
#include <boost/mysql/detail/protocol/deserialize_row.hpp>
#include <boost/mysql/detail/protocol/resultset_encoding.hpp>

#include <boost/asio/buffer.hpp>
#include <boost/test/tools/context.hpp>
#include <boost/test/unit_test.hpp>

#include <cstdint>

#include "buffer_concat.hpp"
#include "create_execution_state.hpp"
#include "create_meta.hpp"
#include "test_common.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql::detail;
using boost::mysql::client_errc;
using boost::mysql::common_server_errc;
using boost::mysql::date;
using boost::mysql::diagnostics;
using boost::mysql::error_code;
using boost::mysql::field_view;
using boost::mysql::metadata;

namespace {

BOOST_AUTO_TEST_SUITE(test_deserialize_row)

BOOST_AUTO_TEST_SUITE(without_state)

std::vector<metadata> make_meta(const std::vector<protocol_field_type>& types)
{
    std::vector<metadata> res;
    for (const auto type : types)
    {
        res.push_back(create_meta(type));
    }
    return res;
}

BOOST_AUTO_TEST_CASE(success)
{
    // clang-format off
    struct
    {
        const char* name;
        resultset_encoding encoding;
        std::vector<std::uint8_t> from;
        std::vector<field_view> expected_with_offsets;  // before offset conversion
        std::vector<field_view> expected;               // after offset conversion
        std::vector<metadata> meta;
    } test_cases [] = {
        // Text
        {
            "text_one_value",
            resultset_encoding::text,
            {0x01, 0x35},
            make_fv_vector(std::int64_t(5)),
            make_fv_vector(std::int64_t(5)),
            make_meta({ protocol_field_type::tiny })
        },
        {
            "text_one_null",
            resultset_encoding::text,
            {0xfb},
            make_fv_vector(nullptr),
            make_fv_vector(nullptr),
            make_meta({ protocol_field_type::tiny })
        },
        {
            "text_several_values",
            resultset_encoding::text,
            {0x03, 0x76, 0x61, 0x6c, 0x02, 0x32, 0x31, 0x03, 0x30, 0x2e, 0x30},
            make_fv_vector(make_svoff_fv(1, 3, false), std::int64_t(21), 0.0f),
            make_fv_vector("val", std::int64_t(21), 0.0f),
            make_meta({ protocol_field_type::var_string, protocol_field_type::long_, protocol_field_type::float_ })
        },
        {
            "text_several_values_one_null",
            resultset_encoding::text,
            {0x03, 0x76, 0x61, 0x6c, 0xfb, 0x03, 0x76, 0x61, 0x6c},
            make_fv_vector(make_svoff_fv(1, 3, false), nullptr, make_svoff_fv(6, 3, false)),
            make_fv_vector("val", nullptr, "val"),
            make_meta({ protocol_field_type::var_string, protocol_field_type::long_, protocol_field_type::var_string })
        },
        {
            "text_several_nulls",
            resultset_encoding::text,
            {0xfb, 0xfb, 0xfb},
            make_fv_vector(nullptr, nullptr, nullptr),
            make_fv_vector(nullptr, nullptr, nullptr),
            make_meta({ protocol_field_type::var_string, protocol_field_type::long_, protocol_field_type::datetime })
        },

        // Binary
        {
            "binary_one_value",
            resultset_encoding::binary,
            {0x00, 0x00, 0x14},
            make_fv_vector(std::int64_t(20)),
            make_fv_vector(std::int64_t(20)),
            make_meta({ protocol_field_type::tiny })
        },
        {
            "binary_one_null",
            resultset_encoding::binary,
            {0x00, 0x04},
            make_fv_vector(nullptr),
            make_fv_vector(nullptr),
            make_meta({ protocol_field_type::tiny })
        },
        {
            "binary_two_values",
            resultset_encoding::binary,
            {0x00, 0x00, 0x03, 0x6d, 0x69, 0x6e, 0x6d, 0x07},
            make_fv_vector(make_svoff_fv(3, 3, false), std::int64_t(1901)),
            make_fv_vector("min", std::int64_t(1901)),
            make_meta({ protocol_field_type::var_string, protocol_field_type::short_ })
        },
        {
            "binary_one_value_one_null",
            resultset_encoding::binary,
            {0x00, 0x08, 0x03, 0x6d, 0x61, 0x78},
            make_fv_vector(make_svoff_fv(3, 3, false), nullptr),
            make_fv_vector("max", nullptr),
            make_meta({ protocol_field_type::var_string, protocol_field_type::tiny })
        },
        {
            "binary_two_nulls",
            resultset_encoding::binary,
            {0x00, 0x0c},
            make_fv_vector(nullptr, nullptr),
            make_fv_vector(nullptr, nullptr),
            make_meta({ protocol_field_type::tiny, protocol_field_type::tiny })
        },
        {
            "binary_six_nulls",
            resultset_encoding::binary,
            {0x00, 0xfc},
            std::vector<field_view>(6, field_view()),
            std::vector<field_view>(6, field_view()),
            make_meta(std::vector<protocol_field_type>(6, protocol_field_type::tiny))
        },
        {
            "binary_seven_nulls",
            resultset_encoding::binary,
            {0x00, 0xfc, 0x01},
            std::vector<field_view>(7, field_view()),
            std::vector<field_view>(7, field_view()),
            make_meta(std::vector<protocol_field_type>(7, protocol_field_type::tiny))
        },
        {
            "binary_several_values",
            resultset_encoding::binary,
            {
                0x00, 0x90, 0x00, 0xfd, 0x03, 0x61, 0x62, 0x63,
                0xc3, 0xf5, 0x48, 0x40, 0x02, 0x61, 0x62, 0x04,
                0xe2, 0x07, 0x0a, 0x05, 0x71, 0x99, 0x6d, 0xe2,
                0x93, 0x4d, 0xf5, 0x3d
            },
            make_fv_vector(
                std::int64_t(-3),
                make_svoff_fv(5, 3, false),
                nullptr,
                3.14f,
                make_svoff_fv(13, 2, false),
                nullptr,
                date(2018u, 10u, 5u),
                3.10e-10
            ),
            make_fv_vector(
                std::int64_t(-3),
                "abc",
                nullptr,
                3.14f,
                "ab",
                nullptr,
                date(2018u, 10u, 5u),
                3.10e-10
            ),
            make_meta({
                protocol_field_type::tiny,
                protocol_field_type::var_string,
                protocol_field_type::long_,
                protocol_field_type::float_,
                protocol_field_type::string,
                protocol_field_type::long_,
                protocol_field_type::date,
                protocol_field_type::double_
            })
        }
    };
    // clang-format on

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            const auto& buffer = tc.from;
            deserialization_context ctx(buffer.data(), buffer.data() + buffer.size(), capabilities());
            std::vector<field_view> actual;
            error_code err;

            deserialize_row(tc.encoding, ctx, tc.meta, tc.from.data(), actual, err);
            BOOST_TEST(err == error_code());
            BOOST_TEST(actual == tc.expected_with_offsets);

            offsets_to_string_views(actual, tc.from.data());
            BOOST_TEST(actual == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(error)
{
    // clang-format off
    struct
    {
        const char* name;
        resultset_encoding encoding;
        std::vector<std::uint8_t> from;
        client_errc expected;
        std::vector<metadata> meta;
    } test_cases [] = {
        // text
        {
            "text_no_space_string_single",
            resultset_encoding::text,
            {0x02, 0x00},
            client_errc::incomplete_message,
            make_meta({ protocol_field_type::short_} )
        },
        {
            "text_no_space_string_final",
            resultset_encoding::text,
            {0x01, 0x35, 0x02, 0x35},
            client_errc::incomplete_message,
            make_meta({ protocol_field_type::tiny, protocol_field_type::short_ }),
        },
        {
            "text_no_space_null_single",
            resultset_encoding::text,
            {},
            client_errc::incomplete_message,
            make_meta({ protocol_field_type::tiny })
        },
        {
            "text_no_space_null_final",
            resultset_encoding::text,
            {0xfb},
            client_errc::incomplete_message,
            make_meta({protocol_field_type::tiny, protocol_field_type::tiny}),
        },
        {
            "text_extra_bytes",
            resultset_encoding::text,
            {0x01, 0x35, 0xfb, 0x00},
            client_errc::extra_bytes,
            make_meta({ protocol_field_type::tiny, protocol_field_type::tiny })
        },
        {
            "text_contained_value_error_single",
            resultset_encoding::text,
            {0x01, 0x00},
            client_errc::protocol_value_error,
            make_meta({ protocol_field_type::date })
        },
        {
            "text_contained_value_error_middle",
            resultset_encoding::text,
            {0xfb, 0x01, 0x00, 0xfb},
            client_errc::protocol_value_error,
            make_meta({ protocol_field_type::date, protocol_field_type::date, protocol_field_type::date })
        },
        {
            "text_row_for_empty_meta",
            resultset_encoding::text,
            {0xfb, 0x01, 0x00, 0xfb},
            client_errc::extra_bytes,
            make_meta({})
        },

        // binary
        {
            "binary_no_space_null_bitmap_1",
            resultset_encoding::binary,
            {0x00},
            client_errc::incomplete_message,
            make_meta({ protocol_field_type::tiny })
        },
        {
            "binary_no_space_null_bitmap_2",
            resultset_encoding::binary,
            {0x00, 0xfc},
            client_errc::incomplete_message,
            make_meta(std::vector<protocol_field_type>(7, protocol_field_type::tiny))
        },
        {
            "binary_no_space_value_single",
            resultset_encoding::binary,
            {0x00, 0x00},
            client_errc::incomplete_message,
            make_meta({ protocol_field_type::tiny })
        },
        {
            "binary_no_space_value_last",
            resultset_encoding::binary,
            {0x00, 0x00, 0x01},
            client_errc::incomplete_message,
            make_meta(std::vector<protocol_field_type>(2, protocol_field_type::tiny))
        },
        {
            "binary_no_space_value_middle",
            resultset_encoding::binary,
            {0x00, 0x00, 0x01},
            client_errc::incomplete_message,
            make_meta(std::vector<protocol_field_type>(3, protocol_field_type::tiny))
        },
        {
            "binary_extra_bytes",
            resultset_encoding::binary,
            {0x00, 0x00, 0x01, 0x02},
            client_errc::extra_bytes,
            make_meta({ protocol_field_type::tiny })
        },
        {
            "binary_row_for_empty_meta",
            resultset_encoding::binary,
            {0xfb, 0x01, 0x00, 0xfb},
            client_errc::extra_bytes,
            make_meta({})
        },
    };
    // clang-format on

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            const auto& buffer = tc.from;
            deserialization_context ctx(buffer.data(), buffer.data() + buffer.size(), capabilities());
            std::vector<field_view> actual;
            error_code err;

            deserialize_row(tc.encoding, ctx, tc.meta, buffer.data(), actual, err);
            BOOST_TEST(err == error_code(tc.expected));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(with_execution_state)

BOOST_AUTO_TEST_CASE(text_rows)
{
    std::vector<std::uint8_t> row1{0x03, 0x76, 0x61, 0x6c, 0x02, 0x32, 0x31, 0x03, 0x30, 0x2e, 0x30};
    std::vector<std::uint8_t> row2{0x03, 0x61, 0x62, 0x63, 0x02, 0x32, 0x30, 0x03, 0x30, 0x2e, 0x30};
    auto buff = concat_copy(row1, row2);
    auto st = create_execution_state(
        resultset_encoding::text,
        {protocol_field_type::var_string, protocol_field_type::long_, protocol_field_type::float_}
    );
    auto expected_fields = make_fv_vector(make_svoff_fv(1, 3, false), std::int64_t(21), 0.0f);
    std::vector<field_view> fields;
    error_code err;
    diagnostics diag;

    // First row
    deserialize_row(
        boost::asio::const_buffer(buff.data(), row1.size()),
        capabilities(),
        db_flavor::mysql,
        buff.data(),
        st,
        fields,
        err,
        diag
    );

    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.server_message() == "");
    BOOST_TEST(!st.complete());
    BOOST_TEST(fields == expected_fields);

    // Second row (fields get appended to existing ones)
    deserialize_row(
        boost::asio::const_buffer(buff.data() + row1.size(), row2.size()),
        capabilities(),
        db_flavor::mysql,
        buff.data(),
        st,
        fields,
        err,
        diag
    );
    expected_fields.emplace_back(make_svoff_fv(12, 3, false));
    expected_fields.emplace_back(20);
    expected_fields.emplace_back(0.0f);

    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.server_message() == "");
    BOOST_TEST(!st.complete());
    BOOST_TEST(fields == expected_fields);

    // Convert offsets to string views
    offsets_to_string_views(fields, buff.data());
    BOOST_TEST(fields == make_fv_vector("val", 21, 0.0f, "abc", 20, 0.0f));
}

BOOST_AUTO_TEST_CASE(binary_rows)
{
    std::vector<std::uint8_t> row1{0x00, 0x00, 0x03, 0x6d, 0x69, 0x6e, 0x6d, 0x07};
    std::vector<std::uint8_t> row2{0x00, 0x08, 0x03, 0x6d, 0x61, 0x78};
    auto buff = concat_copy(row1, row2);
    auto st = create_execution_state(
        resultset_encoding::binary,
        {protocol_field_type::var_string, protocol_field_type::short_}
    );
    auto expected_fields = make_fv_vector(make_svoff_fv(3, 3, false), std::int64_t(1901));
    std::vector<field_view> fields;
    error_code err;
    diagnostics diag;

    // First row
    deserialize_row(
        boost::asio::const_buffer(buff.data(), row1.size()),
        capabilities(),
        db_flavor::mysql,
        buff.data(),
        st,
        fields,
        err,
        diag
    );

    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.server_message() == "");
    BOOST_TEST(!st.complete());
    BOOST_TEST(fields == expected_fields);

    // Second row (fields get appended to existing ones)
    deserialize_row(
        boost::asio::const_buffer(buff.data() + row1.size(), row2.size()),
        capabilities(),
        db_flavor::mysql,
        buff.data(),
        st,
        fields,
        err,
        diag
    );
    expected_fields.emplace_back(make_svoff_fv(11, 3, false));
    expected_fields.emplace_back(nullptr);

    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.server_message() == "");
    BOOST_TEST(!st.complete());
    BOOST_TEST(fields == expected_fields);

    // Convert offsets to string views
    offsets_to_string_views(fields, buff.data());
    BOOST_TEST(fields == make_fv_vector("min", 1901, "max", nullptr));
}

BOOST_AUTO_TEST_CASE(ok_packet)
{
    std::vector<std::uint8_t> buff{0xfe, 0x01, 0x06, 0x02, 0x00, 0x09, 0x00, 0x02, 0x61, 0x62};
    auto st = create_execution_state(
        resultset_encoding::binary,
        {protocol_field_type::var_string, protocol_field_type::short_}
    );
    auto fields_before = make_fv_vector("abc", 20);  // previous row
    auto fields = fields_before;
    error_code err;
    diagnostics diag;

    // First row
    deserialize_row(
        boost::asio::buffer(buff),
        capabilities(),
        db_flavor::mysql,
        buff.data(),
        st,
        fields,
        err,
        diag
    );

    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.server_message() == "");
    BOOST_TEST(st.complete());
    BOOST_TEST(st.affected_rows() == 1u);
    BOOST_TEST(st.last_insert_id() == 6u);
    BOOST_TEST(st.warning_count() == 9u);
    BOOST_TEST(st.info() == "ab");
    BOOST_TEST(fields == fields_before);  // they didn't change
}

BOOST_AUTO_TEST_CASE(error)
{
    // clang-format off
    struct
    {
        const char* name;
        std::vector<std::uint8_t> buffer;
        error_code expected_error;
        const char* expected_info;
    } test_cases [] = {
        {
            "invalid_row",
            { 0x00, 0x00, 0x03, 0x6d, 0x69, 0x6e, 0x6d, }, // 1 byte missing
            client_errc::incomplete_message,
            ""
        },
        {
            "invalid_ok_packet",
            { 0xfe, 0x00, 0x00, 0x02, 0x00, 0x00 }, // 1 byte missing
            client_errc::incomplete_message,
            ""
        },
        {
            "error_packet",
            {
                0xff, 0x19, 0x04, 0x23, 0x34, 0x32, 0x30, 0x30, 0x30, 0x55, 0x6e, 0x6b,
                0x6e, 0x6f, 0x77, 0x6e, 0x20, 0x64, 0x61, 0x74,
                0x61, 0x62, 0x61, 0x73, 0x65, 0x20, 0x27, 0x61, 0x27
            },
            common_server_errc::er_bad_db_error,
            "Unknown database 'a'"
        },
        {
            "invalid_error_packet",
            { 0xff, 0x19 }, // bytes missing
            client_errc::incomplete_message,
            ""
        },
        {
            "empty_message",
            {},
            client_errc::incomplete_message,
            ""
        }
    };
    // clang-format on

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_}
            );
            std::vector<field_view> fields;
            error_code err;
            diagnostics diag;

            // First row
            deserialize_row(
                boost::asio::buffer(tc.buffer),
                capabilities(),
                db_flavor::mysql,
                tc.buffer.data(),
                st,
                fields,
                err,
                diag
            );

            BOOST_TEST(err == tc.expected_error);
            BOOST_TEST(diag.server_message() == tc.expected_info);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
