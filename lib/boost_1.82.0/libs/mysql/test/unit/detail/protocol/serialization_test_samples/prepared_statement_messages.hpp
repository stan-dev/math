//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_PREPARED_STATEMENT_MESSAGES_HPP
#define BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_PREPARED_STATEMENT_MESSAGES_HPP

#include <boost/mysql/field_view.hpp>

#include <boost/mysql/detail/protocol/prepared_statement_messages.hpp>

#include <array>
#include <forward_list>
#include <memory>

#include "../serialization_test.hpp"

namespace boost {
namespace mysql {
namespace test {

// clang-format off
const serialization_test_spec com_stmt_prepare_packet_spec {
    serialization_test_type::serialization, {
        { "com_stmt_prepare_packet", detail::com_stmt_prepare_packet{
            string_eof("SELECT * from three_rows_table WHERE id = ?")
        }, {
            0x16, 0x53, 0x45, 0x4c, 0x45, 0x43, 0x54, 0x20,
            0x2a, 0x20, 0x66, 0x72, 0x6f, 0x6d, 0x20, 0x74,
            0x68, 0x72, 0x65, 0x65, 0x5f, 0x72, 0x6f, 0x77,
            0x73, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x20,
            0x57, 0x48, 0x45, 0x52, 0x45, 0x20, 0x69, 0x64,
            0x20, 0x3d, 0x20, 0x3f
        } }
    }
};

const serialization_test_spec com_stmt_prepare_ok_packet_spec {
    serialization_test_type::deserialization_space, {
        { "com_stmt_prepare_ok_packet", detail::com_stmt_prepare_ok_packet{
            1, // statement id
            2, // number of fields
            3, // number of params
            0, // warnings
        }, {
            0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03, 0x00,
            0x00, 0x00, 0x00
        } }
    }
};

// Helper for composing ComStmtExecute tests
template <std::size_t N, class Collection = std::vector<field_view>>
serialization_sample make_stmt_execute_sample(
    std::uint32_t stmt_id,
    std::uint8_t flags,
    std::uint32_t itercount,
    std::uint8_t new_params_flag,
    const std::array<field_view, N>& params,
    std::vector<std::uint8_t>&& buffer,
    std::string&& test_name
)
{
    auto params_shared = std::make_shared<Collection>(params.begin(), params.end());
    return serialization_sample(
        std::move(test_name),
        detail::com_stmt_execute_packet<typename Collection::const_iterator> {
            stmt_id,
            flags,
            itercount,
            new_params_flag,
            params_shared->begin(),
            params_shared->end()
        },
        std::move(buffer),
        0, // capabilities
        params_shared
    );
}

constexpr std::uint8_t com_stmt_execute_blob_buffer[] = {0x70, 0x00, 0x01, 0xff};

const serialization_test_spec com_stmt_execute_packet_spec {
    serialization_test_type::serialization, {
        make_stmt_execute_sample(1, 0x80, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(std::uint64_t(0xabffffabacadae)), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x80, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0x08, 0x80, 0xae, 0xad,
                0xac, 0xab, 0xff, 0xff, 0xab, 0x00
            },
            "uint64_t"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(std::int64_t(-0xabffffabacadae)), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0x08, 0x00, 0x52, 0x52,
                0x53, 0x54, 0x00, 0x00, 0x54, 0xff
            },
            "int64_t"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(string_view("test")), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0xfe, 0x00, 0x04, 0x74,
                0x65, 0x73, 0x74
            },
            "string"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(blob_view(com_stmt_execute_blob_buffer)), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0xfc, 0x00, 0x04, 0x70,
                0x00, 0x01, 0xff
            },
            "blob"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(3.14e20f), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0x04, 0x00, 0x01, 0x2d,
                0x88, 0x61
            },
            "float"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(2.1e214), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0x05, 0x00, 0x56, 0xc0,
                0xee, 0xa6, 0x95, 0x30, 0x6f, 0x6c
            },
            "double"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(date(2010u, 9u, 3u)), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0x0a, 0x00, 0x04, 0xda,
                0x07, 0x09, 0x03
            },
            "date"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(datetime(2010u, 9u, 3u, 10u, 30u, 59u, 231800u)), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0x0c, 0x00, 0x0b, 0xda,
                0x07, 0x09, 0x03, 0x0a, 0x1e, 0x3b, 0x78, 0x89,
                0x03, 0x00
            },
            "datetime"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(maket(230, 30, 59, 231800)), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0x0b, 0x00, 0x0c, 0x00,
                0x09, 0x00, 0x00, 0x00, 0x0e, 0x1e, 0x3b, 0x78,
                0x89, 0x03, 0x00
            },
            "time"
        ),
        make_stmt_execute_sample(1, 0, 1, 1, // stmt ID, flags, itercount, new params
            make_fv_arr(nullptr), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x01, 0x01, 0x06, 0x00
            },
            "null"
        ),
        make_stmt_execute_sample(2, 0, 1, 1,
            make_fv_arr(
                std::uint64_t(0xabffffabacadae),
                std::int64_t(-0xabffffabacadae),
                string_view("test"),
                nullptr,
                2.1e214,
                date(2010u, 9u, 3u),
                datetime(2010u, 9u, 3u, 10u, 30u, 59u, 231800u),
                maket(230, 30, 59, 231800),
                nullptr
            ), {
                0x17, 0x02, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
                0x00, 0x00, 0x08, 0x01, 0x01, 0x08, 0x80, 0x08, 0x00, 0xfe, 0x00, 0x06,
                0x00, 0x05, 0x00, 0x0a, 0x00, 0x0c, 0x00, 0x0b,
                0x00, 0x06, 0x00, 0xae, 0xad, 0xac, 0xab, 0xff,
                0xff, 0xab, 0x00, 0x52, 0x52, 0x53, 0x54, 0x00,
                0x00, 0x54, 0xff, 0x04, 0x74, 0x65, 0x73, 0x74,
                0x56, 0xc0, 0xee, 0xa6, 0x95, 0x30, 0x6f, 0x6c,
                0x04, 0xda, 0x07, 0x09, 0x03, 0x0b, 0xda, 0x07,
                0x09, 0x03, 0x0a, 0x1e, 0x3b, 0x78, 0x89, 0x03,
                0x00, 0x0c, 0x00, 0x09, 0x00, 0x00, 0x00, 0x0e,
                0x1e, 0x3b, 0x78, 0x89, 0x03, 0x00
            },
            "several_params"
        ),
        make_stmt_execute_sample<1, std::forward_list<field_view>>(1, 0x80, 1, 1,
            make_fv_arr(std::uint64_t(0xabffff)), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x80, 0x01, 0x00,
                0x00, 0x00, 0x00, 0x01, 0x08, 0x80, 0xff, 0xff,
                0xab, 0x00, 0x00, 0x00, 0x00, 0x00
            },
            "forward_list_iterator"
        ),
        make_stmt_execute_sample(1, 0x80, 1, 1,
            make_fv_arr(), {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x80, 0x01, 0x00,
                0x00, 0x00
            },
            "empty"
        )
    }
};

const serialization_test_spec com_stmt_close_packet_spec {
    serialization_test_type::serialization, {
        { "com_stmt_close_packet", detail::com_stmt_close_packet{1}, {0x19, 0x01, 0x00, 0x00, 0x00} }
    }
};
// clang-format on

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
