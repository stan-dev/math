//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_COMMON_MESSAGES_HPP
#define BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_COMMON_MESSAGES_HPP

#include <boost/mysql/mysql_collations.hpp>

#include <boost/mysql/detail/protocol/common_messages.hpp>

#include "../serialization_test.hpp"

namespace boost {
namespace mysql {
namespace test {

// clang-format off
const serialization_test_spec packet_header_spec {
    serialization_test_type::full, {
        { "small_packet_seqnum_0",
            detail::packet_header{int3(3), 0}, {0x03, 0x00, 0x00, 0x00} },
        { "small_packet_seqnum_not_0",
            detail::packet_header{int3(9), 2}, {0x09, 0x00, 0x00, 0x02} },
        { "big_packet_seqnum_0",
            detail::packet_header{int3(0xcacbcc), 0xfa}, {0xcc, 0xcb, 0xca, 0xfa} },
        { "max_packet_max_seqnum",
            detail::packet_header{int3(0xffffff), 0xff}, {0xff, 0xff, 0xff, 0xff} }
    }
};

const serialization_test_spec ok_packet_spec {
    serialization_test_type::deserialization, {
        { "successful_update", detail::ok_packet{
            int_lenenc(4), // affected rows
            int_lenenc(0), // last insert ID
            static_cast<std::uint16_t>(
                detail::SERVER_STATUS_AUTOCOMMIT |
                detail::SERVER_QUERY_NO_INDEX_USED
            ), // server status
            0, // warnings
            string_lenenc("Rows matched: 5  Changed: 4  Warnings: 0")
        }, {
            0x04, 0x00, 0x22, 0x00, 0x00, 0x00, 0x28, 0x52, 0x6f, 0x77, 0x73,
            0x20, 0x6d, 0x61, 0x74, 0x63, 0x68, 0x65, 0x64, 0x3a, 0x20, 0x35, 0x20, 0x20, 0x43, 0x68, 0x61,
            0x6e, 0x67, 0x65, 0x64, 0x3a, 0x20, 0x34, 0x20, 0x20, 0x57, 0x61, 0x72, 0x6e, 0x69, 0x6e, 0x67,
            0x73, 0x3a, 0x20, 0x30
        } },

        { "successful_insert", detail::ok_packet{
            int_lenenc(1), // affected rows
            int_lenenc(6), // last insert ID
            static_cast<std::uint16_t>(detail::SERVER_STATUS_AUTOCOMMIT), // server status
            0, // warnings
            string_lenenc("")  // no message
        },{
            0x01, 0x06, 0x02, 0x00, 0x00, 0x00
        } },

        { "successful_login", detail::ok_packet{
            int_lenenc(0), // affected rows
            int_lenenc(0), // last insert ID
            static_cast<std::uint16_t>(detail::SERVER_STATUS_AUTOCOMMIT), // server status
            0, // warnings
            string_lenenc("")  // no message
        }, {
            0x00, 0x00, 0x02, 0x00, 0x00, 0x00
        } }
    }
};

const serialization_test_spec err_packet_spec {
    serialization_test_type::deserialization, {
        { "wrong_use_database", detail::err_packet{
            1049, // eror code
            string_fixed<1>{{0x23}}, // sql state marker
            makesfixed<5>("42000"), // sql state
            string_eof("Unknown database 'a'") // err msg
        }, {
            0x19, 0x04, 0x23, 0x34, 0x32, 0x30, 0x30, 0x30, 0x55, 0x6e, 0x6b,
            0x6e, 0x6f, 0x77, 0x6e, 0x20, 0x64, 0x61, 0x74,
            0x61, 0x62, 0x61, 0x73, 0x65, 0x20, 0x27, 0x61, 0x27
        } },

        { "unknown_table", detail::err_packet{
            1146, // eror code
            string_fixed<1>{{0x23}}, // sql state marker
            makesfixed<5>("42S02"), // sql state
            string_eof("Table 'awesome.unknown' doesn't exist") // err msg
        }, {
            0x7a, 0x04, 0x23, 0x34, 0x32, 0x53, 0x30, 0x32,
            0x54, 0x61, 0x62, 0x6c, 0x65, 0x20, 0x27, 0x61,
            0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x2e, 0x75,
            0x6e, 0x6b, 0x6e, 0x6f, 0x77, 0x6e, 0x27, 0x20,
            0x64, 0x6f, 0x65, 0x73, 0x6e, 0x27, 0x74, 0x20,
            0x65, 0x78, 0x69, 0x73, 0x74
        } },

        { "failed_login", detail::err_packet{
            1045, // error code
            string_fixed<1>{{0x23}}, // SQL state marker
            makesfixed<5>("28000"), // sql state
            string_eof("Access denied for user 'root'@'localhost' (using password: YES)")
        }, {
          0x15, 0x04, 0x23, 0x32, 0x38, 0x30, 0x30, 0x30,
          0x41, 0x63, 0x63, 0x65, 0x73, 0x73, 0x20, 0x64,
          0x65, 0x6e, 0x69, 0x65, 0x64, 0x20, 0x66, 0x6f,
          0x72, 0x20, 0x75, 0x73, 0x65, 0x72, 0x20, 0x27,
          0x72, 0x6f, 0x6f, 0x74, 0x27, 0x40, 0x27, 0x6c,
          0x6f, 0x63, 0x61, 0x6c, 0x68, 0x6f, 0x73, 0x74,
          0x27, 0x20, 0x28, 0x75, 0x73, 0x69, 0x6e, 0x67,
          0x20, 0x70, 0x61, 0x73, 0x73, 0x77, 0x6f, 0x72,
          0x64, 0x3a, 0x20, 0x59, 0x45, 0x53, 0x29
        } }
    }
};

const serialization_test_spec column_definition_spec {
    serialization_test_type::deserialization_space, {
        { "numeric_auto_increment_primary_key", detail::column_definition_packet{
            string_lenenc("def"), //catalog
            string_lenenc("awesome"), // schema (database)
            string_lenenc("test_table"), // table
            string_lenenc("test_table"), // physical table
            string_lenenc("id"), // field name
            string_lenenc("id"), // physical field name
            mysql_collations::binary,
            11, // length
            detail::protocol_field_type::long_,
            detail::column_flags::not_null | detail::column_flags::pri_key |
                    detail::column_flags::auto_increment | detail::column_flags::part_key,
            0 // decimals
        }, {
            0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65,
            0x73, 0x6f, 0x6d, 0x65, 0x0a, 0x74, 0x65, 0x73,
            0x74, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0a,
            0x74, 0x65, 0x73, 0x74, 0x5f, 0x74, 0x61, 0x62,
            0x6c, 0x65, 0x02, 0x69, 0x64, 0x02, 0x69, 0x64,
            0x0c, 0x3f, 0x00, 0x0b, 0x00, 0x00, 0x00, 0x03,
            0x03, 0x42, 0x00, 0x00, 0x00
        } },
        { "varchar_field_aliased_field_and_table_names_join", detail::column_definition_packet{
            string_lenenc("def"), //catalog
            string_lenenc("awesome"), // schema (database)
            string_lenenc("child"), // table
            string_lenenc("child_table"), // physical table
            string_lenenc("field_alias"), // field name
            string_lenenc("field_varchar"), // physical field name
            mysql_collations::utf8_general_ci,
            765, // length
            detail::protocol_field_type::var_string,
            0, // no column flags
            0 // decimals
        }, {
            0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65,
            0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63, 0x68, 0x69,
            0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64,
            0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0b, 0x66,
            0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c, 0x69,
            0x61, 0x73, 0x0d, 0x66, 0x69, 0x65, 0x6c, 0x64,
            0x5f, 0x76, 0x61, 0x72, 0x63, 0x68, 0x61, 0x72,
            0x0c, 0x21, 0x00, 0xfd, 0x02, 0x00, 0x00, 0xfd,
            0x00, 0x00, 0x00, 0x00, 0x00
        } },

        { "float_field", detail::column_definition_packet{
            string_lenenc("def"), //catalog
            string_lenenc("awesome"), // schema (database)
            string_lenenc("test_table"), // table
            string_lenenc("test_table"), // physical table
            string_lenenc("field_float"), // field name
            string_lenenc("field_float"), // physical field name
            mysql_collations::binary, // binary
            12, // length
            detail::protocol_field_type::float_,
            0, // no column flags
            31 // decimals
        }, {
            0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65,
            0x73, 0x6f, 0x6d, 0x65, 0x0a, 0x74, 0x65, 0x73,
            0x74, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0a,
            0x74, 0x65, 0x73, 0x74, 0x5f, 0x74, 0x61, 0x62,
            0x6c, 0x65, 0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64,
            0x5f, 0x66, 0x6c, 0x6f, 0x61, 0x74, 0x0b, 0x66,
            0x69, 0x65, 0x6c, 0x64, 0x5f, 0x66, 0x6c, 0x6f,
            0x61, 0x74, 0x0c, 0x3f, 0x00, 0x0c, 0x00, 0x00,
            0x00, 0x04, 0x00, 0x00, 0x1f, 0x00, 0x00
        } }
    }
};

const serialization_test_spec quit_packet_spec {
    serialization_test_type::serialization, {
        { "quit_packet", detail::quit_packet(), {0x01} }
    }
};

const serialization_test_spec ping_packet_spec {
    serialization_test_type::serialization, {
        { "ping_packet", detail::ping_packet(), {0x0e} }
    }
};
// clang-format on

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
