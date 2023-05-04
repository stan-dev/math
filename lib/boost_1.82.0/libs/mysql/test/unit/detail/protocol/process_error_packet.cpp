//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/mariadb_server_errc.hpp>
#include <boost/mysql/mysql_server_errc.hpp>

#include <boost/mysql/detail/protocol/capabilities.hpp>
#include <boost/mysql/detail/protocol/db_flavor.hpp>
#include <boost/mysql/detail/protocol/deserialization_context.hpp>
#include <boost/mysql/detail/protocol/process_error_packet.hpp>

#include <boost/test/unit_test.hpp>

#include "create_message.hpp"

using namespace boost::mysql::detail;
using boost::mysql::client_errc;
using boost::mysql::common_server_errc;
using boost::mysql::error_code;
using boost::mysql::get_mariadb_server_category;
using boost::mysql::get_mysql_server_category;
using boost::mysql::test::create_err_packet_body;

namespace {

BOOST_AUTO_TEST_SUITE(test_process_error_packet)

BOOST_AUTO_TEST_CASE(regular)
{
    struct
    {
        const char* name;
        db_flavor flavor;
        std::vector<std::uint8_t> buffer;
        error_code ec;
        const char* msg;
    } test_cases[] = {
        {"bad_error_packet", db_flavor::mariadb, {0xff, 0x00, 0x01}, client_errc::incomplete_message, ""},
        {"code_lt_min",
         db_flavor::mariadb,
         create_err_packet_body(999, "abc"),
         error_code(999, get_mariadb_server_category()),
         "abc"},
        {"code_common",
         db_flavor::mariadb,
         create_err_packet_body(1064, "abc"),
         common_server_errc::er_parse_error,
         "abc"},
        {"code_common_hole_mysql",
         db_flavor::mysql,
         create_err_packet_body(1076),
         error_code(1076, get_mysql_server_category()),
         ""},
        {"code_common_hole_mariadb",
         db_flavor::mariadb,
         create_err_packet_body(1076),
         error_code(1076, get_mariadb_server_category()),
         ""},
        {"code_mysql",
         db_flavor::mysql,
         create_err_packet_body(4004),
         error_code(4004, get_mysql_server_category()),
         ""},
        {"code_mariadb",
         db_flavor::mariadb,
         create_err_packet_body(4004),
         error_code(4004, get_mariadb_server_category()),
         ""},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            deserialization_context ctx{
                tc.buffer.data(),
                tc.buffer.data() + tc.buffer.size(),
                capabilities()};
            ctx.advance(1);  // the generating functions include a 0xff header
            boost::mysql::diagnostics diag;
            auto ec = process_error_packet(ctx, tc.flavor, diag);
            BOOST_TEST(ec == tc.ec);
            BOOST_TEST(diag.server_message() == tc.msg);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
