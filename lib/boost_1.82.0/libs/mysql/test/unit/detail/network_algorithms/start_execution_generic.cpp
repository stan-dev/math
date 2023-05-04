//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/row.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/network_algorithms/start_execution_generic.hpp>
#include <boost/mysql/detail/protocol/capabilities.hpp>
#include <boost/mysql/detail/protocol/common_messages.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>
#include <boost/mysql/detail/protocol/db_flavor.hpp>
#include <boost/mysql/detail/protocol/protocol_types.hpp>
#include <boost/mysql/detail/protocol/resultset_encoding.hpp>

#include <boost/asio/bind_executor.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <functional>
#include <iterator>

#include "assert_buffer_equals.hpp"
#include "create_execution_state.hpp"
#include "create_message.hpp"
#include "network_result.hpp"
#include "printing.hpp"
#include "test_channel.hpp"
#include "test_common.hpp"
#include "test_stream.hpp"
#include "unit_netfun_maker.hpp"

using boost::mysql::blob;
using boost::mysql::client_errc;
using boost::mysql::column_type;
using boost::mysql::common_server_errc;
using boost::mysql::diagnostics;
using boost::mysql::error_code;
using boost::mysql::execution_state;
using boost::mysql::detail::async_start_execution_generic;
using boost::mysql::detail::capabilities;
using boost::mysql::detail::db_flavor;
using boost::mysql::detail::deserialize_execute_response;
using boost::mysql::detail::execute_response;
using boost::mysql::detail::execution_state_access;
using boost::mysql::detail::protocol_field_type;
using boost::mysql::detail::resultset_encoding;
using boost::mysql::detail::start_execution_generic;
using namespace boost::mysql::test;

BOOST_TEST_DONT_PRINT_LOG_VALUE(execute_response::type_t)

namespace {

using serialize_fn = std::function<void(capabilities, std::vector<std::uint8_t>&)>;
using netfun_maker = netfun_maker_fn<
    void,
    resultset_encoding,
    test_channel&,
    const serialize_fn&,
    execution_state&>;

struct fns_t
{
    typename netfun_maker::signature start_execution_generic;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&start_execution_generic),           "sync" },
    {netfun_maker::async_errinfo(&async_start_execution_generic), "async"}
};

// Verify that we reset the state
execution_state create_initial_state()
{
    return create_execution_state(resultset_encoding::text, {protocol_field_type::geometry}, 4);
}

BOOST_AUTO_TEST_SUITE(test_start_execution_generic)

BOOST_AUTO_TEST_SUITE(deserialize_execute_response_)
BOOST_AUTO_TEST_CASE(ok_packet)
{
    std::uint8_t msg[] = {0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00};
    diagnostics diag;
    auto response = deserialize_execute_response(
        boost::asio::buffer(msg),
        capabilities(),
        db_flavor::mariadb,
        diag
    );
    BOOST_TEST(response.type == execute_response::type_t::ok_packet);
    BOOST_TEST(response.data.ok_pack.affected_rows.value == 0u);
    BOOST_TEST(response.data.ok_pack.status_flags == 2u);
}

BOOST_AUTO_TEST_CASE(num_fields)
{
    struct
    {
        const char* name;
        std::vector<std::uint8_t> msg;
        std::size_t num_fields;
    } test_cases[] = {
        {"1",                    {0x01},             1     },
        {"0xfa",                 {0xfa},             0xfa  },
        {"0xfb_no_local_infile", {0xfb},             0xfb  }, // legal when LOCAL INFILE capability not enabled
        {"0xfb_local_infile",    {0xfc, 0xfb, 0x00}, 0xfb  }, // sent LOCAL INFILE capability is enabled
        {"0xff",                 {0xfc, 0xff, 0x00}, 0xff  },
        {"0x01ff",               {0xfc, 0x00, 0x01}, 0x01ff},
        {"max",                  {0xfc, 0xff, 0xff}, 0xffff},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            std::uint8_t msg[] = {0xfc, 0xff, 0x00};
            diagnostics diag;
            auto response = deserialize_execute_response(
                boost::asio::buffer(msg),
                capabilities(),
                db_flavor::mysql,
                diag
            );
            BOOST_TEST(response.type == execute_response::type_t::num_fields);
            BOOST_TEST(response.data.num_fields == 0xffu);
            BOOST_TEST(diag.server_message() == "");
        }
    }
}

BOOST_AUTO_TEST_CASE(error)
{
    struct
    {
        const char* name;
        std::vector<std::uint8_t> msg;
        error_code err;
        const char* expected_info;
    } test_cases[] = {
        {"server_error",
         {0xff, 0x7a, 0x04, 0x23, 0x34, 0x32, 0x53, 0x30, 0x32, 0x54, 0x61, 0x62, 0x6c, 0x65,
          0x20, 0x27, 0x6d, 0x79, 0x74, 0x65, 0x73, 0x74, 0x2e, 0x61, 0x62, 0x63, 0x27, 0x20,
          0x64, 0x6f, 0x65, 0x73, 0x6e, 0x27, 0x74, 0x20, 0x65, 0x78, 0x69, 0x73, 0x74},
         common_server_errc::er_no_such_table,
         "Table 'mytest.abc' doesn't exist"                                                                                   },
        {"bad_server_error", {0xff, 0x00},                                               client_errc::incomplete_message,   ""},
        {"bad_ok_packet",    {0x00, 0xff},                                               client_errc::incomplete_message,   ""},
        {"bad_num_fields",   {0xfc, 0xff, 0x00, 0x01},                                   client_errc::extra_bytes,          ""},
        {"zero_num_fields",  {0xfc, 0x00, 0x00},                                         client_errc::protocol_value_error, ""},
        {"3byte_integer",    {0xfd, 0xff, 0xff, 0xff},                                   client_errc::protocol_value_error, ""},
        {"8byte_integer",
         {0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
         client_errc::protocol_value_error,
         ""                                                                                                                   },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            diagnostics diag;
            auto response = deserialize_execute_response(
                boost::asio::buffer(tc.msg),
                capabilities(),
                db_flavor::mysql,
                diag
            );
            BOOST_TEST(response.type == execute_response::type_t::error);
            BOOST_TEST(response.data.err == tc.err);
            BOOST_TEST(diag.server_message() == tc.expected_info);
        }
    }
}
BOOST_AUTO_TEST_SUITE_END()

struct fixture
{
    execution_state st{create_initial_state()};
    serialize_fn serializer{[](capabilities, std::vector<std::uint8_t>& buff) {
        std::uint8_t message[] = {0x01, 0x02, 0x03};
        buff.assign(std::begin(message), std::end(message));
    }};
    test_channel chan{create_channel()};

    fixture()
    {
        chan.set_meta_mode(boost::mysql::metadata_mode::full);
        chan.shared_sequence_number() = 42;
    }

    void check_written_message()
    {
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(
            chan.lowest_layer().bytes_written(),
            blob({0x03, 0x00, 0x00, 0x00, 0x01, 0x02, 0x03})
        );
        BOOST_TEST(chan.shared_sequence_number() == 42u);  // not used
    }

    network_result<void> call_fn(const fns_t& fns)
    {
        return fns.start_execution_generic(resultset_encoding::binary, chan, serializer, st);
    }
};

BOOST_AUTO_TEST_CASE(success_one_meta)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            auto response = create_message(1, {0x01});
            auto col = create_coldef_message(2, protocol_field_type::var_string);
            fix.chan.lowest_layer().add_message(concat_copy(response, col));

            // Call the function
            fix.call_fn(fns).validate_no_error();

            // We've written the request message
            fix.check_written_message();

            // We've read the response
            BOOST_TEST(execution_state_access::get_encoding(fix.st) == resultset_encoding::binary);
            BOOST_TEST(!fix.st.complete());
            BOOST_TEST(execution_state_access::get_sequence_number(fix.st) == 3u);
            BOOST_TEST_REQUIRE(fix.st.meta().size() == 1u);
            BOOST_TEST(fix.st.meta()[0].column_name() == "mycol");
            BOOST_TEST(fix.st.meta()[0].type() == column_type::varchar);
        }
    }
}

BOOST_AUTO_TEST_CASE(success_one_meta_metadata_minimal)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            auto response = create_message(1, {0x01});
            auto col = create_coldef_message(2, protocol_field_type::var_string);
            fix.chan.lowest_layer().add_message(concat_copy(response, col));
            fix.chan.set_meta_mode(boost::mysql::metadata_mode::minimal);

            // Call the function
            fix.call_fn(fns).validate_no_error();

            // We've written the request message
            fix.check_written_message();

            // We've read the response
            BOOST_TEST(execution_state_access::get_encoding(fix.st) == resultset_encoding::binary);
            BOOST_TEST(!fix.st.complete());
            BOOST_TEST(execution_state_access::get_sequence_number(fix.st) == 3u);
            BOOST_TEST_REQUIRE(fix.st.meta().size() == 1u);
            BOOST_TEST(fix.st.meta()[0].column_name() == "");
            BOOST_TEST(fix.st.meta()[0].type() == column_type::varchar);
        }
    }
}

BOOST_AUTO_TEST_CASE(success_several_meta_separate)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            auto response = create_message(1, {0x02});
            auto col1 = create_coldef_message(2, protocol_field_type::var_string, "f1");
            auto col2 = create_coldef_message(3, protocol_field_type::tiny, "f2");
            fix.chan.lowest_layer().add_message(response);
            fix.chan.lowest_layer().add_message(col1);
            fix.chan.lowest_layer().add_message(col2);

            // Call the function
            fix.call_fn(fns).validate_no_error();

            // We've written the request message
            fix.check_written_message();

            // We've read the response
            BOOST_TEST(execution_state_access::get_encoding(fix.st) == resultset_encoding::binary);
            BOOST_TEST(!fix.st.complete());
            BOOST_TEST(execution_state_access::get_sequence_number(fix.st) == 4u);
            BOOST_TEST_REQUIRE(fix.st.meta().size() == 2u);
            BOOST_TEST(fix.st.meta()[0].column_name() == "f1");
            BOOST_TEST(fix.st.meta()[0].type() == column_type::varchar);
            BOOST_TEST(fix.st.meta()[1].column_name() == "f2");
            BOOST_TEST(fix.st.meta()[1].type() == column_type::tinyint);
        }
    }
}

BOOST_AUTO_TEST_CASE(success_ok_packet)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            auto response = create_ok_packet_message(1, 42, 43, 44, 45, "abc");
            fix.chan.lowest_layer().add_message(response);

            // Call the function
            fix.call_fn(fns).validate_no_error();

            // We've written the request message
            fix.check_written_message();

            // We've read the response
            BOOST_TEST(execution_state_access::get_encoding(fix.st) == resultset_encoding::binary);
            BOOST_TEST(fix.st.meta().size() == 0u);
            BOOST_TEST_REQUIRE(fix.st.complete());
            BOOST_TEST(fix.st.affected_rows() == 42u);
            BOOST_TEST(fix.st.last_insert_id() == 43u);
            BOOST_TEST(fix.st.warning_count() == 45u);
            BOOST_TEST(fix.st.info() == "abc");
        }
    }
}

BOOST_AUTO_TEST_CASE(error_network_error)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            // This covers: error writing the request, error reading
            // the initial response, error reading successive metadata packets
            for (std::size_t i = 0; i <= 2; ++i)
            {
                BOOST_TEST_CONTEXT(i)
                {
                    fixture fix;
                    auto response = create_message(1, {0x02});
                    auto col1 = create_coldef_message(2, protocol_field_type::var_string, "f1");
                    auto col2 = create_coldef_message(3, protocol_field_type::tiny, "f2");
                    fix.chan.lowest_layer().add_message(response);
                    fix.chan.lowest_layer().add_message(col1);
                    fix.chan.lowest_layer().add_message(col2);
                    fix.chan.lowest_layer().set_fail_count(fail_count(i, client_errc::server_unsupported));

                    // Call the function
                    fix.call_fn(fns).validate_error_exact(client_errc::server_unsupported);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(error_metadata_packets_seqnum_mismatch)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            auto response = create_message(1, {0x02});
            auto col1 = create_coldef_message(2, protocol_field_type::var_string, "f1");
            auto col2 = create_coldef_message(4, protocol_field_type::tiny, "f2");
            fix.chan.lowest_layer().add_message(response);
            fix.chan.lowest_layer().add_message(col1);
            fix.chan.lowest_layer().add_message(col2);

            // Call the function
            fix.call_fn(fns).validate_error_exact(client_errc::sequence_number_mismatch);
        }
    }
}

// All cases where the deserialization of the execution_response
// yields an error are handled uniformly, so it's enough with this test
BOOST_AUTO_TEST_CASE(error_deserialize_execution_response)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            auto response = create_err_packet_message(1, common_server_errc::er_bad_db_error, "no_db");
            fix.chan.lowest_layer().add_message(response);

            // Call the function
            fix.call_fn(fns).validate_error_exact(common_server_errc::er_bad_db_error, "no_db");
        }
    }
}

BOOST_AUTO_TEST_CASE(error_deserialize_metadata)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            auto response = create_message(1, {0x01});
            auto col = create_message(2, {0x08, 0x03});
            fix.chan.lowest_layer().add_message(response);
            fix.chan.lowest_layer().add_message(col);

            // Call the function
            fix.call_fn(fns).validate_error_exact(client_errc::incomplete_message);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace