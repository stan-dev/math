//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/buffer.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <functional>

#include "assert_buffer_equals.hpp"
#include "printing.hpp"
#include "serialization_test.hpp"
#include "serialization_test_samples/basic_types.hpp"
#include "serialization_test_samples/binary_serialization.hpp"
#include "serialization_test_samples/common_messages.hpp"
#include "serialization_test_samples/handshake_messages.hpp"
#include "serialization_test_samples/prepared_statement_messages.hpp"
#include "serialization_test_samples/query_messages.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(test_serialization)

struct samples_by_type
{
    std::vector<serialization_sample> serialization_samples;
    std::vector<serialization_sample> deserialization_samples;
    std::vector<serialization_sample> space_samples;
};

void add_samples(const serialization_test_spec& spec, std::vector<serialization_sample>& to)
{
    to.insert(to.end(), spec.samples.begin(), spec.samples.end());
}

samples_by_type make_all_samples()
{
    const serialization_test_spec* all_specs[]{
        &int_spec,
        &enum_spec,
        &string_fixed_spec,
        &string_null_spec,
        &string_lenenc_spec,
        &string_eof_spec,

        &packet_header_spec,
        &ok_packet_spec,
        &err_packet_spec,
        &column_definition_spec,
        &quit_packet_spec,
        &ping_packet_spec,

        &handshake_packet_spec,
        &handshake_response_packet_spec,
        &auth_switch_request_packet_spec,
        &auth_switch_response_packet_spec,
        &ssl_request_spec,
        &auth_more_data_packet_spec,

        &com_query_packet_spec,

        &binary_serialization_spec,

        &com_stmt_prepare_packet_spec,
        &com_stmt_prepare_ok_packet_spec,
        &com_stmt_execute_packet_spec,
        &com_stmt_close_packet_spec,
    };

    samples_by_type res;

    for (const auto* spec : all_specs)
    {
        switch (spec->type)
        {
        case serialization_test_type::serialization: add_samples(*spec, res.serialization_samples); break;
        case serialization_test_type::deserialization: add_samples(*spec, res.deserialization_samples); break;
        case serialization_test_type::deserialization_space:
            add_samples(*spec, res.deserialization_samples);
            add_samples(*spec, res.space_samples);
            break;
        case serialization_test_type::full_no_space:
            add_samples(*spec, res.serialization_samples);
            add_samples(*spec, res.deserialization_samples);
            break;
        case serialization_test_type::full:
            add_samples(*spec, res.serialization_samples);
            add_samples(*spec, res.deserialization_samples);
            add_samples(*spec, res.space_samples);
            break;
        default: assert(false);
        }
    }
    return res;
}

samples_by_type all_samples = make_all_samples();

BOOST_DATA_TEST_CASE(get_size, data::make(all_samples.serialization_samples))
{
    serialization_context ctx(sample.caps, nullptr);
    auto size = sample.value->get_size(ctx);
    BOOST_TEST(size == sample.expected_buffer.size());
}

BOOST_DATA_TEST_CASE(serialize, data::make(all_samples.serialization_samples))
{
    auto expected_size = sample.expected_buffer.size();
    std::vector<uint8_t> buffer(expected_size + 8, 0x7a);  // buffer overrun detector
    serialization_context ctx(sample.caps, buffer.data());
    sample.value->serialize(ctx);

    // Buffer
    boost::asio::const_buffer expected_populated(sample.expected_buffer.data(), expected_size);
    boost::asio::const_buffer actual_populated(buffer.data(), sample.value->get_size(ctx));
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(expected_populated, actual_populated);

    // Iterator
    BOOST_TEST(ctx.first() == buffer.data() + expected_size, "Iterator not updated correctly");

    // Check for buffer overruns
    std::vector<std::uint8_t> expected_clean(8, 0x7a);
    boost::asio::const_buffer actual_clean(buffer.data() + expected_size, 8);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(boost::asio::buffer(expected_clean), actual_clean);
}

BOOST_DATA_TEST_CASE(deserialize, data::make(all_samples.deserialization_samples))
{
    auto first = sample.expected_buffer.data();
    auto size = sample.expected_buffer.size();
    deserialization_context ctx(first, first + size, sample.caps);
    auto actual_value = sample.value->default_construct();
    auto err = actual_value->deserialize(ctx);

    // No error
    BOOST_TEST(err == deserialize_errc::ok);

    // Iterator advanced
    BOOST_TEST(ctx.first() == first + size);

    // Actual value
    BOOST_TEST(*actual_value == *sample.value);
}

BOOST_DATA_TEST_CASE(deserialize_extra_space, data::make(all_samples.space_samples))
{
    std::vector<uint8_t> buffer(sample.expected_buffer);
    buffer.push_back(0xff);
    auto first = buffer.data();
    deserialization_context ctx(first, first + buffer.size(), sample.caps);
    auto actual_value = sample.value->default_construct();
    auto err = actual_value->deserialize(ctx);

    // No error
    BOOST_TEST(err == deserialize_errc::ok);

    // Iterator advanced
    BOOST_TEST(ctx.first() == first + sample.expected_buffer.size());

    // Actual value
    BOOST_TEST(*actual_value == *sample.value);
}

BOOST_DATA_TEST_CASE(deserialize_not_enough_space, data::make(all_samples.space_samples))
{
    std::vector<uint8_t> buffer(sample.expected_buffer);
    buffer.back() = 0x7a;  // try to detect any overruns
    deserialization_context ctx(buffer.data(), buffer.data() + buffer.size() - 1, sample.caps);
    auto actual_value = sample.value->default_construct();
    auto err = actual_value->deserialize(ctx);
    BOOST_TEST(err == deserialize_errc::incomplete_message);
}

BOOST_AUTO_TEST_SUITE_END()
