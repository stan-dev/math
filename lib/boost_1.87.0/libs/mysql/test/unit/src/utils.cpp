//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_address.hpp>
#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/character_set.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/coldef_view.hpp>
#include <boost/mysql/detail/engine_impl.hpp>
#include <boost/mysql/detail/engine_stream_adaptor.hpp>
#include <boost/mysql/detail/next_action.hpp>
#include <boost/mysql/detail/pipeline.hpp>
#include <boost/mysql/detail/results_iterator.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>

#include <boost/mysql/impl/internal/connection_pool/sansio_connection_node.hpp>
#include <boost/mysql/impl/internal/protocol/capabilities.hpp>
#include <boost/mysql/impl/internal/protocol/db_flavor.hpp>
#include <boost/mysql/impl/internal/protocol/deserialization.hpp>
#include <boost/mysql/impl/internal/protocol/frame_header.hpp>
#include <boost/mysql/impl/internal/protocol/impl/protocol_field_type.hpp>
#include <boost/mysql/impl/internal/protocol/impl/protocol_types.hpp>
#include <boost/mysql/impl/internal/protocol/impl/serialization_context.hpp>

#include <boost/asio/buffer.hpp>
#include <boost/asio/compose.hpp>
#include <boost/asio/error.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/post.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <utility>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/printing.hpp"
#include "test_unit/algo_test.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_prepare_statement_response.hpp"
#include "test_unit/create_query_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/mock_timer.hpp"
#include "test_unit/printing.hpp"
#include "test_unit/serialize_to_vector.hpp"
#include "test_unit/test_any_connection.hpp"
#include "test_unit/test_stream.hpp"

using std::chrono::steady_clock;
using namespace boost::mysql;
using namespace boost::mysql::test;

//
// algo_test.hpp
//
void boost::mysql::test::algo_test::handle_read(detail::connection_state_data& st, const step_t& op)
{
    if (!op.result)
    {
        std::size_t bytes_transferred = 0;
        while (!st.reader.done() && bytes_transferred < op.bytes.size())
        {
            auto ec = st.reader.prepare_buffer();
            BOOST_TEST_REQUIRE(ec == error_code());
            auto buff = st.reader.buffer();
            std::size_t size_to_copy = (std::min)(op.bytes.size() - bytes_transferred, buff.size());
            std::memcpy(buff.data(), op.bytes.data() + bytes_transferred, size_to_copy);
            bytes_transferred += size_to_copy;
            st.reader.resume(size_to_copy);
        }
        BOOST_TEST_REQUIRE(st.reader.done());
        BOOST_TEST_REQUIRE(st.reader.error() == error_code());
    }
}

detail::next_action boost::mysql::test::algo_test::run_algo_until_step(
    detail::connection_state_data& st,
    any_algo_ref algo,
    std::size_t num_steps_to_run
) const
{
    BOOST_ASSERT(num_steps_to_run <= num_steps());

    // Start the op
    auto act = algo.resume(st, error_code());

    // Go through the requested steps
    for (std::size_t i = 0; i < num_steps_to_run; ++i)
    {
        BOOST_TEST_CONTEXT("Step " << i)
        {
            const auto& step = steps_[i];
            BOOST_TEST_REQUIRE(act.type() == step.type);
            if (step.type == detail::next_action_type::read)
                handle_read(st, step);
            else if (step.type == detail::next_action_type::write)
                BOOST_MYSQL_ASSERT_BUFFER_EQUALS(act.write_args().buffer, step.bytes);
            // Other actions don't need any handling

            act = algo.resume(st, step.result);
        }
    }

    return act;
}

void boost::mysql::test::algo_test::check_network_errors_impl(
    detail::connection_state_data& st,
    any_algo_ref algo,
    std::size_t step_number,
    const diagnostics& actual_diag,
    source_location loc
) const
{
    BOOST_TEST_CONTEXT("Called from " << loc << " at step " << step_number)
    {
        BOOST_ASSERT(step_number < num_steps());

        // Run all the steps that shouldn't cause an error
        auto act = run_algo_until_step(st, algo, step_number);
        BOOST_TEST_REQUIRE(act.type() == steps_[step_number].type);

        // Trigger an error in the requested step
        act = algo.resume(st, asio::error::bad_descriptor);

        // The operation finished and returned the network error
        BOOST_TEST_REQUIRE(act.type() == detail::next_action_type::none);
        BOOST_TEST(act.error() == error_code(asio::error::bad_descriptor));
        BOOST_TEST(actual_diag == diagnostics());
    }
}

boost::mysql::test::algo_test& boost::mysql::test::algo_test::add_step(
    detail::next_action_type act_type,
    std::vector<std::uint8_t> bytes,
    error_code ec
)
{
    steps_.push_back(step_t{act_type, std::move(bytes), ec});
    return *this;
}

void boost::mysql::test::algo_test::check_impl(
    detail::connection_state_data& st,
    any_algo_ref algo,
    const diagnostics& actual_diag,
    error_code expected_ec,
    const diagnostics& expected_diag,
    source_location loc
) const
{
    BOOST_TEST_CONTEXT("Called from " << loc)
    {
        // Run the op until completion
        auto act = run_algo_until_step(st, algo, steps_.size());

        // Check that we've finished
        BOOST_TEST_REQUIRE(act.type() == detail::next_action_type::none);

        // Check results
        BOOST_TEST(act.error() == expected_ec);
        BOOST_TEST(actual_diag == expected_diag);
    }
}

//
// create_frame.hpp
//
std::vector<std::uint8_t> boost::mysql::test::create_frame(std::uint8_t seqnum, span<const std::uint8_t> body)
{
    BOOST_ASSERT(body.size() <= 0xffffff);  // it should fit in a single frame

    // Compose the frame header
    std::array<std::uint8_t, detail::frame_header_size> frame_header{};
    detail::serialize_frame_header(
        frame_header,
        detail::frame_header{static_cast<std::uint32_t>(body.size()), seqnum}
    );

    // Compose the frame.
    // Inserting the header separately (instead of using the range constructor)
    // avoids spurious gcc warnings
    std::vector<std::uint8_t> res;
    res.insert(res.end(), frame_header.begin(), frame_header.end());
    res.insert(res.end(), body.begin(), body.end());
    return res;
}

//
// create_ok_frame.hpp
//
std::vector<std::uint8_t> boost::mysql::test::serialize_ok_impl(
    const detail::ok_view& pack,
    std::uint8_t header
)
{
    return serialize_to_vector([=](detail::serialization_context& ctx) {
        ctx.serialize(
            detail::int1{header},
            detail::int_lenenc{pack.affected_rows},
            detail::int_lenenc{pack.last_insert_id},
            detail::int2{pack.status_flags},
            detail::int2{pack.warnings}
        );
        // When info is empty, it's actually omitted in the ok_packet
        if (!pack.info.empty())
        {
            detail::string_lenenc{pack.info}.serialize(ctx);
        }
    });
}

//
// create_coldef_frame.hpp
//
std::vector<std::uint8_t> boost::mysql::test::create_coldef_body(const detail::coldef_view& pack)
{
    auto to_protocol_type = [](column_type t) {
        // Note: we perform an approximate mapping, good enough for unit tests.
        // The actual mapping is not one to one and depends on flags
        using detail::protocol_field_type;
        switch (t)
        {
        case column_type::tinyint: return protocol_field_type::tiny;
        case column_type::smallint: return protocol_field_type::short_;
        case column_type::mediumint: return protocol_field_type::int24;
        case column_type::int_: return protocol_field_type::long_;
        case column_type::bigint: return protocol_field_type::longlong;
        case column_type::float_: return protocol_field_type::float_;
        case column_type::double_: return protocol_field_type::double_;
        case column_type::decimal: return protocol_field_type::newdecimal;
        case column_type::bit: return protocol_field_type::bit;
        case column_type::year: return protocol_field_type::year;
        case column_type::time: return protocol_field_type::time;
        case column_type::date: return protocol_field_type::date;
        case column_type::datetime: return protocol_field_type::datetime;
        case column_type::timestamp: return protocol_field_type::timestamp;
        case column_type::char_: return protocol_field_type::string;
        case column_type::varchar: return protocol_field_type::var_string;
        case column_type::binary: return protocol_field_type::string;
        case column_type::varbinary: return protocol_field_type::var_string;
        case column_type::text: return protocol_field_type::blob;
        case column_type::blob: return protocol_field_type::blob;
        case column_type::enum_: return protocol_field_type::enum_;
        case column_type::set: return protocol_field_type::set;
        case column_type::json: return protocol_field_type::json;
        case column_type::geometry: return protocol_field_type::geometry;
        default: BOOST_ASSERT(false); return protocol_field_type::var_string;  // LCOV_EXCL_LINE
        }
    };

    return serialize_to_vector([=](detail::serialization_context& ctx) {
        ctx.serialize(
            detail::string_lenenc{"def"},
            detail::string_lenenc{pack.database},
            detail::string_lenenc{pack.table},
            detail::string_lenenc{pack.org_table},
            detail::string_lenenc{pack.name},
            detail::string_lenenc{pack.org_name},
            detail::int_lenenc{0x0c},  // length of fixed fields
            detail::int2{pack.collation_id},
            detail::int4{pack.column_length},
            detail::int1{static_cast<std::uint8_t>(to_protocol_type(pack.type))},
            detail::int2{pack.flags},
            detail::int1{pack.decimals},
            detail::int2{0}  // padding
        );
    });
}

//
// create_err.hpp
//
std::vector<std::uint8_t> boost::mysql::test::serialize_err_impl(detail::err_view pack, bool with_header)
{
    return serialize_to_vector([=](detail::serialization_context& ctx) {
        if (with_header)
            ctx.add(0xff);  // header
        ctx.serialize(
            detail::int2{pack.error_code},
            detail::string_fixed<1>{},  // SQL state marker
            detail::string_fixed<5>{},  // SQL state
            detail::string_eof{pack.error_message}
        );
    });
}

//
// create_prepare_statement_response.hpp
//
std::vector<std::uint8_t> boost::mysql::test::prepare_stmt_response_builder::build() const
{
    auto body = serialize_to_vector([this](detail::serialization_context& ctx) {
        ctx.serialize(
            detail::int1{0u},             // OK header
            detail::int4{statement_id_},  // statement_id
            detail::int2{num_columns_},   // num columns
            detail::int2{num_params_},    // num_params
            detail::int1{0u},             // reserved
            detail::int2{90u}             // warning_count
        );
    });
    return create_frame(seqnum_, body);
}

//
// create_query_frame
//
std::vector<std::uint8_t> boost::mysql::test::create_query_body_impl(std::uint8_t command_id, string_view sql)
{
    return serialize_to_vector([=](detail::serialization_context& ctx) {
        ctx.add(command_id);
        ctx.add(detail::to_span(sql));
    });
}

//
// create_row_message.hpp
//
std::vector<std::uint8_t> boost::mysql::test::serialize_text_row_impl(span<const field_view> fields)
{
    return serialize_to_vector([=](detail::serialization_context& ctx) {
        for (field_view f : fields)
        {
            std::string s;
            switch (f.kind())
            {
            case field_kind::int64: s = std::to_string(f.get_int64()); break;
            case field_kind::uint64: s = std::to_string(f.get_uint64()); break;
            case field_kind::float_: s = std::to_string(f.get_float()); break;
            case field_kind::double_: s = std::to_string(f.get_double()); break;
            case field_kind::string: s = f.get_string(); break;
            case field_kind::blob: s.assign(f.get_blob().begin(), f.get_blob().end()); break;
            case field_kind::null: ctx.add(std::uint8_t(0xfb)); continue;
            default: throw std::runtime_error("create_text_row_message: type not implemented");
            }
            detail::string_lenenc{s}.serialize(ctx);
        }
    });
}

//
// mock_timer.hpp
//
// The current time
static thread_local steady_clock::time_point g_mock_now{};

boost::mysql::test::mock_clock::time_point boost::mysql::test::mock_clock::now() { return g_mock_now; }

void boost::mysql::test::mock_clock::advance_time_by(steady_clock::duration dur) { g_mock_now += dur; }

//
// test_stream.hpp
//

std::size_t boost::mysql::test::test_stream::get_size_to_read(std::size_t buffer_size) const
{
    auto it = read_break_offsets_.upper_bound(num_bytes_read_);
    std::size_t max_bytes_by_break = it == read_break_offsets_.end() ? std::size_t(-1)
                                                                     : *it - num_bytes_read_;
    return (std::min)({num_unread_bytes(), buffer_size, max_bytes_by_break});
}

std::size_t boost::mysql::test::test_stream::do_read(asio::mutable_buffer buff, error_code& ec)
{
    // Fail count
    error_code err = fail_count_.maybe_fail();
    if (err)
    {
        ec = err;
        return 0;
    }

    // If the user requested some bytes but we don't have any,
    // fail. In the real world, the stream would block until more
    // bytes are received, but this is a test, and this condition
    // indicates an error.
    if (num_unread_bytes() == 0 && buff.size() != 0)
    {
        ec = boost::asio::error::eof;
        return 0;
    }

    // Actually read
    std::size_t bytes_to_transfer = get_size_to_read(buff.size());
    if (bytes_to_transfer)
    {
        std::memcpy(buff.data(), bytes_to_read_.data() + num_bytes_read_, bytes_to_transfer);
        num_bytes_read_ += bytes_to_transfer;
    }

    // Clear errors
    ec = error_code();

    return bytes_to_transfer;
}

std::size_t boost::mysql::test::test_stream::do_write(asio::const_buffer buff, error_code& ec)
{
    // Fail count
    error_code err = fail_count_.maybe_fail();
    if (err)
    {
        ec = err;
        return 0;
    }

    // Actually write
    std::size_t num_bytes_to_transfer = (std::min)(buff.size(), write_break_size_);
    span<const std::uint8_t> span_to_transfer(
        static_cast<const std::uint8_t*>(buff.data()),
        num_bytes_to_transfer
    );
    bytes_written_.insert(bytes_written_.end(), span_to_transfer.begin(), span_to_transfer.end());

    // Clear errors
    ec = error_code();

    return num_bytes_to_transfer;
}

struct boost::mysql::test::test_stream::read_op
{
    test_stream& stream_;
    asio::mutable_buffer buff_;
    bool has_posted_{};

    read_op(test_stream& stream, asio::mutable_buffer buff) noexcept : stream_(stream), buff_(buff){};

    template <class Self>
    void operator()(Self& self)
    {
        if (!has_posted_)
        {
            // Post
            has_posted_ = true;
            asio::post(stream_.get_executor(), std::move(self));
        }
        else
        {
            // Complete
            error_code err;
            std::size_t bytes_read = stream_.do_read(buff_, err);
            self.complete(err, bytes_read);
        }
    }
};

struct boost::mysql::test::test_stream::write_op
{
    test_stream& stream_;
    asio::const_buffer buff_;
    bool has_posted_{};

    write_op(test_stream& stream, asio::const_buffer buff) noexcept : stream_(stream), buff_(buff){};

    template <class Self>
    void operator()(Self& self)
    {
        // Post
        if (!has_posted_)
        {
            has_posted_ = true;
            asio::post(stream_.get_executor(), std::move(self));
        }
        else
        {
            error_code err;
            std::size_t bytes_written = stream_.do_write(buff_, err);
            self.complete(err, bytes_written);
        }
    }
};

std::size_t boost::mysql::test::test_stream::read_some(asio::mutable_buffer buff, error_code& ec)
{
    return do_read(buff, ec);
}
void boost::mysql::test::test_stream::async_read_some(
    asio::mutable_buffer buff,
    asio::any_completion_handler<void(error_code, std::size_t)> handler
)
{
    asio::async_compose<
        asio::any_completion_handler<void(error_code, std::size_t)>,
        void(error_code, std::size_t)>(read_op(*this, buff), handler, get_executor());
}

std::size_t boost::mysql::test::test_stream::write_some(boost::asio::const_buffer buff, error_code& ec)
{
    return do_write(buff, ec);
}

void boost::mysql::test::test_stream::async_write_some(
    asio::const_buffer buff,
    asio::any_completion_handler<void(error_code, std::size_t)> handler
)
{
    asio::async_compose<
        asio::any_completion_handler<void(error_code, std::size_t)>,
        void(error_code, std::size_t)>(write_op(*this, buff), handler, get_executor());
}

test_stream& boost::mysql::test::test_stream::add_bytes(span<const std::uint8_t> bytes)
{
    bytes_to_read_.insert(bytes_to_read_.end(), bytes.begin(), bytes.end());
    return *this;
}

test_stream& boost::mysql::test::test_stream::add_break(std::size_t byte_num)
{
    BOOST_ASSERT(byte_num <= bytes_to_read_.size());
    read_break_offsets_.insert(byte_num);
    return *this;
}

//
// printing.hpp
//

// address_type
static const char* to_string(address_type v)
{
    switch (v)
    {
    case address_type::host_and_port: return "address_type::host_and_port";
    case address_type::unix_path: return "address_type::unix_path";
    default: return "<unknown address_type>";
    }
}

std::ostream& boost::mysql::operator<<(std::ostream& os, address_type v) { return os << ::to_string(v); }

// capabilities
std::ostream& boost::mysql::detail::operator<<(std::ostream& os, const capabilities& v)
{
    return os << "capabilities{" << v.get() << "}";
}

// db_flavor
static const char* to_string(detail::db_flavor v)
{
    switch (v)
    {
    case detail::db_flavor::mysql: return "db_flavor::mysql";
    case detail::db_flavor::mariadb: return "db_flavor::mariadb";
    default: return "<unknown db_flavor>";
    }
}

std::ostream& boost::mysql::detail::operator<<(std::ostream& os, db_flavor v) { return os << ::to_string(v); }

// resultset_encoding
static const char* to_string(detail::resultset_encoding v)
{
    switch (v)
    {
    case detail::resultset_encoding::text: return "resultset_encoding::text";
    case detail::resultset_encoding::binary: return "resultset_encoding::binary";
    default: return "<unknown resultset_encoding>";
    }
}

std::ostream& boost::mysql::detail::operator<<(std::ostream& os, detail::resultset_encoding v)
{
    return os << ::to_string(v);
}

// results_iterator
std::ostream& boost::mysql::detail::operator<<(std::ostream& os, const results_iterator& it)
{
    return os << "results_iterator{ .self = " << static_cast<const void*>(it.obj())
              << ", .index = " << it.index() << "}";
}

// next_action_type;
static const char* to_string(detail::next_action_type v)
{
    switch (v)
    {
    case detail::next_action_type::none: return "next_action_type::none";
    case detail::next_action_type::read: return "next_action_type::read";
    case detail::next_action_type::write: return "next_action_type::write";
    case detail::next_action_type::ssl_handshake: return "next_action_type::ssl_handshake";
    case detail::next_action_type::ssl_shutdown: return "next_action_type::ssh_shutdown";
    default: return "<unknown next_action_type>";
    }
}

std::ostream& boost::mysql::detail::operator<<(std::ostream& os, next_action_type v)
{
    return os << ::to_string(v);
}

// pipeline_stage_kind
static const char* to_string(detail::pipeline_stage_kind v)
{
    switch (v)
    {
    case detail::pipeline_stage_kind::execute: return "pipeline_stage_kind::execute";
    case detail::pipeline_stage_kind::prepare_statement: return "pipeline_stage_kind::prepare_statement";
    case detail::pipeline_stage_kind::close_statement: return "pipeline_stage_kind::close_statement";
    case detail::pipeline_stage_kind::reset_connection: return "pipeline_stage_kind::reset_connection";
    case detail::pipeline_stage_kind::set_character_set: return "pipeline_stage_kind::set_character_set";
    case detail::pipeline_stage_kind::ping: return "pipeline_stage_kind::ping";
    default: return "<unknown pipeline_stage_kind>";
    }
}

std::ostream& boost::mysql::detail::operator<<(std::ostream& os, pipeline_stage_kind v)
{
    return os << ::to_string(v);
}

// pipeline_request_stage
bool boost::mysql::detail::operator==(const pipeline_request_stage& lhs, const pipeline_request_stage& rhs)
{
    if (lhs.kind != rhs.kind || lhs.seqnum != rhs.seqnum)
        return false;
    switch (lhs.kind)
    {
    case pipeline_stage_kind::execute: return lhs.stage_specific.enc == rhs.stage_specific.enc;
    case pipeline_stage_kind::set_character_set:
        return lhs.stage_specific.charset == rhs.stage_specific.charset;
    default: return true;
    }
}

std::ostream& boost::mysql::detail::operator<<(std::ostream& os, pipeline_request_stage v)
{
    os << "pipeline_request_stage{ .kind = " << v.kind << ", .seqnum = " << +v.seqnum;
    switch (v.kind)
    {
    case pipeline_stage_kind::execute: os << ", .enc = " << v.stage_specific.enc; break;
    case pipeline_stage_kind::set_character_set: os << ", .charset = " << v.stage_specific.charset; break;
    default: break;
    }
    return os << " }";
}

// connection_status
static const char* to_string(detail::connection_status v)
{
    switch (v)
    {
    case detail::connection_status::initial: return "connection_status::initial";
    case detail::connection_status::connect_in_progress: return "connection_status::connect_in_progress";
    case detail::connection_status::sleep_connect_failed_in_progress:
        return "connection_status::sleep_connect_failed_in_progress";
    case detail::connection_status::reset_in_progress: return "connection_status::reset_in_progress";
    case detail::connection_status::ping_in_progress: return "connection_status::ping_in_progress";
    case detail::connection_status::idle: return "connection_status::idle";
    case detail::connection_status::in_use: return "connection_status::in_use";
    default: return "<unknown connection_status>";
    }
}

std::ostream& boost::mysql::detail::operator<<(std::ostream& os, connection_status v)
{
    return os << ::to_string(v);
}

// collection_state
static const char* to_string(detail::collection_state v)
{
    switch (v)
    {
    case detail::collection_state::needs_collect: return "collection_state::needs_collect";
    case detail::collection_state::needs_collect_with_reset:
        return "collection_state::needs_collect_with_reset";
    case detail::collection_state::none: return "collection_state::none";
    default: return "<unknown collection_state>";
    }
}

std::ostream& boost::mysql::detail::operator<<(std::ostream& os, collection_state v)
{
    return os << ::to_string(v);
}

static const char* to_string(detail::next_connection_action v)
{
    switch (v)
    {
    case detail::next_connection_action::none: return "next_connection_action::none";
    case detail::next_connection_action::connect: return "next_connection_action::connect";
    case detail::next_connection_action::sleep_connect_failed:
        return "next_connection_action::sleep_connect_failed";
    case detail::next_connection_action::idle_wait: return "next_connection_action::idle_wait";
    case detail::next_connection_action::reset: return "next_connection_action::reset";
    case detail::next_connection_action::ping: return "next_connection_action::ping";
    default: return "<unknown next_connection_action>";
    }
}

std::ostream& boost::mysql::detail::operator<<(std::ostream& os, next_connection_action v)
{
    return os << ::to_string(v);
}

//
// test_any_connection.hpp
//
boost::mysql::any_connection boost::mysql::test::create_test_any_connection(
    asio::io_context& ctx,
    any_connection_params params
)
{
    return any_connection(detail::access::construct<any_connection>(
        std::unique_ptr<detail::engine>(
            new detail::engine_impl<detail::engine_stream_adaptor<test_stream>>(ctx.get_executor())
        ),
        params
    ));
}

test_stream& boost::mysql::test::get_stream(any_connection& conn)
{
    return detail::stream_from_engine<test_stream>(detail::access::get_impl(conn).get_engine());
}