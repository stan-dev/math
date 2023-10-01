//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include "test_unit/serialization.hpp"

#include <boost/mysql/column_type.hpp>

#include <boost/mysql/impl/internal/protocol/basic_types.hpp>
#include <boost/mysql/impl/internal/protocol/protocol.hpp>
#include <boost/mysql/impl/internal/protocol/protocol_field_type.hpp>
#include <boost/mysql/impl/internal/protocol/serialization.hpp>

#include <stdexcept>

using namespace boost::mysql::detail;
using boost::mysql::column_type;

template <class... Args>
void serialize_to_vector_inplace(std::vector<std::uint8_t>& res, const Args&... args)
{
    std::size_t size = get_size(args...);
    std::size_t old_size = res.size();
    res.resize(old_size + size);
    serialization_context ctx(res.data() + old_size);
    serialize(ctx, args...);
}

template <class... Args>
static std::vector<std::uint8_t> serialize_to_vector(const Args&... args)
{
    std::vector<std::uint8_t> res;
    serialize_to_vector_inplace(res, args...);
    return res;
}

static std::vector<std::uint8_t> serialize_ok_impl(const ok_view& pack, std::uint8_t header)
{
    auto res = serialize_to_vector(
        std::uint8_t(header),
        int_lenenc{pack.affected_rows},
        int_lenenc{pack.last_insert_id},
        pack.status_flags,
        pack.warnings
    );
    // When info is empty, it's actually omitted in the ok_packet
    if (!pack.info.empty())
    {
        serialize_to_vector_inplace(res, string_lenenc{pack.info});
    }
    return res;
}

std::vector<std::uint8_t> boost::mysql::test::serialize_ok(const ok_view& pack)
{
    return serialize_ok_impl(pack, 0x00);
}

std::vector<std::uint8_t> boost::mysql::test::serialize_eof(const ok_view& pack)
{
    return serialize_ok_impl(pack, 0xfe);
}

// Trying to merge common code from the two functions below hits a gcc-13 codegen bug
std::vector<std::uint8_t> boost::mysql::test::serialize_err_without_header(const err_view& pack)
{
    return serialize_to_vector(
        pack.error_code,
        string_fixed<1>{},  // SQL state marker
        string_fixed<5>{},  // SQL state
        string_eof{pack.error_message}
    );
}

std::vector<std::uint8_t> boost::mysql::test::serialize_err(const err_view& pack)
{
    return serialize_to_vector(
        std::uint8_t(0xff),  // header
        pack.error_code,
        string_fixed<1>{},  // SQL state marker
        string_fixed<5>{},  // SQL state
        string_eof{pack.error_message}
    );
}

static protocol_field_type to_protocol_type(column_type t) noexcept
{
    // Note: we perform an approximate mapping, good enough for unit tests.
    // The actual mapping is not one to one and depends on flags
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
    default: BOOST_ASSERT(false); return protocol_field_type::var_string;
    }
}

std::vector<std::uint8_t> boost::mysql::test::serialize_coldef(const detail::coldef_view& pack)
{
    return serialize_to_vector(
        string_lenenc{"def"},
        string_lenenc{pack.database},
        string_lenenc{pack.table},
        string_lenenc{pack.org_table},
        string_lenenc{pack.name},
        string_lenenc{pack.org_name},
        int_lenenc{0x0c},  // length of fixed fields
        pack.collation_id,
        pack.column_length,
        to_protocol_type(pack.type),
        pack.flags,
        pack.decimals,
        std::uint16_t(0)  // padding
    );
}

std::vector<std::uint8_t> boost::mysql::test::serialize_text_row(span<const field_view> fields)
{
    std::vector<std::uint8_t> res;
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
        case field_kind::null: serialize_to_vector_inplace(res, std::uint8_t(0xfb)); continue;
        default: throw std::runtime_error("create_text_row_message: type not implemented");
        }
        detail::string_lenenc slenenc{s};
        serialize_to_vector_inplace(res, slenenc);
    }
    return res;
}