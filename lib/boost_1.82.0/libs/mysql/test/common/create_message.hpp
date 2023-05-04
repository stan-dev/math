//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_CREATE_MESSAGE_HPP
#define BOOST_MYSQL_TEST_COMMON_CREATE_MESSAGE_HPP

#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/mysql_collations.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/protocol/capabilities.hpp>
#include <boost/mysql/detail/protocol/common_messages.hpp>
#include <boost/mysql/detail/protocol/protocol_types.hpp>
#include <boost/mysql/detail/protocol/serialization.hpp>
#include <boost/mysql/detail/protocol/serialization_context.hpp>

#include <cassert>
#include <cstdint>
#include <cstring>

#include "buffer_concat.hpp"

namespace boost {
namespace mysql {
namespace test {

inline std::vector<std::uint8_t> create_message(std::uint8_t seqnum, std::vector<std::uint8_t> body)
{
    assert(body.size() <= std::numeric_limits<std::uint32_t>::max());
    auto body_size = static_cast<std::uint32_t>(body.size());
    boost::mysql::detail::packet_header header{boost::mysql::detail::int3{body_size}, seqnum};
    body.resize(body_size + 4);
    std::memmove(body.data() + 4, body.data(), body_size);
    boost::mysql::detail::serialization_context ctx(boost::mysql::detail::capabilities(), body.data());
    boost::mysql::detail::serialize(ctx, header);
    return body;
}

inline std::vector<std::uint8_t> create_message(
    std::uint8_t seqnum1,
    std::vector<std::uint8_t> body1,
    std::uint8_t seqnum2,
    std::vector<std::uint8_t> body2
)
{
    return concat_copy(create_message(seqnum1, std::move(body1)), create_message(seqnum2, std::move(body2)));
}

inline std::vector<std::uint8_t> create_message(
    std::uint8_t seqnum1,
    std::vector<std::uint8_t> body1,
    std::uint8_t seqnum2,
    std::vector<std::uint8_t> body2,
    std::uint8_t seqnum3,
    std::vector<std::uint8_t> body3
)
{
    return concat_copy(
        create_message(seqnum1, std::move(body1)),
        create_message(seqnum2, std::move(body2)),
        create_message(seqnum3, std::move(body3))
    );
}

template <class... Args>
std::vector<std::uint8_t> serialize_to_vector(const Args&... args)
{
    std::vector<std::uint8_t> res;
    detail::serialization_context ctx(detail::capabilities(0));
    std::size_t size = detail::get_size(ctx, args...);
    res.resize(size);
    ctx.set_first(res.data());
    detail::serialize(ctx, args...);
    return res;
}

inline detail::ok_packet create_ok_packet(
    std::uint64_t affected_rows = 0,
    std::uint64_t last_insert_id = 0,
    std::uint16_t status_flags = 0,
    std::uint16_t warnings = 0,
    string_view info = ""
)
{
    return detail::ok_packet{
        detail::int_lenenc{affected_rows},
        detail::int_lenenc{last_insert_id},
        status_flags,
        warnings,
        detail::string_lenenc{info},
    };
}

inline std::vector<std::uint8_t> create_ok_packet_body(
    std::uint64_t affected_rows = 0,
    std::uint64_t last_insert_id = 0,
    std::uint16_t status_flags = 0,
    std::uint16_t warnings = 0,
    string_view info = "",
    std::uint8_t header = 0x00
)
{
    auto pack = create_ok_packet(affected_rows, last_insert_id, status_flags, warnings, info);
    auto res = serialize_to_vector(
        std::uint8_t(header),
        pack.affected_rows,
        pack.last_insert_id,
        pack.status_flags,
        pack.warnings
    );
    // When info is empty, it's actually omitted in the ok_packet
    if (!info.empty())
    {
        auto vinfo = serialize_to_vector(pack.info);
        concat(res, vinfo);
    }
    return res;
}

inline std::vector<std::uint8_t> create_ok_packet_message(
    std::uint8_t seqnum,
    std::uint64_t affected_rows = 0,
    std::uint64_t last_insert_id = 0,
    std::uint16_t status_flags = 0,
    std::uint16_t warnings = 0,
    string_view info = "",
    std::uint8_t header = 0x00
)
{
    return create_message(
        seqnum,
        create_ok_packet_body(affected_rows, last_insert_id, status_flags, warnings, info, header)
    );
}

inline std::vector<std::uint8_t> create_eof_packet_message(
    std::uint8_t seqnum,
    std::uint64_t affected_rows = 0,
    std::uint64_t last_insert_id = 0,
    std::uint16_t status_flags = 0,
    std::uint16_t warnings = 0,
    string_view info = ""
)
{
    return create_ok_packet_message(
        seqnum,
        affected_rows,
        last_insert_id,
        status_flags,
        warnings,
        info,
        0xfe
    );
}

inline std::vector<std::uint8_t> create_err_packet_body(std::uint16_t code, string_view message = "")
{
    detail::err_packet pack{
        code,
        detail::string_fixed<1>{},
        detail::string_fixed<5>{},
        detail::string_eof{message},
    };
    return serialize_to_vector(
        std::uint8_t(0xff),
        pack.error_code,
        pack.sql_state_marker,
        pack.sql_state,
        pack.error_message
    );
}

inline std::vector<std::uint8_t> create_err_packet_body(common_server_errc code, string_view message = "")
{
    return create_err_packet_body(static_cast<std::uint16_t>(code), message);
}

inline std::vector<std::uint8_t> create_err_packet_message(
    std::uint8_t seqnum,
    common_server_errc code,
    string_view message = ""
)
{
    return create_message(seqnum, create_err_packet_body(code, message));
}

inline std::vector<std::uint8_t> create_coldef_message(
    std::uint8_t seqnum,
    detail::protocol_field_type type,
    string_view name = "mycol"
)
{
    boost::mysql::detail::column_definition_packet pack{
        detail::string_lenenc("def"),
        detail::string_lenenc("mydb"),
        detail::string_lenenc("mytable"),
        detail::string_lenenc("mytable"),
        detail::string_lenenc(name),
        detail::string_lenenc(name),
        mysql_collations::utf8_general_ci,
        10,  // column_length
        type,
        0,  // flags
        0,  // decimals
    };
    return create_message(
        seqnum,
        serialize_to_vector(
            pack.catalog,
            pack.schema,
            pack.table,
            pack.org_table,
            pack.name,
            pack.org_name,
            boost::mysql::detail::int_lenenc(0x0c),  // length of fixed fields
            pack.character_set,
            pack.column_length,
            pack.type,
            pack.flags,
            pack.decimals,
            std::uint16_t(0)  // padding
        )
    );
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
