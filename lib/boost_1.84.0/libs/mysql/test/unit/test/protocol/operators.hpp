//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_TEST_PROTOCOL_OPERATORS_HPP
#define BOOST_MYSQL_TEST_UNIT_TEST_PROTOCOL_OPERATORS_HPP

#include <boost/mysql/impl/internal/protocol/basic_types.hpp>
#include <boost/mysql/impl/internal/protocol/protocol_field_type.hpp>
#include <boost/mysql/impl/internal/protocol/serialization.hpp>

#include <cstring>
#include <ostream>

namespace boost {
namespace mysql {
namespace detail {

// deserialize_errc
inline std::ostream& operator<<(std::ostream& os, deserialize_errc v)
{
    switch (v)
    {
    case deserialize_errc::ok: return os << "ok";
    case deserialize_errc::incomplete_message: return os << "incomplete_message";
    case deserialize_errc::protocol_value_error: return os << "protocol_value_error";
    case deserialize_errc::server_unsupported: return os << "server_unsupported";
    default: return os << "unknown";
    }
}

// int3
inline bool operator==(int3 lhs, int3 rhs) noexcept { return lhs.value == rhs.value; }
inline std::ostream& operator<<(std::ostream& os, int3 value) { return os << value.value; }

// int_lenenc
inline bool operator==(int_lenenc lhs, int_lenenc rhs) noexcept { return lhs.value == rhs.value; }
inline std::ostream& operator<<(std::ostream& os, int_lenenc value) { return os << value.value; }

// protocol_field_type
inline std::ostream& operator<<(std::ostream& os, protocol_field_type t)
{
    switch (t)
    {
    case detail::protocol_field_type::decimal: return os << "decimal";
    case detail::protocol_field_type::tiny: return os << "tiny";
    case detail::protocol_field_type::short_: return os << "short_";
    case detail::protocol_field_type::long_: return os << "long_";
    case detail::protocol_field_type::float_: return os << "float_";
    case detail::protocol_field_type::double_: return os << "double_";
    case detail::protocol_field_type::null: return os << "null";
    case detail::protocol_field_type::timestamp: return os << "timestamp";
    case detail::protocol_field_type::longlong: return os << "longlong";
    case detail::protocol_field_type::int24: return os << "int24";
    case detail::protocol_field_type::date: return os << "date";
    case detail::protocol_field_type::time: return os << "time";
    case detail::protocol_field_type::datetime: return os << "datetime";
    case detail::protocol_field_type::year: return os << "year";
    case detail::protocol_field_type::varchar: return os << "varchar";
    case detail::protocol_field_type::bit: return os << "bit";
    case detail::protocol_field_type::newdecimal: return os << "newdecimal";
    case detail::protocol_field_type::enum_: return os << "enum_";
    case detail::protocol_field_type::set: return os << "set";
    case detail::protocol_field_type::tiny_blob: return os << "tiny_blob";
    case detail::protocol_field_type::medium_blob: return os << "medium_blob";
    case detail::protocol_field_type::long_blob: return os << "long_blob";
    case detail::protocol_field_type::blob: return os << "blob";
    case detail::protocol_field_type::var_string: return os << "var_string";
    case detail::protocol_field_type::string: return os << "string";
    case detail::protocol_field_type::geometry: return os << "geometry";
    default: return os << "unknown";
    }
}

// string_fixed
template <std::size_t N>
inline bool operator==(string_fixed<N> lhs, string_fixed<N> rhs) noexcept
{
    return std::memcmp(lhs.value.data(), rhs.value.data(), N) == 0;
}

template <std::size_t N>
inline std::ostream& operator<<(std::ostream& os, string_fixed<N> value)
{
    return os << string_view(value.value.data(), N);
}

// string_null
inline bool operator==(string_null lhs, string_null rhs) noexcept { return lhs.value == rhs.value; }
inline std::ostream& operator<<(std::ostream& os, string_null value) { return os << value.value; }

// string_lenenc
inline bool operator==(string_lenenc lhs, string_lenenc rhs) noexcept { return lhs.value == rhs.value; }
inline std::ostream& operator<<(std::ostream& os, string_lenenc value) { return os << value.value; }

// string_eof
inline bool operator==(string_eof lhs, string_eof rhs) noexcept { return lhs.value == rhs.value; }
inline std::ostream& operator<<(std::ostream& os, string_eof value) { return os << value.value; }

}  // namespace detail
}  // namespace mysql
}  // namespace boost

#endif
