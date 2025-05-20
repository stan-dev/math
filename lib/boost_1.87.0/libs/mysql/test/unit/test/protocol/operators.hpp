//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_TEST_PROTOCOL_OPERATORS_HPP
#define BOOST_MYSQL_TEST_UNIT_TEST_PROTOCOL_OPERATORS_HPP

#include <boost/mysql/impl/internal/protocol/impl/deserialization_context.hpp>
#include <boost/mysql/impl/internal/protocol/impl/protocol_types.hpp>

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

// int_holder
template <class T>
inline bool operator==(int_holder<T> lhs, int_holder<T> rhs) noexcept
{
    return lhs.value == rhs.value;
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, int_holder<T> value)
{
    return os << value.value;
}

// int3
inline bool operator==(int3 lhs, int3 rhs) noexcept { return lhs.value == rhs.value; }
inline std::ostream& operator<<(std::ostream& os, int3 value) { return os << value.value; }

// int_lenenc
inline bool operator==(int_lenenc lhs, int_lenenc rhs) noexcept { return lhs.value == rhs.value; }
inline std::ostream& operator<<(std::ostream& os, int_lenenc value) { return os << value.value; }

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
