//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_PRINTING_HPP
#define BOOST_MYSQL_TEST_COMMON_PRINTING_HPP

#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>

#include <boost/mysql/detail/auxiliar/rows_iterator.hpp>
#include <boost/mysql/detail/auxiliar/static_string.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>
#include <boost/mysql/detail/protocol/deserialize_errc.hpp>
#include <boost/mysql/detail/protocol/resultset_encoding.hpp>

#include <boost/test/unit_test.hpp>

#include <ostream>

namespace boost {
namespace mysql {

inline std::ostream& operator<<(std::ostream& os, client_errc v) { return os << detail::error_to_string(v); }

inline std::ostream& operator<<(std::ostream& os, common_server_errc v)
{
    return os << detail::error_to_string(v);
}

inline std::ostream& operator<<(std::ostream& os, const diagnostics& diag)
{
    return os << diag.server_message();
}

inline std::ostream& operator<<(std::ostream& os, const row_view& value)
{
    os << '{';
    if (!value.empty())
    {
        os << value[0];
        for (auto it = std::next(value.begin()); it != value.end(); ++it)
        {
            os << ", " << *it;
        }
    }
    return os << '}';
}

inline std::ostream& operator<<(std::ostream& os, const row& r) { return os << row_view(r); }

inline std::ostream& operator<<(std::ostream& os, metadata_mode v)
{
    switch (v)
    {
    case metadata_mode::full: return os << "full";
    case metadata_mode::minimal: return os << "minimal";
    default: return os << "<unknown metadata_mode>";
    }
}

namespace detail {

template <std::size_t max_size>
std::ostream& operator<<(std::ostream& os, const static_string<max_size>& value)
{
    return os << value.value();
}

inline std::ostream& operator<<(std::ostream& os, resultset_encoding t)
{
    switch (t)
    {
    case resultset_encoding::binary: return os << "binary";
    case resultset_encoding::text: return os << "text";
    default: return os << "unknown";
    }
}

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

}  // namespace detail
}  // namespace mysql
}  // namespace boost

BOOST_TEST_DONT_PRINT_LOG_VALUE(boost::mysql::time)
BOOST_TEST_DONT_PRINT_LOG_VALUE(boost::mysql::detail::rows_iterator)

#endif
