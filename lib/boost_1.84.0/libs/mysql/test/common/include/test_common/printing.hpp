//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_PRINTING_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_PRINTING_HPP

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>

#include <boost/test/unit_test.hpp>

#include <ostream>

namespace boost {
namespace mysql {

inline std::ostream& operator<<(std::ostream& os, client_errc v)
{
    return os << get_client_category().message(static_cast<int>(v));
}

inline std::ostream& operator<<(std::ostream& os, common_server_errc v)
{
    return os << get_common_server_category().message(static_cast<int>(v));
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

}  // namespace mysql
}  // namespace boost

BOOST_TEST_DONT_PRINT_LOG_VALUE(boost::mysql::time)

#endif
