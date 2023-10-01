//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_PRINTING_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_PRINTING_HPP

#include <boost/mysql/detail/results_iterator.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>
#include <boost/mysql/detail/rows_iterator.hpp>

#include <boost/mysql/impl/internal/protocol/capabilities.hpp>
#include <boost/mysql/impl/internal/protocol/db_flavor.hpp>

#include <boost/test/unit_test.hpp>

#include <ostream>

namespace boost {
namespace mysql {
namespace detail {

inline std::ostream& operator<<(std::ostream& os, capabilities caps)
{
    return os << "capabilities(" << caps.get() << ")";
}

inline std::ostream& operator<<(std::ostream& os, db_flavor value)
{
    switch (value)
    {
    case db_flavor::mysql: return os << "mysql";
    case db_flavor::mariadb: return os << "mariadb";
    default: return os << "<unknown db_flavor>";
    }
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

inline std::ostream& operator<<(std::ostream& os, results_iterator it)
{
    return os << "results_iterator(" << static_cast<const void*>(it.obj()) << ", index=" << it.index() << ")";
}

}  // namespace detail
}  // namespace mysql
}  // namespace boost

#endif
