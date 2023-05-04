//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_CREATE_STATEMENT_HPP
#define BOOST_MYSQL_TEST_COMMON_CREATE_STATEMENT_HPP

#include <boost/mysql/statement.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>

#include <cstdint>

namespace boost {
namespace mysql {
namespace test {

inline statement create_statement(std::uint16_t num_params, std::uint32_t stmt_id = 1)
{
    statement stmt;
    detail::statement_access::reset(
        stmt,
        boost::mysql::detail::com_stmt_prepare_ok_packet{stmt_id, 2, num_params, 0}
    );
    return stmt;
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
