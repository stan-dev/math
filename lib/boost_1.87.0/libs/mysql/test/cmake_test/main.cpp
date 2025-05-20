//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/tcp.hpp>

#include <boost/asio/system_executor.hpp>

int main()
{
    boost::mysql::tcp_connection conn(boost::asio::system_executor{});
    return static_cast<int>(conn.uses_ssl());  // should be false for a non-connected connection
}
