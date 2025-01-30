//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/ssl/context.hpp>
#include <boost/asio/system_executor.hpp>

int main()
{
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    boost::mysql::tcp_ssl_connection conn(boost::asio::system_executor{}, ssl_ctx);
    return static_cast<int>(conn.uses_ssl());  // should be false for a non-connected connection
}
