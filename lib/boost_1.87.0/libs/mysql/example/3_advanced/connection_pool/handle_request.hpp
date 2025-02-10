//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_EXAMPLE_3_ADVANCED_CONNECTION_POOL_HANDLE_REQUEST_HPP
#define BOOST_MYSQL_EXAMPLE_3_ADVANCED_CONNECTION_POOL_HANDLE_REQUEST_HPP

//[example_connection_pool_handle_request_hpp
//
// File: handle_request.hpp
//

#include <boost/asio/error.hpp>
#include <boost/asio/spawn.hpp>
#include <boost/beast/http/message.hpp>
#include <boost/beast/http/string_body.hpp>

#include "repository.hpp"

namespace notes {

// Handles an individual HTTP request, producing a response.
boost::beast::http::response<boost::beast::http::string_body> handle_request(
    const boost::beast::http::request<boost::beast::http::string_body>& request,
    note_repository repo,
    boost::asio::yield_context yield
);

}  // namespace notes

//]

#endif
