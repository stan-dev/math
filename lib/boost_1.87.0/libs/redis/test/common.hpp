#pragma once

#include <boost/system/error_code.hpp>
#include <boost/asio/redirect_error.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/redis/connection.hpp>
#include <boost/redis/operation.hpp>
#include <memory>

#ifdef BOOST_ASIO_HAS_CO_AWAIT
inline
auto redir(boost::system::error_code& ec)
   { return boost::asio::redirect_error(boost::asio::use_awaitable, ec); }
auto start(boost::asio::awaitable<void> op) -> int;
#endif // BOOST_ASIO_HAS_CO_AWAIT

boost::redis::config make_test_config();
std::string get_server_hostname();

void
run(
   std::shared_ptr<boost::redis::connection> conn,
   boost::redis::config cfg = make_test_config(),
   boost::system::error_code ec = boost::asio::error::operation_aborted,
   boost::redis::operation op = boost::redis::operation::receive,
   boost::redis::logger::level l = boost::redis::logger::level::disabled);

