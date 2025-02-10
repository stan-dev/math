
/* Copyright (c) 2018-2024 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/redis/request.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/use_future.hpp>
#include <thread>
#include <chrono>

using namespace std::chrono_literals;

namespace boost::redis
{

class sync_connection {
public:
   sync_connection()
   : ioc_{1}
   , conn_{std::make_shared<connection>(ioc_)}
   { }

   ~sync_connection()
   {
      thread_.join();
   }

   void run(config cfg)
   {
      // Starts a thread that will can io_context::run on which the
      // connection will run.
      thread_ = std::thread{[this, cfg]() {
         conn_->async_run(cfg, {}, asio::detached);
         ioc_.run();
      }};
   }

   void stop()
   {
      asio::dispatch(ioc_, [this]() { conn_->cancel(); });
   }

   template <class Response>
   auto exec(request const& req, Response& resp)
   {
      asio::dispatch(
         conn_->get_executor(),
         asio::deferred([this, &req, &resp]() { return conn_->async_exec(req, resp, asio::deferred); }))
         (asio::use_future).get();
   }

private:
   asio::io_context ioc_{1};
   std::shared_ptr<connection> conn_;
   std::thread thread_;
};

}
