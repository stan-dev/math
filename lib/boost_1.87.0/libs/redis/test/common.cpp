#include "common.hpp"
#include <iostream>
#include <cstdlib>
#include <boost/asio/consign.hpp>
#include <boost/asio/co_spawn.hpp>

namespace net = boost::asio;

struct run_callback {
   std::shared_ptr<boost::redis::connection> conn;
   boost::redis::operation op;
   boost::system::error_code expected;

   void operator()(boost::system::error_code const& ec) const
   {
      std::cout << "async_run: " << ec.message() << std::endl;
      conn->cancel(op);
   }
};

void
run(
   std::shared_ptr<boost::redis::connection> conn,
   boost::redis::config cfg,
   boost::system::error_code ec,
   boost::redis::operation op,
   boost::redis::logger::level l)
{
   conn->async_run(cfg, {l}, run_callback{conn, op, ec});
}

static std::string safe_getenv(const char* name, const char* default_value)
{
    // MSVC doesn't like getenv
#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
    const char* res = std::getenv(name);
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif
    return res ? res : default_value;
}

std::string get_server_hostname()
{
   return safe_getenv("BOOST_REDIS_TEST_SERVER", "localhost");
}

boost::redis::config make_test_config()
{
   boost::redis::config cfg;
   cfg.addr.host = get_server_hostname();
   return cfg;
}

#ifdef BOOST_ASIO_HAS_CO_AWAIT
auto start(net::awaitable<void> op) -> int
{
   try {
      net::io_context ioc;
      net::co_spawn(ioc, std::move(op), [](std::exception_ptr p) {
         if (p)
            std::rethrow_exception(p);
      });
      ioc.run();

      return 0;

   } catch (std::exception const& e) {
      std::cerr << "start> " << e.what() << std::endl;
   }

   return 1;
}
#endif // BOOST_ASIO_HAS_CO_AWAIT
