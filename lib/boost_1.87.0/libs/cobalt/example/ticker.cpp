 // Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt.hpp>
#include <boost/cobalt/main.hpp>
#include <boost/cobalt/join.hpp>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl.hpp>
#include <boost/asio/as_tuple.hpp>
#include <boost/beast/http.hpp>
#include <boost/beast/websocket.hpp>
#include <boost/beast/core/flat_buffer.hpp>
#include <boost/beast/websocket/ssl.hpp>
#include <boost/json.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/unordered_map.hpp>
// cause cmake
#include <boost/url.hpp>

#include <list>
#include <iostream>

// https://blockchain.info/ticker
// https://api.coingecko.com/api/v3/coins/list

using namespace boost;

// tag::decls[]
using executor_type = cobalt::use_op_t::executor_with_default<cobalt::executor>;
using socket_type   = typename asio::ip::tcp::socket::rebind_executor<executor_type>::other;
using acceptor_type = typename asio::ip::tcp::acceptor::rebind_executor<executor_type>::other;
using websocket_type = beast::websocket::stream<asio::ssl::stream<socket_type>>;
namespace http = beast::http;
// end::decls[]

// tag::connect[]
cobalt::promise<asio::ssl::stream<socket_type>> connect(
        std::string host, boost::asio::ssl::context & ctx)
{
    asio::ip::tcp::resolver res{cobalt::this_thread::get_executor()};
    auto ep = co_await res.async_resolve(host, "https", cobalt::use_op); // <1>

    asio::ssl::stream<socket_type> sock{cobalt::this_thread::get_executor(), ctx};
    co_await sock.next_layer().async_connect(*ep.begin()); // <2>
    co_await sock.async_handshake(asio::ssl::stream_base::client); // <3>

    co_return sock; // <4>
}
// end::connect[]

// tag::ws_upgrade[]
cobalt::promise<void> connect_to_blockchain_info(websocket_type & ws)
{
 ws.set_option(beast::websocket::stream_base::decorator(
     [](beast::websocket::request_type& req)
     {
       req.set(http::field::user_agent,
               std::string(BOOST_BEAST_VERSION_STRING) + " cobalt-ticker");
       req.set(http::field::origin,
               "https://exchange.blockchain.com"); // <1>
     }));

 co_await ws.async_handshake("ws.blockchain.info", "/mercury-gateway/v1/ws"); // <2>
}
// end::ws_upgrade[]

// tag::json_reader[]
cobalt::generator<json::object> json_reader(websocket_type & ws)
try
{
    beast::flat_buffer buf;
    while (ws.is_open()) // <1>
    {
        auto sz = co_await ws.async_read(buf); // <2>
        json::string_view data{static_cast<const char*>(buf.cdata().data()), sz};
        auto obj = json::parse(data);
        co_yield obj.as_object(); // <3>
        buf.consume(sz);
    }
    co_return {};
}
catch (std::exception & e)
{
  std::cerr << "Error reading: " << e.what() << std::endl;
  throw;
}
// end::json_reader[]

// tag::subscription_types[]
using subscription = std::pair<std::string, std::weak_ptr<cobalt::channel<json::object>>>;
using subscription_channel = std::weak_ptr<cobalt::channel<json::object>>;
using subscription_map = boost::unordered_multimap<std::string, subscription_channel>;
// end::subscription_types[]

cobalt::promise<void> handle_rejections(
    std::list<std::string> & unconfirmed,
    subscription_map & subs,
    const json::object & ms)
{
  if (unconfirmed.empty())
    co_return;
  auto rej = unconfirmed.front();
  unconfirmed.pop_front();
  auto r = subs.equal_range(rej);
  for (const auto & [k, chn] : boost::make_iterator_range(r))
    if (auto ptr = chn.lock())
    {
      co_await ptr->write(ms);
      ptr->close();
    }

  subs.erase(r.first, r.second);
}

cobalt::promise<void> handle_update(
   std::list<std::string> & unconfirmed,
   subscription_map & subs,
   const json::object & ms,
   websocket_type & ws)
{

  const auto & sym = json::value_to<std::string>(ms.at("symbol"));

  if (!unconfirmed.empty() && sym == unconfirmed.front())
    unconfirmed.pop_front();

  bool has_expired = false;
  auto r = boost::make_iterator_range(subs.equal_range(sym));
  for (const auto & [k, chn] : r)
    if (auto ptr = chn.lock())
      co_await ptr->write(ms);
    else
      has_expired = true;

  if (has_expired)
    erase_if(subs, [](const auto & p) {return p.second.expired();});

  if (r.empty() && ms.at("event") != "unsubscribed") //
  {
    auto msg = json::serialize(
        json::object{
            {"action", "unsubscribe"},
            {"channel", "ticker"},
            {"symbol", sym}});

    co_await ws.async_write(asio::buffer(msg));
  }
}

cobalt::promise<void> handle_new_subscription(
    std::list<std::string> & unconfirmed,
    subscription_map & subs,
    subscription msg,
    websocket_type & ws)
{
  auto sym = msg.first;
  if (!subs.contains(sym))
  {
    auto msg = json::serialize(
        json::object{
            {"action", "subscribe"},
            {"channel", "ticker"},
            {"symbol", sym}})
    ;
    unconfirmed.push_back(sym);
    co_await ws.async_write(asio::buffer(msg));
  }
  subs.emplace(std::move(msg));
}

// tag::run_blockchain_info[]
cobalt::promise<void> run_blockchain_info(cobalt::channel<subscription> & subc)
try
{
    asio::ssl::context ctx{asio::ssl::context_base::tls_client};
    websocket_type ws{co_await connect("blockchain.info", ctx)};
    co_await connect_to_blockchain_info(ws); // <1>

    subscription_map subs;
    std::list<std::string> unconfirmed;

    auto rd = json_reader(ws); // <2>
    while (ws.is_open()) // <3>
    {
      switch (auto msg = co_await cobalt::race(rd, subc.read()); msg.index()) // <4>
      {
        case 0: // <5>
          if (auto ms = get<0>(msg);
              ms.at("event") == "rejected") // invalid sub, cancel however subbed
            co_await handle_rejections(unconfirmed, subs, ms);
          else
            co_await handle_update(unconfirmed, subs, ms, ws);
        break;
        case 1: // //<6>
            co_await handle_new_subscription(
                unconfirmed, subs,
                std::move(get<1>(msg)), ws);
        break;
      }
    }

    for (auto & [k ,c] : subs)
    {
        if (auto ptr = c.lock())
            ptr->close();
    }
}
catch(std::exception & e)
{
  std::cerr << "Exception: " << e.what() << std::endl;
  throw;
}
// end::run_blockchain_info[]

// tag::read_and_close[]
cobalt::promise<void> read_and_close(beast::websocket::stream<socket_type> & st, beast::flat_buffer buf)
{
    system::error_code ec;
    co_await st.async_read(buf, asio::as_tuple(cobalt::use_op));
    co_await st.async_close(beast::websocket::close_code::going_away, asio::as_tuple(cobalt::use_op));
    st.next_layer().close(ec);
}
// end::read_and_close[]

// tag::run_session[]
cobalt::promise<void> run_session(beast::websocket::stream<socket_type> st,
                                 cobalt::channel<subscription> & subc)
try
{
    http::request<http::empty_body> req;
    beast::flat_buffer buf;
    co_await http::async_read(st.next_layer(), buf, req); // <1>
    // check the target
    auto r = urls::parse_uri_reference(req.target());
    if (r.has_error() || (r->segments().size() != 2u)) // <2>
    {
        http::response<http::string_body> res{http::status::bad_request, 11};
        res.body() = r.has_error() ? r.error().message() :
                    "url needs two segments, e.g. /btc/usd";
        co_await http::async_write(st.next_layer(), res);
        st.next_layer().close();
        co_return ;
    }

    co_await st.async_accept(req); // <3>

    auto sym = std::string(r->segments().front()) + "-" +
               std::string(r->segments().back());
    boost::algorithm::to_upper(sym);
    // close when data gets sent
    auto p = read_and_close(st, std::move(buf)); // <4>

    auto ptr = std::make_shared<cobalt::channel<json::object>>(1u); // <5>
    co_await subc.write(subscription{sym, ptr}); // <6>

    while (ptr->is_open() && st.is_open()) // <7>
    {
      auto bb = json::serialize(co_await ptr->read());
      co_await st.async_write(asio::buffer(bb));
    }

    co_await st.async_close(beast::websocket::close_code::going_away,
                            asio::as_tuple(cobalt::use_op)); // <8>
    st.next_layer().close();
    co_await p; // <9>

}
catch(std::exception & e)
{
    std::cerr << "Session ended with exception: " << e.what() << std::endl;
}
// end::run_session[]

// tag::main[]
cobalt::main co_main(int argc, char * argv[])
{
    acceptor_type acc{co_await cobalt::this_coro::executor,
                      asio::ip::tcp::endpoint (asio::ip::tcp::v4(), 8080)};
    std::cout << "Listening on localhost:8080" << std::endl;

    constexpr int limit = 10; // allow 10 ongoing sessions
    cobalt::channel<subscription> sub_manager; // <1>

    co_await join( // <2>
      run_blockchain_info(sub_manager),
      cobalt::with( // <3>
        cobalt::wait_group(
            asio::cancellation_type::all,
            asio::cancellation_type::all),
        [&](cobalt::wait_group & sessions) -> cobalt::promise<void>
        {
          while (!co_await cobalt::this_coro::cancelled) // <4>
          {
            if (sessions.size() >= limit) // <5>
              co_await sessions.wait_one();

            auto conn = co_await acc.async_accept(); // <6>
            sessions.push_back( // <7>
                run_session(
                    beast::websocket::stream<socket_type>{std::move(conn)},
                    sub_manager));
          }
        })
      );

    co_return 0;
}
// end::main[]
