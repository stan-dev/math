//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// tag::example_router[]

/*
    This example defines a router for URL paths.
    Each path is associated with a callback
    function.
*/

#ifndef BOOST_URL_SOURCE
#define BOOST_URL_SOURCE
#endif

#include "router.hpp"

#include <boost/beast/core.hpp>
#include <boost/beast/http.hpp>
#include <boost/beast/version.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/config.hpp>

#include <iostream>
#include <functional>

namespace urls = boost::urls;
namespace core = boost::core;
namespace asio = boost::asio;
namespace beast = boost::beast;
namespace http = beast::http;
using string_view = core::string_view;
using request_t = http::request<http::string_body>;
struct connection;
using handler = std::function<void(connection&, urls::matches)>;

int
serve(
    urls::router<handler> const& r,
    asio::ip::address const& a,
    unsigned short port,
    std::string const& doc_root);

struct connection
{
    connection(asio::io_context& ioc)
        : socket(ioc) {}

    void
    string_reply(core::string_view msg);

    void
    file_reply(core::string_view path);

    void
    error_reply(http::status, core::string_view msg);

    beast::error_code ec;
    asio::ip::tcp::socket socket;
    std::string doc_root;
    request_t req;
};

int
main(int argc, char **argv)
{
    /*
     * Parse cmd-line params
     */
    if (argc != 4)
    {
        core::string_view exec = argv[0];
        auto file_pos = exec.find_last_of("/\\");
        if (file_pos != core::string_view::npos)
            exec = exec.substr(file_pos + 1);
        std::cerr
            << "Usage: " << exec
            << " <address> <port> <doc_root>\n"
               "Example: " << exec << " 0.0.0.0 8080 .\n"
                       "Default values:\n"
                       "- address: 0.0.0.0\n"
                       "- port: 8080\n"
                       "- doc_root: ./\n";
    }
    auto const address = asio::ip::make_address(argc > 1 ? argv[1] : "0.0.0.0");
    auto const port = static_cast<unsigned short>(argc > 2 ? std::atoi(argv[2]) : 8080);
    auto const doc_root = std::string(argc > 3 ? argv[3] : ".");

    /*
     * Create router
     */
    urls::router<handler> r;

    r.insert("/", [&](connection& c, urls::matches const&) {
        c.string_reply("Hello!");
    });

    r.insert("/user/{name}", [&](connection& c, urls::matches const& m) {
        std::string msg = "Hello, ";
        urls::pct_string_view(m[0]).decode({}, urls::string_token::append_to(msg));
        msg += "!";
        c.string_reply(msg);
    });

    r.insert("/user", [&](connection& c, urls::matches const&) {
        std::string msg = "Users: ";
        auto names = {"johndoe", "maria", "alice"};
        for (auto name: names) {
            msg += "<a href=\"/user/";
            msg += name;
            msg += "\">";
            msg += name;
            msg += "</a> ";
        }
        c.string_reply(msg);
    });

    r.insert("/public/{path+}", [&](connection& c, urls::matches m) {
        c.file_reply(m["path"]);
    });

    return serve(r, address, port, doc_root);
}

#define ROUTER_CHECK(cond)       if(!(cond)) { break; }
#define ROUTER_CHECK_EC(ec, cat) if(ec.failed()) { std::cerr << #cat << ": " << ec.message() << "\n"; break; }

int
serve(
    urls::router<handler> const& r,
    asio::ip::address const& address,
    unsigned short port,
    std::string const& doc_root)
{
    /*
     * Serve the routes with a simple synchronous
     * server. This is an implementation detail
     * in the context of this example.
     */
    std::cout << "Listening on http://" << address << ":" << port << "\n";
    asio::io_context ioc(1);
    asio::ip::tcp::acceptor acceptor(ioc, {address, port});
    urls::matches m;
    for(;;)
    {
        connection c(ioc);
        c.doc_root = doc_root;
        acceptor.accept(c.socket);
        beast::flat_buffer buffer;
        for(;;)
        {
            // Read a request
            http::read(c.socket, buffer, c.req, c.ec);
            ROUTER_CHECK(c.ec != http::error::end_of_stream)
            ROUTER_CHECK_EC(c.ec, read)
            // Handle request
            auto rpath = urls::parse_path(c.req.target());
            if (c.req.method() != http::verb::get &&
                c.req.method() != http::verb::head)
                c.error_reply(
                    http::status::bad_request,
                    std::string("Unknown HTTP-method: ") +
                        std::string(c.req.method_string()));
            else if (!rpath)
                c.error_reply(http::status::bad_request, "Illegal request-target");
            else if (auto h = r.find(*rpath, m))
                (*h)(c, m);
            else
                c.error_reply(
                    http::status::not_found,
                    "The resource '" +
                        std::string(rpath->buffer()) +
                        "' was not found.");
            ROUTER_CHECK_EC(c.ec, write)
            ROUTER_CHECK(c.req.keep_alive())
        }
        c.socket.shutdown(asio::ip::tcp::socket::shutdown_send, c.ec);
    }
    return EXIT_SUCCESS;
}

#undef ROUTER_CHECK_EC
#undef ROUTER_CHECK

void
connection::
error_reply(http::status s, core::string_view msg)
{
    // invalid route
    http::response<http::string_body> res{s, req.version()};
    res.set(http::field::server, BOOST_BEAST_VERSION_STRING);
    res.set(http::field::content_type, "text/html");
    res.keep_alive(req.keep_alive());
    res.body() = msg;
    res.prepare_payload();
    http::write(socket, res, ec);
}


void
connection::
string_reply(core::string_view msg)
{
    http::response<http::string_body> res{http::status::ok, req.version()};
    res.set(http::field::server, BOOST_BEAST_VERSION_STRING);
    res.set(http::field::content_type, "text/html");
    res.keep_alive(req.keep_alive());
    res.body() = msg;
    res.prepare_payload();
    http::write(socket, res, ec);
}

core::string_view
mime_type(core::string_view path);

std::string
path_cat(
    beast::string_view base,
    beast::string_view path);

void
connection::
file_reply(core::string_view path)
{
    http::file_body::value_type body;
    std::string jpath = path_cat(doc_root, path);
    body.open(jpath.c_str(), beast::file_mode::scan, ec);
    if(ec == beast::errc::no_such_file_or_directory)
    {
        error_reply(
            http::status::not_found,
            "The resource '" + std::string(path) +
                "' was not found in " + jpath);
        return;
    }
    auto const size = body.size();
    http::response<http::file_body> res{
        std::piecewise_construct,
        std::make_tuple(std::move(body)),
        std::make_tuple(http::status::ok, req.version())};
    res.set(http::field::server, BOOST_BEAST_VERSION_STRING);
    res.set(http::field::content_type, mime_type(path));
    res.content_length(size);
    res.keep_alive(req.keep_alive());
    http::write(socket, res, ec);
}

// Append an HTTP rel-path to a local filesystem path.
// The returned path is normalized for the platform.
std::string
path_cat(
    core::string_view base,
    core::string_view path)
{
    if (base.empty())
        return std::string(path);
    std::string result(base);
#ifdef BOOST_MSVC
    char constexpr path_separator = '\\';
#else
    char constexpr path_separator = '/';
#endif
    if( result.back() == path_separator &&
        path.starts_with(path_separator))
        result.resize(result.size() - 1);
    else if (result.back() != path_separator &&
             !path.starts_with(path_separator))
    {
        result.push_back(path_separator);
    }
    result.append(path.data(), path.size());
#ifdef BOOST_MSVC
    for(auto& c : result)
        if(c == '/')
            c = path_separator;
#endif
    return result;
}

core::string_view
mime_type(core::string_view path)
{
    using beast::iequals;
    auto const ext = [&path]
    {
        auto const pos = path.rfind(".");
        if(pos == beast::string_view::npos)
            return beast::string_view{};
        return path.substr(pos);
    }();
    if(iequals(ext, ".htm"))  return "text/html";
    if(iequals(ext, ".html")) return "text/html";
    if(iequals(ext, ".php"))  return "text/html";
    if(iequals(ext, ".css"))  return "text/css";
    if(iequals(ext, ".txt"))  return "text/plain";
    if(iequals(ext, ".js"))   return "application/javascript";
    if(iequals(ext, ".json")) return "application/json";
    if(iequals(ext, ".xml"))  return "application/xml";
    if(iequals(ext, ".swf"))  return "application/x-shockwave-flash";
    if(iequals(ext, ".flv"))  return "video/x-flv";
    if(iequals(ext, ".png"))  return "image/png";
    if(iequals(ext, ".jpe"))  return "image/jpeg";
    if(iequals(ext, ".jpeg")) return "image/jpeg";
    if(iequals(ext, ".jpg"))  return "image/jpeg";
    if(iequals(ext, ".gif"))  return "image/gif";
    if(iequals(ext, ".bmp"))  return "image/bmp";
    if(iequals(ext, ".ico"))  return "image/vnd.microsoft.icon";
    if(iequals(ext, ".tiff")) return "image/tiff";
    if(iequals(ext, ".tif"))  return "image/tiff";
    if(iequals(ext, ".svg"))  return "image/svg+xml";
    if(iequals(ext, ".svgz")) return "image/svg+xml";
    return "application/text";
}

// end::example_router[]
