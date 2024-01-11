//
// Copyright (c) 2022 alandefreitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//

//[example_finicky

/*
    This example shows how to classify URLs
    according to a set of rules. This example is
    inspired by Finicky. The URLs are classified
    and redirected to a browser according to their
    category. See the example config.json file.
    https://github.com/johnste/finicky
*/

#include <boost/url/url.hpp>
#include <boost/url/parse.hpp>
#include <boost/system/result.hpp>
#include <boost/json/stream_parser.hpp>
#include <boost/core/detail/string_view.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>

namespace urls = boost::urls;
namespace json = boost::json;
namespace core = boost::core;

json::value
read_json( std::istream& is, json::error_code& ec )
{
    json::parse_options opt;
    opt.allow_comments = true;
    json::stream_parser p(json::storage_ptr(), opt);
    std::string line;
    while( std::getline( is, line ) )
    {
        p.write( line, ec );
        if( ec )
            return nullptr;
    }
    p.finish( ec );
    if( ec )
        return nullptr;
    return p.release();
}

bool
glob_match(
    core::string_view pattern,
    core::string_view str)
{
    // regex
    if (str.starts_with("/") &&
        str.ends_with("/"))
    {
        const boost::regex pr(pattern.begin() + 1, pattern.end() - 1);
        return boost::regex_match(std::string(str), pr);
    }

    // literal
    if (!pattern.contains('*'))
    {
        return pattern == str;
    }

    // glob
    std::string p = pattern;
    std::size_t i = p.find('*');
    while (i != std::string::npos)
    {
        auto e = std::min(p.find_first_not_of('*', i), p.size());
        std::size_t n = e - i;
        if (n == 1)
        {
            p.replace(i, e, "[^/]*");
            i += 5;
        }
        else
        {
            p.replace(i, e, ".*");
            i += 2;
        }
        i = p.find('*', i);
    }
    const boost::regex pr(p);
    return boost::regex_match(std::string(str), pr);
}

bool
url_match(
    json::value& mv,
    urls::url const& u)
{
    if (mv.is_string())
    {
        json::string& p = mv.as_string();
        return glob_match(u.buffer(), p);
    }
    else if (mv.is_array())
    {
        json::array& m = mv.as_array();
        for (auto& mi: m)
        {
            if (!mi.is_string())
                throw std::invalid_argument(
                    "handle match is not a string");
            if (glob_match(mi.as_string(), u.buffer()))
                return true;
        }
    }
    else if (mv.is_object())
    {
        json::object& m = mv.as_object();
        std::pair<core::string_view, core::string_view>
            field_values[] = {
                {"protocol",  u.scheme()},
                {"authority", u.encoded_authority()},
                {"username",  u.encoded_user()},
                {"user",      u.encoded_user()},
                {"password",  u.encoded_password()},
                {"userinfo",  u.encoded_userinfo()},
                {"host",      u.encoded_host()},
                {"port",      u.port()},
                {"path",      u.encoded_path()},
                {"pathname",  u.encoded_path()},
                {"query",     u.encoded_query()},
                {"search",    u.encoded_query()},
                {"fragment",  u.encoded_fragment()},
                {"hash",      u.encoded_fragment()},
            };
        for (auto& p: field_values)
        {
            auto it = m.find(p.first);
            if (it != m.end())
            {
                if (!it->value().is_string())
                    throw std::invalid_argument(
                        "match fields should be a strings");
                if (glob_match(p.second, p.first))
                    return true;
            }
        }
    }
    return false;
}

#define CHECK(c, msg)             \
    if (!(c))                     \
    {                             \
        std::cerr << msg << "\n"; \
        return EXIT_FAILURE;      \
    }

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cout << argv[0] << "\n";
        std::cout << "Usage: finicky <config> <url>\n"
                     "options:\n"
                     "    <config>: Configuration file\n"
                     "    <url>:    The url to open\n"
                     "examples:\n"
                     "    finicky config.json \"http://www.example.com\"\n";
        return EXIT_FAILURE;
    }

    // Parse url
    boost::system::result<urls::url> ru = urls::parse_uri(argv[2]);
    CHECK(ru, "Invalid URL");
    urls::url u = *ru;

    // Open config file
    std::fstream fin(argv[1]);
    CHECK(fin.good(), "Cannot open configuration file");
    json::error_code ec;
    json::value c = read_json(fin, ec);
    CHECK(!ec.failed(), "Cannot parse configuration file");
    CHECK(c.is_object(), "Configuration file is not an object");
    json::object& o = c.as_object();

    // Set initial browser
    auto bit = o.find("defaultBrowser");
    CHECK(
        bit != o.end(),
        "Configuration file has no defaultBrowser");
    CHECK(
        bit->value().is_string(),
        "defaultBrowser should be a string");
    json::string& browser = bit->value().as_string();

    // Apply rewrites to the input string
    auto rsit = o.find("rewrite");
    if (rsit != o.end())
    {
        CHECK(
            rsit->value().is_array(),
            "rewrite rules should be an array");
        auto& rs = rsit->value().as_array();
        for (auto& rv: rs)
        {
            CHECK(
                rv.is_object(),
                "individual rewrite rule should be an object");
            json::object& r = rv.as_object();

            // Look for match
            auto mit = r.find("match");
            CHECK(
                mit != r.end(),
                "rewrite rule should have a match field");
            CHECK(
                mit->value().is_object() || mit->value().is_string(),
                "rewrite match field is not an object");
            if (!url_match(mit->value(), u))
                continue;

            // Apply replacement rule
            auto uit = r.find("url");
            CHECK(
                uit != r.end(),
                "rewrite rule should have a url field");
            CHECK(
                uit->value().is_object() ||
                uit->value().is_string(),
                "url field must be an object or string");

            if (uit->value().is_string())
            {
                json::string& uo = uit->value().as_string();
                auto ru1 = urls::parse_uri(uo);
                CHECK(ru1, "url " << uo.c_str() << " is invalid");
                u = *ru;
            }
            else
            {
                json::object& uo = uit->value().as_object();
                auto it = uo.find("protocol");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "protocol field should be a string");
                    u.set_scheme(it->value().as_string());
                }

                it = uo.find("authority");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "authority field should be a string");
                    u.set_encoded_authority(
                        it->value().as_string().subview());
                }

                it = uo.find("username");
                if (it == uo.end())
                    it = uo.find("user");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "username field should be a string");
                    u.set_encoded_user(
                        it->value().as_string().subview());
                }

                it = uo.find("password");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "password field should be a string");
                    u.set_encoded_password(
                        it->value().as_string().subview());
                }

                it = uo.find("userinfo");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "userinfo field should be a string");
                    u.set_encoded_userinfo(
                        it->value().as_string().subview());
                }

                it = uo.find("host");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "host field should be a string");
                    u.set_encoded_host(
                        it->value().as_string().subview());
                }

                it = uo.find("port");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "port field should be a string");
                    u.set_port(
                        it->value().as_string().subview());
                }

                it = uo.find("path");
                if (it == uo.end())
                    it = uo.find("pathname");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "path field should be a string");
                    u.set_encoded_path(
                        it->value().as_string().subview());
                }

                it = uo.find("query");
                if (it == uo.end())
                    it = uo.find("search");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "query field should be a string");
                    u.set_encoded_query(
                        it->value().as_string().subview());
                }

                it = uo.find("fragment");
                if (it == uo.end())
                    it = uo.find("hash");
                if (it != uo.end())
                {
                    CHECK(
                        it->value().is_string(),
                        "fragment field should be a string");
                    u.set_encoded_fragment(
                        it->value().as_string().subview());
                }
            }
        }
    }

    // Determine which browser should handle the url
    auto hsit = o.find("handlers");
    if (hsit != o.end())
    {
        CHECK(
            hsit->value().is_array(),
            "handler rules should be an array");
        auto& hs = hsit->value().as_array();
        for (auto& hv: hs)
        {
            CHECK(
                hv.is_object(),
                "individual handlers should be an object");
            json::object& h = hv.as_object();

            auto mit = h.find("match");
            CHECK(
                mit != h.end(),
                "handle rule should have a match field");
            CHECK(
                mit->value().is_string() || mit->value().is_array(),
                "handle match field must be an array or a string");

            auto hbit = h.find("browser");
            CHECK(
                hbit != h.end(),
                "handle rule should have a browser field");
            CHECK(
                hbit->value().is_string(),
                "browser field is not a string");

            // Look for match and change browser
            if (url_match(mit->value(), u))
            {
                browser = hbit->value().as_string().subview();
                break;
            }
        }
    }

    // Print command finicky would run
    std::cout << "\"" << browser.c_str() << "\" " << u << '\n';

    return EXIT_SUCCESS;
}

//]
