//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/scheme.hpp>

#include <boost/url/grammar/ci_string.hpp>
#include "test_suite.hpp"

namespace boost {
namespace urls {

class scheme_test
{
public:
    void
    check(
        std::string s0,
        scheme sc0 = scheme::unknown)
    {
        auto sc1 =
            string_to_scheme(s0);
        if(! BOOST_TEST(sc1 == sc0))
            return;
        if(sc0 == scheme::unknown)
            return;
        auto s1 = to_string(sc1);
        for(auto& c : s0)
            c = grammar::to_lower(c);
        BOOST_TEST_EQ(s1, s0);
        BOOST_TEST(s1 != to_string(
            scheme::unknown));
    }

    void
    run()
    {
        // (none)
        check("", scheme::none);

        // ftp
        check("f");
        check("ft");
        check("f0");
        check("ftp0");
        check("ft0");
        check("f00");
        check("ftp", scheme::ftp);
        check("Ftp", scheme::ftp);
        check("FTP", scheme::ftp);

        // file
        check("f");
        check("fi");
        check("fil");
        check("f0");
        check("f00");
        check("f000");
        check("fi00");
        check("fil0");
        check("file0");
        check("file", scheme::file);
        check("File", scheme::file);
        check("FILE", scheme::file);

        // http
        check("h");
        check("ht");
        check("htt");
        check("h0");
        check("h00");
        check("h000");
        check("ht00");
        check("htt");
        check("http0");
        check("http", scheme::http);
        check("htTp", scheme::http);
        check("HTTP", scheme::http);

        // https
        check("h");
        check("ht");
        check("htt");
        check("h0");
        check("h00");
        check("h000");
        check("h0000");
        check("ht000");
        check("htt00");
        check("http0");
        check("https", scheme::https);
        check("htTps", scheme::https);
        check("HTTPS", scheme::https);

        // ws
        check("w");
        check("w0");
        check("w00");
        check("ws0");
        check("ws", scheme::ws);
        check("Ws", scheme::ws);
        check("wS", scheme::ws);
        check("WS", scheme::ws);

        // wss
        check("w");
        check("w0");
        check("w00");
        check("w000");
        check("ws0");
        check("wss0");
        check("wss", scheme::wss);
        check("Wss", scheme::wss);
        check("wSs", scheme::wss);
        check("WSS", scheme::wss);

        // unknown
        check("gopher");
        check("magnet");
        check("mailto");

        BOOST_TEST_EQ(default_port(scheme::none), 0);
        BOOST_TEST_EQ(default_port(scheme::unknown), 0);
        BOOST_TEST_EQ(default_port(scheme::ftp), 21);
        BOOST_TEST_EQ(default_port(scheme::file), 0);
        BOOST_TEST_EQ(default_port(scheme::http), 80);
        BOOST_TEST_EQ(default_port(scheme::https), 443);
        BOOST_TEST_EQ(default_port(scheme::ws), 80);
        BOOST_TEST_EQ(default_port(scheme::wss), 443);
    }
};

TEST_SUITE(
    scheme_test,
    "boost.url.scheme");

} // urls
} // boost
