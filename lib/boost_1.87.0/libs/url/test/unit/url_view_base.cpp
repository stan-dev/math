//
// Copyright (c) 2022 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/url_view_base.hpp>

#include <boost/url/url_view.hpp>
#include <boost/core/ignore_unused.hpp>

#include "test_suite.hpp"

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

namespace boost {
namespace urls {

struct url_view_base_test
{
    //--------------------------------------------
    //
    // Host
    //
    //--------------------------------------------

    static
    std::string
    bracketed(
        std::string s)
    {
        return
            std::string("[") + s +
            std::string("]");
    }

    void
    testHost()
    {
        auto const ipv4 = [](
            core::string_view s)
        {
            std::string sa = std::string(s);
            std::string su = "https://" + sa + "/";
            url_view u;
            BOOST_TEST_NO_THROW(u = url_view(su));
            BOOST_TEST_EQ(u.host_type(), host_type::ipv4);
            BOOST_TEST_EQ(u.host(),
                u.encoded_host().decode());
            BOOST_TEST_EQ(u.encoded_host(), s);
            BOOST_TEST_EQ(u.host_address(),
                u.encoded_host_address().decode());
            BOOST_TEST_EQ(u.encoded_host_address(), s);
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address(s));
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
            BOOST_TEST_EQ(u.host_ipvfuture(), "");
            BOOST_TEST_EQ(u.host_name(), "");
            BOOST_TEST_EQ(u.encoded_host_name(), "");
            BOOST_TEST_EQ(u.authority().buffer(), sa);
        };

        auto const ipv6 = [](
            core::string_view s)
        {
            std::string sa = bracketed(s);
            std::string su = "https://" + sa + "/";
            url_view u;
            BOOST_TEST_NO_THROW(u = url_view(su));
            BOOST_TEST_EQ(u.host_type(), host_type::ipv6);
            BOOST_TEST_EQ(u.host(),
                u.encoded_host().decode());
            BOOST_TEST_EQ(u.encoded_host(), sa);
            BOOST_TEST_EQ(u.host_address(),
                u.encoded_host_address().decode());
            BOOST_TEST_EQ(u.encoded_host_address(), s);
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address(s));
            BOOST_TEST_EQ(u.host_ipvfuture(), "");
            BOOST_TEST_EQ(u.host_name(), "");
            BOOST_TEST_EQ(u.encoded_host_name(), "");
            BOOST_TEST_EQ(u.authority().buffer(), sa);
        };

        auto const ipvfut = [](
            core::string_view s)
        {
            std::string sa = bracketed(s);
            std::string su = "https://" + sa + "/";
            url_view u;
            BOOST_TEST_NO_THROW(u = url_view(su));
            BOOST_TEST_EQ(u.host_type(), host_type::ipvfuture);
            BOOST_TEST_EQ(u.host(),
                u.encoded_host().decode());
            BOOST_TEST_EQ(u.encoded_host(), sa);
            BOOST_TEST_EQ(u.host_address(),
                u.encoded_host_address().decode());
            BOOST_TEST_EQ(u.encoded_host_address(), s);
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
            BOOST_TEST_EQ(u.host_ipvfuture(), s);
            BOOST_TEST_EQ(u.host_name(), "");
            BOOST_TEST_EQ(u.encoded_host_name(), "");
            BOOST_TEST_EQ(u.authority().buffer(), sa);
        };

        auto const name = [](
            core::string_view s)
        {
            std::string sa = std::string(s);
            std::string su = "https://" + sa + "/";
            url_view u;
            BOOST_TEST_NO_THROW(u = url_view(su));
            BOOST_TEST_EQ(u.host_type(), host_type::name);
            BOOST_TEST_EQ(u.host(),
                u.encoded_host().decode());
            BOOST_TEST_EQ(u.encoded_host(), s);
            BOOST_TEST_EQ(u.host_address(),
                u.encoded_host_address().decode());
            BOOST_TEST_EQ(u.encoded_host_address(), s);
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
            BOOST_TEST_EQ(u.host_ipvfuture(), "");
            BOOST_TEST_EQ(u.host_name(),
                u.encoded_host_name().decode());
            BOOST_TEST_EQ(u.encoded_host_name(), s);
            BOOST_TEST_EQ(u.authority().buffer(), sa);
        };

        ipv4("0.0.0.0");
        ipv4("127.0.0.1");
        ipv4("192.168.0.1");
        ipv4("255.255.255.255");

        ipv6("1::6:192.168.0.1");

        ipvfut("v1.x");

        name("www.example.com");
        name("www%2eexample%2ecom");

        BOOST_TEST(url_view().encoded_host_address().empty());
    }

    void
    testJavadocs()
    {
        //----------------------------------------
        //
        // Observers
        //
        //----------------------------------------

        // size
        {
        assert( url_view( "file:///Program%20Files" ).size() == 23 );
        }

        // empty
        {
        assert( url_view( "" ).empty() );
        }

        // persist
        {
            std::shared_ptr< url_view const > sp;
            {
                std::string s( "http://example.com" );
                url_view u( s );                        // u references characters in s

                assert( u.data() == s.data() );         // same buffer

                sp = u.persist();

                assert( sp->data() != s.data() );       // different buffer
                assert( sp->buffer() == s);             // same contents

                // s is destroyed and thus u
                // becomes invalid, but sp remains valid.
            }
        }

        //----------------------------------------
        //
        // Scheme
        //
        //----------------------------------------

        // has_scheme
        assert( url_view( "http://www.example.com" ).has_scheme() );

        // scheme
        assert( url_view( "http://www.example.com" ).scheme() == "http" );

        // scheme_id
        assert( url_view( "wss://www.example.com/crypto.cgi" ).scheme_id() == scheme::wss );

        //----------------------------------------
        //
        // Authority
        //
        //----------------------------------------

        // has_authority
        assert( url_view( "http://www.example.com/index.htm" ).has_authority() );

        // authority
        {
        authority_view a = url_view( "https://www.example.com:8080/index.htm" ).authority();

        ignore_unused(a);
        }

        // encoded_authority
        assert( url_view( "file://Network%20Drive/My%2DFiles" ).encoded_authority() == "Network%20Drive" );

        //----------------------------------------
        //
        // Userinfo
        //
        //----------------------------------------

        // has_userinfo
        assert( url_view( "http://jane%2Ddoe:pass@example.com" ).has_userinfo() );

        // userinfo
        assert( url_view( "http://jane%2Ddoe:pass@example.com" ).userinfo() == "jane-doe:pass" );

        // encoded_userinfo
        assert( url_view( "http://jane%2Ddoe:pass@example.com" ).encoded_userinfo() == "jane%2Ddoe:pass" );

        // user
        assert( url_view( "http://jane%2Ddoe:pass@example.com" ).user() == "jane-doe" );

        // encoded_user
        assert( url_view( "http://jane%2Ddoe:pass@example.com" ).encoded_user() == "jane%2Ddoe" );

        // has_password
        assert( url_view( "http://jane%2Ddoe:pass@example.com" ).has_password() );

        // password
        assert( url_view( "http://jane%2Ddoe:pass@example.com" ).password() == "pass" );

        // encoded_password
        assert( url_view( "http://jane%2Ddoe:pass@example.com" ).encoded_password() == "pass" );

        //----------------------------------------
        //
        // Host
        //
        //----------------------------------------

        // host_type
        assert( url_view( "https://192.168.0.1/local.htm" ).host_type() == host_type::ipv4 );

        // host
        assert( url_view( "https://www%2droot.example.com/" ).host() == "www-root.example.com" );

        // encoded_host
        assert( url_view( "https://www%2droot.example.com/" ).encoded_host() == "www%2droot.example.com" );

        // host_address
        assert( url_view( "https://[1::6:c0a8:1]/" ).host_address() == "1::6:c0a8:1" );

        // encoded_host_address
        assert( url_view( "https://www%2droot.example.com/" ).encoded_host_address() == "www%2droot.example.com" );

        // ipv4_address
        assert( url_view( "http://127.0.0.1/index.htm?user=win95" ).host_ipv4_address() == ipv4_address( "127.0.0.1" ) );

        // ipv6_address
        assert( url_view( "ftp://[::1]/" ).host_ipv6_address() == ipv6_address( "::1" ) );

        // ipvfuture
        assert( url_view( "http://[v1fe.d:9]/index.htm" ).host_ipvfuture() == "v1fe.d:9" );

        // host_name
        assert( url_view( "https://www%2droot.example.com/" ).host_name() == "www-root.example.com" );

        // encoded_host_name
        assert( url_view( "https://www%2droot.example.com/" ).encoded_host_name() == "www%2droot.example.com" );

        //----------------------------------------
        //
        // Port
        //
        //----------------------------------------

        // has_port
        assert( url_view( "wss://www.example.com:443" ).has_port() );

        // port
        assert( url_view( "http://localhost.com:8080" ).port() == "8080" );

        // port_number
        assert( url_view( "http://localhost.com:8080" ).port_number() == 8080 );

        //----------------------------------------
        //
        // Path
        //
        //----------------------------------------

        // is_path_absolute
        assert( url_view( "/path/to/file.txt" ).is_path_absolute() );

        // path
        assert( url_view( "file:///Program%20Files/Games/config.ini" ).path() == "/Program Files/Games/config.ini" );

        // encoded_path
        assert( url_view( "file:///Program%20Files/Games/config.ini" ).encoded_path() == "/Program%20Files/Games/config.ini" );

        // segments
        {
        segments_view sv = url_view( "/path/to/file.txt" ).segments();

        ignore_unused(sv);
        }

        // encoded_segments
        {
        segments_encoded_view sv = url_view( "/path/to/file.txt" ).encoded_segments();

        ignore_unused(sv);
        }

        //----------------------------------------
        //
        // Query
        //
        //----------------------------------------

        // has_query
        assert( url_view( "/sql?id=42&col=name&page-size=20" ).has_query() );

        // query
        assert( url_view( "/sql?id=42&name=jane%2Ddoe&page+size=20" ).query() == "id=42&name=jane-doe&page+size=20" );

        // encoded_query
        assert( url_view( "/sql?id=42&name=jane%2Ddoe&page+size=20" ).encoded_query() == "id=42&name=jane%2Ddoe&page+size=20" );

        // params
        {
        params_view pv = url_view( "/sql?id=42&name=jane%2Ddoe&page+size=20" ).params();

        ignore_unused(pv);
        }

        // encoded_params
        {
        params_encoded_view pv = url_view( "/sql?id=42&name=jane%2Ddoe&page+size=20" ).encoded_params();

        ignore_unused(pv);
        }

        //----------------------------------------
        //
        // Fragment
        //
        //----------------------------------------

        // has_fragment
        assert( url_view( "http://www.example.com/index.htm#a%2D1" ).has_fragment() );

        // fragment
        assert( url_view( "http://www.example.com/index.htm#a%2D1" ).fragment() == "a-1" );

        // encoded_fragment
        assert( url_view( "http://www.example.com/index.htm#a%2D1" ).encoded_fragment() == "a%2D1" );

        //----------------------------------------
        //
        // Compound Fields
        //
        //----------------------------------------

        // encoded_host_and_port
        assert( url_view( "http://www.example.com:8080/index.htm" ).encoded_host_and_port() == "www.example.com:8080" );

        // encoded_origin
        assert( url_view( "http://www.example.com:8080/index.htm?text=none#h1" ).encoded_origin() == "http://www.example.com:8080" );

        // encoded_target
        assert( url_view( "http://www.example.com/index.html?query#frag" ).encoded_target() == "/index.html?query" );

        // encoded_resource
        assert( url_view( "http://www.example.com/index.html?query#frag" ).encoded_resource() == "/index.html?query#frag" );

        //----------------------------------------
        //
        // Comparison
        //
        //----------------------------------------

        // url_view_base::op<<
        {
            url_view u( "http://www.example.com/index.htm" );
            std::stringstream ss;
            ss << u;
            assert( ss.str() == "http://www.example.com/index.htm" );
            ignore_unused(ss);
        }

        // url_view_base::op==
        {
            url_view u0( "http://www.a.com/index.htm" );
            url_view u1( "http://www.a.com/in%64ex.htm" );
            assert( u0 == u1 );
            ignore_unused(u0, u1);
        }

        // url_view_base::op!=
        {
            url_view u0( "http://www.a.com/index.htm" );
            url_view u1( "http://www.b.com/index.htm" );
            assert( u0 != u1 );
            ignore_unused(u0, u1);
        }

        // url_view_base::op<
        {
            url_view u0( "http://www.a.com/index.htm" );
            url_view u1( "http://www.b.com/index.htm" );
            assert( u0 < u1 );
            ignore_unused(0, u1);
        }

        // url_view_base::op<=
        {
            url_view u0( "http://www.b.com/index.htm" );
            url_view u1( "http://www.b.com/index.htm" );
            assert( u0 <= u1 );
            ignore_unused(u0, u1);
        }

        // url_view_base::op>
        {
            url_view u0( "http://www.b.com/index.htm" );
            url_view u1( "http://www.a.com/index.htm" );
            assert( u0 > u1 );
            ignore_unused(u0, u1);
        }

        // url_view_base::op>=
        {
            url_view u0( "http://www.a.com/index.htm" );
            url_view u1( "http://www.a.com/index.htm" );
            assert( u0 >= u1 );
            ignore_unused(u0, u1);
        }

    }

    void
    run()
    {
        testHost();
        testJavadocs();

        test_suite::log <<
            "    sizeof(url_impl) == " <<
            sizeof(detail::url_impl) << "\n";
    }
};

TEST_SUITE(
    url_view_base_test,
    "boost.url.url_view_base");

} // urls
} // boost
