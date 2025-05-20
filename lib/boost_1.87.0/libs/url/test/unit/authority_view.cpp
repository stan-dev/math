//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/authority_view.hpp>

#include <boost/url/grammar/parse.hpp>
#include "test_rule.hpp"
#include <sstream>

namespace boost {
namespace urls {

class authority_view_test
{
public:
    void
    testSpecialMembers()
    {
        // authority_view()
        // ~authority_view()
        {
            authority_view a;
            BOOST_TEST(a.empty());
            BOOST_TEST_EQ(a.size(), 0u);
        }
    }

    void
    testObservers()
    {
        // size()
        {
            authority_view a;
            BOOST_TEST_EQ(a.size(), 0u);
            a = parse_authority("x").value();
            BOOST_TEST_EQ(a.size(), 1u);
        }

        // empty()
        {
            authority_view a;
            BOOST_TEST(a.empty());
            a = parse_authority("x").value();
            BOOST_TEST(! a.empty());
        }

        // data()
        {
            core::string_view s = "x";
            authority_view a = parse_authority(s).value();
            BOOST_TEST_NE(a.data(), nullptr);
            BOOST_TEST_EQ(a.data(), s.data());
        }

        // string()
        {
            core::string_view s = "xyz";
            authority_view a = parse_authority(s).value();
            BOOST_TEST_EQ(a.buffer(), s);
            BOOST_TEST_EQ(a.buffer().data(), s.data());
        }
    }

    void
    testUserinfo()
    {
        auto const yes =
            []( core::string_view s,
                core::string_view m1,
                core::string_view m2)
        {
            BOOST_TEST_NO_THROW(
            [&]{
                auto a = parse_authority(s).value();
                BOOST_TEST(a.has_userinfo());
                BOOST_TEST(
                    a.encoded_userinfo() == m1);
                BOOST_TEST(
                    a.userinfo() == m2);
            }());
        };

        yes("a@", "a", "a");
        yes(":@", ":", ":");
        yes("@", "", "");
        yes("@x", "", "");
        yes("%61@x", "%61", "a");
        yes(":%61@x", ":%61", ":a");
        yes("%61%3a%62@x", "%61%3a%62", "a:b");

        // issue 828
        {
            auto a = parse_authority("").value();
            BOOST_TEST_NOT(a.has_userinfo());
            BOOST_TEST(a.encoded_userinfo() == "");
            BOOST_TEST(a.userinfo() == "");
        }

        {
            auto a = parse_authority("@").value();
            BOOST_TEST(a.has_userinfo());
            BOOST_TEST_EQ(a.encoded_userinfo(), "");
            BOOST_TEST_EQ(a.userinfo(), "");
            BOOST_TEST_EQ(a.encoded_user(), "");
            BOOST_TEST_EQ(a.user(), "");
            BOOST_TEST_EQ(a.has_password(), false);
            BOOST_TEST_EQ(a.encoded_password(), "");
            BOOST_TEST_EQ(a.password(), "");
        }
        {
            auto a = parse_authority(":@").value();
            BOOST_TEST(a.has_userinfo());
            BOOST_TEST_EQ(a.encoded_userinfo(), ":");
            BOOST_TEST_EQ(a.userinfo(), ":");
            BOOST_TEST_EQ(a.encoded_user(), "");
            BOOST_TEST_EQ(a.user(), "");
            BOOST_TEST_EQ(a.has_password(), true);
            BOOST_TEST_EQ(a.encoded_password(), "");
            BOOST_TEST_EQ(a.password(), "");
        }
        {
            auto a = parse_authority("a%41:@").value();
            BOOST_TEST(a.has_userinfo());
            BOOST_TEST_EQ(a.encoded_userinfo(), "a%41:");
            BOOST_TEST_EQ(a.encoded_user(), "a%41");
            BOOST_TEST_EQ(a.user(), "aA");
            BOOST_TEST_EQ(a.has_password(), true);
            BOOST_TEST_EQ(a.encoded_password(), "");
            BOOST_TEST_EQ(a.password(), "");
        }
        {
            auto a = parse_authority(":b%42@").value();
            BOOST_TEST(a.has_userinfo());
            BOOST_TEST_EQ(a.encoded_userinfo(), ":b%42");
            BOOST_TEST_EQ(a.encoded_user(), "");
            BOOST_TEST_EQ(a.user(), "");
            BOOST_TEST_EQ(a.has_password(), true);
            BOOST_TEST_EQ(a.encoded_password(), "b%42");
            BOOST_TEST_EQ(a.password(), "bB");
        }
        {
            auto a = parse_authority("a:b@").value();
            BOOST_TEST(a.has_userinfo());
            BOOST_TEST_EQ(a.encoded_userinfo(), "a:b");
            BOOST_TEST_EQ(a.encoded_user(), "a");
            BOOST_TEST_EQ(a.has_password(), true);
            BOOST_TEST_EQ(a.encoded_password(), "b");
        }
        {
            auto a = parse_authority("%3a:%3a@").value();
            BOOST_TEST(a.has_userinfo());
            BOOST_TEST_EQ(a.encoded_userinfo(), "%3a:%3a");
            BOOST_TEST_EQ(a.userinfo(), ":::");
            BOOST_TEST_EQ(a.encoded_user(), "%3a");
            BOOST_TEST_EQ(a.user(), ":");
            BOOST_TEST_EQ(a.has_password(), true);
            BOOST_TEST_EQ(a.encoded_password(), "%3a");
            BOOST_TEST_EQ(a.password(), ":");
        }
        {
            auto a = parse_authority("%2525@").value();
            BOOST_TEST(a.has_userinfo());
            BOOST_TEST_EQ(a.encoded_userinfo(), "%2525");
            BOOST_TEST_EQ(a.userinfo(), "%25");
            BOOST_TEST_EQ(a.encoded_user(), "%2525");
            BOOST_TEST_EQ(a.user(), "%25");
            BOOST_TEST_EQ(a.has_password(), false);
            BOOST_TEST_EQ(a.encoded_password(), "");
            BOOST_TEST_EQ(a.password(), "");
        }
    }

    void
    testHost()
    {
        {
            auto a = parse_authority("").value();
            BOOST_TEST(a.host_type() ==
                host_type::name);
            BOOST_TEST(a.encoded_host() ==
                "");
            BOOST_TEST(a.encoded_host_address() ==
                "");
            BOOST_TEST(a.encoded_host_name() == "");
            BOOST_TEST(a.host_ipv4_address()
                == ipv4_address());
            BOOST_TEST(a.host_ipv6_address()
                == ipv6_address());
            BOOST_TEST(
                a.host_ipvfuture() == "");
        }
        {
            auto a = parse_authority("").value();
            BOOST_TEST(a.host_type() ==
                host_type::name);
            BOOST_TEST(a.encoded_host() ==
                "");
            BOOST_TEST(a.encoded_host_address() ==
                "");
            BOOST_TEST(a.encoded_host_name() == "");
        }
        {
            auto a = parse_authority("").value();
            BOOST_TEST(a.host_type() ==
                host_type::name);
            BOOST_TEST(a.encoded_host() ==
                "");
            BOOST_TEST(a.encoded_host_address() ==
                "");
        }
        {
            auto a = parse_authority("www.example.com").value();
            BOOST_TEST(a.host_type() ==
                host_type::name);
            BOOST_TEST(a.encoded_host() ==
                "www.example.com");
            BOOST_TEST(a.encoded_host_address() ==
                "www.example.com");
            BOOST_TEST(a.host() ==
                "www.example.com");
        }
        {
            auto a = parse_authority("192.168.0.1").value();
            BOOST_TEST(a.host_type() ==
                host_type::ipv4);
            BOOST_TEST(a.encoded_host() ==
                "192.168.0.1");
            BOOST_TEST(a.encoded_host_address() ==
                "192.168.0.1");
            BOOST_TEST(a.encoded_host_name() == "");
            BOOST_TEST(a.host() ==
                "192.168.0.1");
            BOOST_TEST(
                a.host_ipv4_address().to_uint() ==
                    0xc0a80001);
        }
        {
            auto a = parse_authority(
                "[1::6:192.168.0.1]:8080").value();
            BOOST_TEST(a.host_type() ==
                host_type::ipv6);
            BOOST_TEST(a.encoded_host() ==
                "[1::6:192.168.0.1]");
            BOOST_TEST(a.encoded_host_address() ==
                "1::6:192.168.0.1");
            BOOST_TEST(a.host() ==
                "[1::6:192.168.0.1]");
            BOOST_TEST(a.host_ipv6_address() ==
                ipv6_address("1::6:c0a8:1"));
        }
        {
            auto a = parse_authority(
                "[v1.x]:8080").value();
            BOOST_TEST(a.host_type() ==
                host_type::ipvfuture);
            BOOST_TEST(a.encoded_host() ==
                "[v1.x]");
            BOOST_TEST(a.encoded_host_address() ==
                "v1.x");
            BOOST_TEST(a.host() ==
                "[v1.x]");
            BOOST_TEST(a.host_ipvfuture() == "v1.x");
        }
    }

    void
    testPort()
    {
        {
            auto a = parse_authority("").value();
            BOOST_TEST(! a.has_port());
            BOOST_TEST_EQ(a.port(), "");
            BOOST_TEST_EQ(a.port_number(), 0);
        }
        {
            auto a = parse_authority("www").value();
            BOOST_TEST(! a.has_port());
            BOOST_TEST_EQ(a.port(), "");
            BOOST_TEST_EQ(a.port_number(), 0);
        }
        {
            auto a = parse_authority(":").value();
            BOOST_TEST(a.has_port());
            BOOST_TEST_EQ(a.port(), "");
            BOOST_TEST_EQ(a.port_number(), 0);
        }
        {
            auto a = parse_authority(":0").value();
            BOOST_TEST(a.has_port());
            BOOST_TEST_EQ(a.port(), "0");
            BOOST_TEST_EQ(a.port_number(), 0);
        }
        {
            auto a = parse_authority(":42").value();
            BOOST_TEST(a.has_port());
            BOOST_TEST_EQ(a.port(), "42");
            BOOST_TEST_EQ(a.port_number(), 42);
        }
        {
            auto a = parse_authority(":00000").value();
            BOOST_TEST(a.has_port());
            BOOST_TEST_EQ(a.port(), "00000");
            BOOST_TEST_EQ(a.port_number(), 0);
        }
        {
            auto a = parse_authority(":000001").value();
            BOOST_TEST(a.has_port());
            BOOST_TEST_EQ(a.port(), "000001");
            BOOST_TEST_EQ(a.port_number(), 1);
        }
        {
            auto a = parse_authority(":65535").value();
            BOOST_TEST(a.has_port());
            BOOST_TEST_EQ(a.port(), "65535");
            BOOST_TEST_EQ(a.port_number(), 65535);
        }
        {
            auto a = parse_authority(":65536").value();
            BOOST_TEST(a.has_port());
            BOOST_TEST_EQ(a.port(), "65536");
            BOOST_TEST_EQ(a.port_number(), 0);
        }
    }

    void
    testHostAndPort()
    {
        {
            auto a = parse_authority("x:1").value();
            BOOST_TEST(a.encoded_host_and_port() ==
                "x:1");
        }
        {
            auto a = parse_authority("x%3a:1").value();
            BOOST_TEST(a.encoded_host_and_port() ==
                "x%3a:1");
        }
        {
            auto a = parse_authority(":1").value();
            BOOST_TEST(a.encoded_host_and_port() ==
                ":1");
        }
        {
            auto a = parse_authority(":000001").value();
            BOOST_TEST(a.encoded_host_and_port() ==
                ":000001");
        }
        {
            auto a = parse_authority("xyz:99999").value();
            BOOST_TEST(a.encoded_host_and_port() ==
                "xyz:99999");
        }
    }

    void
    testOStream()
    {
        authority_view a( "user:pass@www.example.com:8080" );
        std::ostringstream os;
        os << a;
        BOOST_TEST_EQ( os.str(), "user:pass@www.example.com:8080" );
    }

    void
    run()
    {
        // javadocs
        // authority_view
        {
            authority_view a( "user:pass@www.example.com:8080" );
            (void)a;
        }
        {
            system::result< authority_view > rv = parse_authority( "user:pass@www.example.com:8080" );
            (void)rv;
        }
        {
            BOOST_TEST( parse_authority( "www.example.com" ).value().buffer() == "www.example.com" );
        }

        testSpecialMembers();
        testObservers();
        testUserinfo();
        testHost();
        testPort();
        testHostAndPort();
        testOStream();
    }
};

TEST_SUITE(
    authority_view_test,
    "boost.url.authority_view");

} // urls
} // boost
