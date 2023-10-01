//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/ipv4_address.hpp>

#include "test_suite.hpp"
#include <sstream>

namespace boost {
namespace urls {

class ipv4_address_test
{
public:
    void
    testMembers()
    {
        // ipv4_address()
        {
            ipv4_address a;
            BOOST_TEST(a.is_unspecified());
            BOOST_TEST_EQ(a, ipv4_address());
        }

        // ipv4_address(ipv4_address const&)
        {
            ipv4_address a1(1);
            ipv4_address a2(a1);
            BOOST_TEST_EQ(a2, a1);
        }

        // operator=(ipv4_address const&)
        {
            ipv4_address a1(1);
            ipv4_address a2;
            BOOST_TEST_NE(a2, a1);
            a2 = a1;
            BOOST_TEST_EQ(a2, a1);
        }

        // ipv4_address(array)
        {
            {
                ipv4_address a(
                    ipv4_address::bytes_type{{
                        1,2,3,4}});
                BOOST_TEST(
                    a.to_uint()==0x01020304);
            }
            {
                ipv4_address a(
                    ipv4_address::bytes_type{
                        {1,2,3,4}});
                BOOST_TEST(
                    a.to_uint()==0x01020304);
            }
        }

        // ipv4_address(uint_type)
        {
            ipv4_address a(0x01020304);
            BOOST_TEST(
                a.to_uint() == 0x01020304);
        }

        // ipv4_address(string_view)
        {
            ipv4_address a("1.2.3.4");
            BOOST_TEST(
                a.to_uint() == 0x01020304);
            BOOST_TEST_THROWS(
                ipv4_address("x"),
                system_error);
        }

        // to_bytes
        {
            ipv4_address a(0x01020304);
            ipv4_address::bytes_type b =
                {{1, 2, 3, 4}};
            BOOST_TEST(a.to_bytes() == b);
        }

        // to_uint
        {
            ipv4_address a(0x01020304);
            BOOST_TEST(
                a.to_uint() == 0x01020304);
        }

        // to_string
        {
            ipv4_address a(0x01020304);
            BOOST_TEST(
                a.to_string() == "1.2.3.4");
        }

        // to_buffer
        {
            ipv4_address a(0x01020304);
            char buf[ipv4_address::max_str_len];
            BOOST_TEST(a.to_buffer(buf,
                sizeof(buf)) == "1.2.3.4");
            char buf2[10];
            BOOST_TEST_THROWS(
                a.to_buffer(buf2, sizeof(buf2)),
                system_error);
        }

        // is_loopback
        {
            BOOST_TEST(ipv4_address(
                "127.0.0.1").is_loopback());
        }

        // is_unspecified
        {
            BOOST_TEST(
                ipv4_address().is_unspecified());
            BOOST_TEST(
                ipv4_address(0).is_unspecified());
            BOOST_TEST(ipv4_address("0.0.0.0"
                ).is_unspecified());
        }

        // is_multicast
        {
            BOOST_TEST(ipv4_address(
                "224.0.0.1").is_multicast());
        }

        // operator==
        // operator!=
        {
            ipv4_address a1(1);
            ipv4_address a2(2);
            ipv4_address a3(1);
            BOOST_TEST_EQ(a1, a1);
            BOOST_TEST_NE(a1, a2);
            BOOST_TEST_EQ(a1, a3);
            BOOST_TEST_NE(a2, a3);
        }

        // static any()
        {
            BOOST_TEST(ipv4_address::any(
                ).is_unspecified());
            BOOST_TEST(! ipv4_address::any(
                ).is_loopback());
            BOOST_TEST(! ipv4_address::any(
                ).is_multicast());
        }

        // static loopback()
        {
            BOOST_TEST(ipv4_address::loopback(
                ).is_loopback());
            BOOST_TEST(! ipv4_address::loopback(
                ).is_unspecified());
            BOOST_TEST(! ipv4_address::loopback(
                ).is_multicast());
        }

        // static broadcast()
        {
            BOOST_TEST(! ipv4_address::broadcast(
                ).is_loopback());
            BOOST_TEST(! ipv4_address::broadcast(
                ).is_unspecified());
            BOOST_TEST(! ipv4_address::broadcast(
                ).is_multicast());
        }

        // operator<<
        {
            std::stringstream ss;
            ss << ipv4_address(0x01020304);
            BOOST_TEST_EQ(ss.str(), "1.2.3.4");
        }
    }

    void
    testParse()
    {
        auto const bad = [](
            string_view s)
        {
            auto r = parse_ipv4_address(s);
            BOOST_TEST(r.has_error());
        };

        auto const good = [](
            string_view s)
        {
            auto r = parse_ipv4_address(s);
            BOOST_TEST(r.has_value());
        };

        auto const check = [](
            string_view s,
            ipv4_address::uint_type v)
        {
            auto r = parse_ipv4_address(s);
            if(! BOOST_TEST(r.has_value()))
                return;
            BOOST_TEST_EQ(r->to_uint(), v);
        };

        bad("0");
        bad("0.");
        bad("0.0");
        bad("0.0.");
        bad("0.0.0");
        bad("0.0.0.");
        bad("0.0.0.256");
        bad("00.0.0.0");
        bad("1.2.3.4.");
        bad("1.2.3.4x");
        bad("1.2.3.300");

        good("0.0.0.0");
        good("1.2.3.4");
        good("1.2.3.42");

        check("0.0.0.0", 0x00000000);
        check("1.2.3.4", 0x01020304);
        check("32.64.128.1", 0x20408001);
        check("255.255.255.255", 0xffffffff);
    }

    void
    run()
    {
        testMembers();
        testParse();
    }
};

TEST_SUITE(
    ipv4_address_test,
    "boost.url.host_ipv4_address");

} // urls
} // boost
