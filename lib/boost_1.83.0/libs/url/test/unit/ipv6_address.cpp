//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/ipv6_address.hpp>

#include <boost/url/ipv4_address.hpp>
#include "test_suite.hpp"
#include <sstream>

namespace boost {
namespace urls {

class ipv6_address_test
{
public:
    static
    void
    testMembers()
    {
        // ipv6_address()
        {
            ipv6_address a;
            BOOST_TEST(
                a == ipv6_address());
            ipv6_address::bytes_type b{{}};
            BOOST_TEST_EQ(a.to_bytes(), b);
        }

        // ipv6_address(ipv6_address const&)
        {
            ipv6_address::bytes_type b = { {
                1, 0, 2, 0, 3, 0, 4, 0,
                5, 0, 6, 0, 7, 0, 8, 0 } };
            ipv6_address a1(b);
            ipv6_address a2(a1);
            ipv6_address a3;
            BOOST_TEST_EQ(a1, a2);
            BOOST_TEST_NE(a1, a3);
            BOOST_TEST_NE(a2, a3);
        }

        // operator=(ipv6_address const&)
        {
            ipv6_address::bytes_type b = { {
                1, 0, 2, 0, 3, 0, 4, 0,
                5, 0, 6, 0, 7, 0, 8, 0 } };
            ipv6_address a1(b);
            ipv6_address a2;
            a2 = a1;
            BOOST_TEST_EQ(a2, a1);
        }

        // ipv6_address(bytes_type)
        {
            ipv6_address::bytes_type b = { {
                1, 0, 2, 0, 3, 0, 4, 0,
                5, 0, 6, 0, 7, 0, 8, 0 } };
            ipv6_address a(b);
            BOOST_TEST_EQ(a.to_bytes(), b);
        }

        // ipv6_address(ipv4_address)
        {
            ipv4_address v4("1.2.3.4");
            ipv6_address a(v4);
            BOOST_TEST(a.is_v4_mapped());
            BOOST_TEST(
                a.to_string() == "::ffff:1.2.3.4");
        }

        // ipv6_address(core::string_view)
        {
            ipv6_address a("::");
            BOOST_TEST_EQ(a, ipv6_address());
            BOOST_TEST_THROWS(ipv6_address("x"),
                system::system_error);
        }

        // to_bytes
        {
            ipv6_address::bytes_type b = { {
                1, 0, 2, 0, 3, 0, 4, 0,
                5, 0, 6, 0, 7, 0, 8, 0 } };
            ipv6_address a;
            a = b;
            BOOST_TEST_EQ(a.to_bytes(), b);
        }

        // to_string
        {
            ipv6_address::bytes_type b = { {
                0, 1, 0, 2, 0, 3, 0, 4,
                0, 5, 0, 6, 0, 7, 0, 8 } };
            ipv6_address a(b);
            BOOST_TEST(a.to_string() ==
                "1:2:3:4:5:6:7:8");
        }

        // to_buffer
        {
            ipv6_address::bytes_type b = { {
                0, 1, 0, 2, 0, 3, 0, 4,
                0, 5, 0, 6, 0, 7, 0, 8 } };
            ipv6_address a(b);
            char buf[ipv6_address::max_str_len];
            BOOST_TEST(a.to_buffer(buf, sizeof(buf)) ==
                "1:2:3:4:5:6:7:8");
            char buf2[10];
            BOOST_TEST_THROWS(
                a.to_buffer(buf2, sizeof(buf2)),
                system::system_error);
        }

        // is_unspecified
        {
            BOOST_TEST(ipv6_address().is_unspecified());
            BOOST_TEST(ipv6_address("::").is_unspecified());
            BOOST_TEST(ipv6_address(
                "0:0:0:0:0:0:0:0").is_unspecified());
        }

        // is_loopback
        {
            BOOST_TEST(! ipv6_address().is_loopback());
            BOOST_TEST(ipv6_address("::1").is_loopback());
        }

        // is_v4_mapped
        {
            BOOST_TEST(ipv6_address(
                ipv4_address("1.2.3.4")).is_v4_mapped());
            BOOST_TEST(! ipv6_address().is_v4_mapped());
        }

        // operator==
        // operator!=
        {
            ipv6_address a1("1::");
            ipv6_address a2("2::");
            ipv6_address a3("1::");
            BOOST_TEST_NE(a1, a2);
            BOOST_TEST_EQ(a1, a3);
            BOOST_TEST_NE(a2, a3);
        }

        // static loopback()
        {
            BOOST_TEST(ipv6_address::loopback() ==
                ipv6_address("::1"));
            BOOST_TEST(
                ipv6_address::loopback().is_loopback());
        }

        // operator<<
        {
            std::stringstream ss;
            ss << ipv6_address("1:0:0:0:0:0:0:1");
            BOOST_TEST_EQ(ss.str(), "1::1");
        }
    }

    static
    void
    bad(core::string_view s)
    {
        BOOST_TEST(parse_ipv6_address(
            s).has_error());
    }

    static
    void
    good(core::string_view s)
    {
        BOOST_TEST_EQ(ipv6_address(
            s).to_string(), s);
    }

    static
    void
    trip(core::string_view s0, core::string_view s1)
    {
        auto ip0 = ipv6_address(s0);
        BOOST_TEST_EQ(ip0.to_string(), s1);
        // round-trip
        auto ip1 = ipv6_address(
            ip0.to_string());
        BOOST_TEST_EQ(ip1, ip0);
    }

    static
    std::uint64_t
    get_u64(std::uint8_t const* p)
    {
        return
            (static_cast<std::uint64_t>(p[0]) << 56) +
            (static_cast<std::uint64_t>(p[1]) << 48) +
            (static_cast<std::uint64_t>(p[2]) << 40) +
            (static_cast<std::uint64_t>(p[3]) << 32) +
            (static_cast<std::uint64_t>(p[4]) << 24) +
            (static_cast<std::uint64_t>(p[5]) << 16) +
            (static_cast<std::uint64_t>(p[6]) <<  8) +
             static_cast<std::uint64_t>(p[7]);
    }

    static
    void
    check(
        core::string_view s,
        std::uint64_t u0,
        std::uint64_t u1)
    {
        ipv6_address a;
        auto r = parse_ipv6_address(s);
        if(! BOOST_TEST(! r.has_error()))
            return;
        a = r.value();
        auto const bytes = a.to_bytes();
        BOOST_TEST(
            get_u64(&bytes[0]) == u0);
        BOOST_TEST(
            get_u64(&bytes[8]) == u1);
        return;
    }

    //--------------------------------------------

    void
    testIO()
    {
        bad("");
        bad(":");
        bad("0");
        bad("0:");
        bad(":0");
        bad(":::");
        bad("x::");
        bad(":0::");
        bad("0:12");
        bad("0:123");
        bad("::0::");
        bad("0::0:x");
        bad("0:1.2.3.4");
        bad("::FFFF:999.2.3.4");
        bad("0:0:0:0:0:0:0:1.2.3.4");
        bad("0:0:0:0:0:0:0::1.2.3.4");

        bad("::1.");
        bad("::1.2");
        bad("::1.2");
        bad("::1.2x");
        bad("::1.2.");
        bad("::1.2.3");
        bad("::1.2.3");
        bad("::1.2.3x");
        bad("::1.2.3.");
        bad("::1.2.3.4x");

        // coverage
        bad("::ffff:260.168.0.1");
        bad("::ffff:999.168.0.1");
        bad("::ffff:1a0.168.0.1");
        bad("::ffff:a.168.0.1");
        bad("::ffff:0a.168.0.1");
        bad("::ffff:0000a.168.0.1");
        bad("::ffff:192.168.0.");
        bad("0:0:0:0:0:0:0:");
        bad("0");
        bad("0:");
        bad("0:0:0:0:ffff");

        good("1::");
        good("12::");
        good("123::");
        good("1234::");
        good("abcd::");
        good("::ffff:1.2.3.4");
        good("1234:1234:1234:1234:1234:1234:1234:1234");

        trip("ABCD::", "abcd::");
        trip("0:0:0:0:0:0:0:0", "::");
        trip("1:0:0:0:0:0:0:0", "1::");
        trip("0:1:0:0:0:0:0:0", "0:1::");
        trip("0:0:1:0:0:0:0:0", "0:0:1::");
        trip("0:0:0:1:0:0:0:0", "0:0:0:1::");
        trip("0:0:0:0:1:0:0:0", "::1:0:0:0");
        trip("0:0:0:0:0:1:0:0", "::1:0:0");
        trip("0:0:0:0:0:0:1:0", "::1:0");
        trip("0:0:0:0:0:0:0:1", "::1");
        trip("1234:1234:1234:1234:1234:1234:255.255.255.255",
              "1234:1234:1234:1234:1234:1234:ffff:ffff");
        trip("0:0:0:0:0:ffff:1.2.3.4", "::ffff:1.2.3.4");

        check("1:2:3:4:5:6:7:8", 0x0001000200030004, 0x0005000600070008);
        check("::2:3:4:5:6:7:8", 0x0000000200030004, 0x0005000600070008);
        check("1::3:4:5:6:7:8",  0x0001000000030004, 0x0005000600070008);
        check("1:2::4:5:6:7:8",  0x0001000200000004, 0x0005000600070008);
        check("1:2:3::5:6:7:8",  0x0001000200030000, 0x0005000600070008);
        check("1:2:3:4::6:7:8",  0x0001000200030004, 0x0000000600070008);
        check("1:2:3:4:5::7:8",  0x0001000200030004, 0x0005000000070008);
        check("1:2:3:4:5:6::8",  0x0001000200030004, 0x0005000600000008);
        check("1:2:3:4:5:6:7::", 0x0001000200030004, 0x0005000600070000);
        check("::3:4:5:6:7:8", 0x0000000000030004, 0x0005000600070008);
        check("1::4:5:6:7:8",  0x0001000000000004, 0x0005000600070008);
        check("1:2::5:6:7:8",  0x0001000200000000, 0x0005000600070008);
        check("1:2:3::6:7:8",  0x0001000200030000, 0x0000000600070008);
        check("1:2:3:4::7:8",  0x0001000200030004, 0x0000000000070008);
        check("1:2:3:4:5::8",  0x0001000200030004, 0x0005000000000008);
        check("1:2:3:4:5:6::", 0x0001000200030004, 0x0005000600000000);
        check("::4:5:6:7:8", 0x0000000000000004, 0x0005000600070008);
        check("1::5:6:7:8",  0x0001000000000000, 0x0005000600070008);
        check("1:2::6:7:8",  0x0001000200000000, 0x0000000600070008);
        check("1:2:3::7:8",  0x0001000200030000, 0x0000000000070008);
        check("1:2:3:4::8",  0x0001000200030004, 0x0000000000000008);
        check("1:2:3:4:5::", 0x0001000200030004, 0x0005000000000000);
        check("::5:6:7:8", 0x0000000000000000, 0x0005000600070008);
        check("1::6:7:8",  0x0001000000000000, 0x0000000600070008);
        check("1:2::7:8",  0x0001000200000000, 0x0000000000070008);
        check("1:2:3::8",  0x0001000200030000, 0x0000000000000008);
        check("1:2:3:4::", 0x0001000200030004, 0x0000000000000000);
        check("::6:7:8", 0x0000000000000000, 0x0000000600070008);
        check("1::7:8",  0x0001000000000000, 0x0000000000070008);
        check("1:2::8",  0x0001000200000000, 0x0000000000000008);
        check("1:2:3::", 0x0001000200030000, 0x0000000000000000);
        check("::7:8", 0x0000000000000000, 0x0000000000070008);
        check("1::8",  0x0001000000000000, 0x0000000000000008);
        check("1:2::", 0x0001000200000000, 0x0000000000000000);
        check("::8", 0x0000000000000000, 0x0000000000000008);
        check("1::", 0x0001000000000000, 0x0000000000000000);

        check("::0", 0, 0);
        check("::1", 0, 1);
        check("0:0:0::1", 0, 1);
        check("0:0:0:0:0:0:0:0", 0, 0);
        check("0:0:0:0:0:0:0.0.0.0", 0, 0);
        check("::1.2.3.4", 0, 0x01020304);
        check("::1234:5678", 0, 0x0000000012345678);
        check("::FFFF:1.2.3.4", 0, 0x0000ffff01020304);
        check("1:2::3:4:5", 0x0001000200000000, 0x0000000300040005);
        check("::1:2:3:4:5", 0x0000000000000001, 0x0002000300040005);
        check("1:2:3:4:5::", 0x0001000200030004, 0x0005000000000000);
        check("1:2:0:0:0:3:4:5", 0x0001000200000000, 0x0000000300040005);
        check("1:2:3:4:5:0:0:0", 0x0001000200030004, 0x0005000000000000);
        check("0:0:0:1:2:3:4:5", 0x0000000000000001, 0x0002000300040005);
        check("0:0:0:0:0:FFFF:102:405", 0, 0x0000ffff01020405);
        check("0000:0000:0000:0000:0000:0000:0000:0000", 0, 0);
        check("1234:5678:9ABC:DEF0:0000:0000:0000:0000", 0x123456789abcdef0, 0);
        check("2001:0DB8:0A0B:12F0:0000:0000:0000:0001", 0x20010db80a0b12f0, 0x0000000000000001);
        check("2001:DB8:3333:4444:5555:6666:7777:8888", 0x20010db833334444, 0x5555666677778888);
        check("2001:DB8:3333:4444:CCCC:DDDD:EEEE:FFFF", 0x20010db833334444, 0xccccddddeeeeffff);
        check("2001:db8::", 0x20010db800000000, 0);
        check("2001:DB8::", 0x20010db800000000, 0);
        check("2001:DB8:A0B:12F0::1", 0x20010db80a0b12f0, 1);
        check("2001:db8::1234:5678", 0x20010db800000000, 0x12345678);
        check("2001:DB8::1234:5678", 0x20010db800000000, 0x0000000012345678);
        check("2001:DB8:1::AB9:C0A8:102", 0x20010db800010000, 0x00000ab9c0a80102);
        check("2001:db8:1::ab9:C0A8:102", 0x20010db800010000, 0x00000ab9c0a80102);
        check("2001:0DB8:0A0B:12F0:0:0:0:1", 0x20010db80a0b12f0, 1);
        check("2001:db8:3333:4444:5555:6666:7777:8888", 0x20010db833334444, 0x5555666677778888);
        check("2001:db8:3333:4444:CCCC:DDDD:EEEE:FFFF", 0x20010db833334444, 0xccccddddeeeeffff);
        check("2001:0DB8:0001:0000:0000:0AB9:C0A8:0102", 0x20010db800010000, 0x00000ab9c0a80102);
        check("2001:0db8:0a0b:12f0:0000:0000:0000:0001", 0x20010db80a0b12f0, 1);
        check("2001:0db8:0001:0000:0000:0ab9:C0A8:0102", 0x20010db800010000, 0x0ab9c0a80102);
        check("3FFE:1900:4545:3:200:F8FF:FE21:67CF", 0x3ffe190045450003, 0x0200f8fffe2167cf);
        check("684D:1111:222:3333:4444:5555:6:77", 0x684d111102223333, 0x4444555500060077);
        check("fe80:0:0:0:200:f8ff:fe21:67cf", 0xfe80000000000000, 0x0200f8fffe2167cf);
        check("FE80:0:0:0:200:F8FF:FE21:67CF", 0xfe80000000000000, 0x0200f8fffe2167cf);
        check("FFFF:0:0:0:0:0:0:1", 0xffff000000000000, 1);
        check("FFFF::1", 0xffff000000000000, 1);
    }

    void
    testIpv4()
    {
        // ipv4_mapped
        BOOST_TEST(ipv6_address(
            "::ffff:1.2.3.4").is_v4_mapped());
        BOOST_TEST(! ipv6_address(
            "1::ffff:1.2.3.4").is_v4_mapped());
        BOOST_TEST(ipv6_address("::ffff:c0a8:0001"
            ).to_string() == "::ffff:192.168.0.1");
        BOOST_TEST(ipv6_address("::ffff:192.168.0.1"
            ).to_string() == "::ffff:192.168.0.1");
        BOOST_TEST(ipv6_address("1::ffff:c0a8:0001"
            ).to_string() != "::ffff:192.168.0.1");
        BOOST_TEST(ipv6_address("1::ffff:192.168.0.1"
            ).to_string() != "::ffff:192.168.0.1");

        ipv4_address a(0x7f000001);
        BOOST_TEST(ipv6_address(a).to_string() ==
            "::ffff:127.0.0.1");

        char buf[ipv6_address::max_str_len];
        BOOST_TEST(ipv6_address(a
            ).to_buffer(buf, sizeof(buf)) ==
                "::ffff:127.0.0.1");
    }

    void
    run()
    {
        testMembers();
        testIO();
        testIpv4();
    }
};

TEST_SUITE(
    ipv6_address_test,
    "boost.url.host_ipv6_address");

} // urls
} // boost
