//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_address.hpp>

#include <boost/test/tools/context.hpp>
#include <boost/test/unit_test.hpp>

#include <memory>
#include <string>

#include "test_unit/printing.hpp"

using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_any_address)

// Helper for C++11
static host_and_port make_hport(std::string host, unsigned short port)
{
    host_and_port res;
    res.host = std::move(host);
    res.port = port;
    return res;
}

static std::unique_ptr<any_address> make_address_ptr(any_address from)
{
    return std::unique_ptr<any_address>{new any_address(std::move(from))};
}

BOOST_AUTO_TEST_CASE(default_ctor)
{
    // Default constructed addresses are empty host and ports
    any_address addr;
    BOOST_TEST(addr.type() == address_type::host_and_port);
    BOOST_TEST(addr.hostname() == "");
    BOOST_TEST(addr.port() == (unsigned short)3306);
}

BOOST_AUTO_TEST_CASE(ctor_from_host_and_port)
{
    // Setup
    host_and_port hp = make_hport("abcd", 2000);

    // Construct
    any_address addr(hp);
    hp = make_hport("0000", 1111);  // lifetime check

    // Check
    BOOST_TEST(addr.type() == address_type::host_and_port);
    BOOST_TEST(addr.hostname() == "abcd");
    BOOST_TEST(addr.port() == (unsigned short)2000);
}

BOOST_AUTO_TEST_CASE(ctor_from_unix_path)
{
    // Setup
    unix_path p{"/var/sock"};

    // Construct
    any_address addr(p);
    p.path = "/aaa/bbbb";  // lifetime check

    // Check
    BOOST_TEST(addr.type() == address_type::unix_path);
    BOOST_TEST(addr.unix_socket_path() == "/var/sock");
}

BOOST_AUTO_TEST_CASE(copy_ctor)
{
    // Setup
    auto addr = make_address_ptr(make_hport("abcd", 2000));

    // Construct
    any_address addr2(*addr);
    addr.reset();  // lifetime check

    // Check
    BOOST_TEST(addr2.type() == address_type::host_and_port);
    BOOST_TEST(addr2.hostname() == "abcd");
    BOOST_TEST(addr2.port() == (unsigned short)2000);
}

BOOST_AUTO_TEST_CASE(move_ctor)
{
    // Setup
    auto addr = make_address_ptr(make_hport("abcd", 2000));

    // Construct
    any_address addr2(std::move(*addr));
    addr.reset();  // lifetime check

    // Check
    BOOST_TEST(addr2.type() == address_type::host_and_port);
    BOOST_TEST(addr2.hostname() == "abcd");
    BOOST_TEST(addr2.port() == (unsigned short)2000);
}

BOOST_AUTO_TEST_CASE(copy_assign)
{
    // Setup
    auto addr = make_address_ptr(unix_path{"/var/blah"});
    any_address addr2(make_hport("blah", 9999));

    // Assign
    addr2 = *addr;
    addr.reset();  // lifetime check

    // Check
    BOOST_TEST(addr2.type() == address_type::unix_path);
    BOOST_TEST(addr2.unix_socket_path() == "/var/blah");
}

BOOST_AUTO_TEST_CASE(move_assign)
{
    //  Setup
    auto addr = make_address_ptr(make_hport("abcd", 2000));
    any_address addr2(unix_path{"/var/sock"});

    // Assign
    addr2 = std::move(*addr);
    addr.reset();  // lifetime check

    // Check
    BOOST_TEST(addr2.type() == address_type::host_and_port);
    BOOST_TEST(addr2.hostname() == "abcd");
    BOOST_TEST(addr2.port() == (unsigned short)2000);
}

BOOST_AUTO_TEST_CASE(const_accessors_host_and_port)
{
    // host and port accessors can be called on const objects
    const any_address addr{make_hport("abcd", 2000)};
    BOOST_TEST(addr.type() == address_type::host_and_port);
    BOOST_TEST(addr.hostname() == "abcd");
    BOOST_TEST(addr.port() == (unsigned short)2000);
}

BOOST_AUTO_TEST_CASE(const_accessors_unix_socket)
{
    // UNIX socket accessors can be called on const objects
    const any_address addr{unix_path{"/var/sock"}};
    BOOST_TEST(addr.type() == address_type::unix_path);
    BOOST_TEST(addr.unix_socket_path() == "/var/sock");
}

BOOST_AUTO_TEST_CASE(emplace_host_and_port)
{
    // Changing type
    any_address addr{unix_path{"/var/sock"}};
    addr.emplace_host_and_port("abcd", 2000);
    BOOST_TEST(addr.type() == address_type::host_and_port);
    BOOST_TEST(addr.hostname() == "abcd");
    BOOST_TEST(addr.port() == (unsigned short)2000);

    // Without changing type
    addr.emplace_host_and_port("def", 3000);
    BOOST_TEST(addr.type() == address_type::host_and_port);
    BOOST_TEST(addr.hostname() == "def");
    BOOST_TEST(addr.port() == (unsigned short)3000);

    // Default port value
    addr.emplace_host_and_port("aaa");
    BOOST_TEST(addr.type() == address_type::host_and_port);
    BOOST_TEST(addr.hostname() == "aaa");
    BOOST_TEST(addr.port() == (unsigned short)3306);
}

BOOST_AUTO_TEST_CASE(emplace_unix_path)
{
    // Changing type
    any_address addr{make_hport("abcd", 2000)};
    addr.emplace_unix_path("/var/sock");
    BOOST_TEST(addr.type() == address_type::unix_path);
    BOOST_TEST(addr.unix_socket_path() == "/var/sock");

    // Without changing type
    addr.emplace_unix_path("/var/blah");
    BOOST_TEST(addr.type() == address_type::unix_path);
    BOOST_TEST(addr.unix_socket_path() == "/var/blah");
}

BOOST_AUTO_TEST_CASE(operator_eq_ne)
{
    // For regression check: UNIX socket paths should compare equal
    // whether they were created directly or via emplace
    any_address addr_unix(make_hport("abcd", 3306));
    addr_unix.emplace_unix_path("abcd");

    struct
    {
        const char* name;
        any_address addr1;
        any_address addr2;
        bool equals;
    } test_cases[] = {
        {"host_and_port_eq", make_hport("abc", 2000), make_hport("abc", 2000), true},
        {"host_and_port_eq_default", host_and_port{}, host_and_port{}, true},
        {"host_and_port_ne_host", make_hport("abcd", 2000), make_hport("abc", 2000), false},
        {"host_and_port_ne_port", make_hport("abcd", 2001), make_hport("abcd", 2000), false},
        {"host_and_port_ne_all", make_hport("abc", 2001), make_hport("abcd", 2000), false},

        {"unix_eq", unix_path{"/var/sock"}, unix_path{"/var/sock"}, true},
        {"unix_eq_default", unix_path{}, unix_path{}, true},
        {"unix_ne", unix_path{"/sock1"}, unix_path{"/sock2"}, false},
        {"unix_emplace_regression", addr_unix, unix_path{"abcd"}, true},

        {"type_ne", make_hport("abcd", 0), unix_path{"abcd"}, false},
        {"type_ne_empty", host_and_port{}, unix_path{}, false},
        {"all_ne", make_hport("abcd", 2000), unix_path{"/var/sock"}, false},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            BOOST_TEST((tc.addr1 == tc.addr2) == tc.equals);
            BOOST_TEST((tc.addr2 == tc.addr1) == tc.equals);
            BOOST_TEST((tc.addr1 != tc.addr2) == !tc.equals);
            BOOST_TEST((tc.addr2 != tc.addr1) == !tc.equals);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()