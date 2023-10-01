//
// Copyright (c) 2022 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/url_base.hpp>

#include <boost/url/decode_view.hpp>
#include <boost/url/url.hpp>
#include "test_suite.hpp"

/*  Legend

    '#' 0x23    '=' 0x3d
    '%' 0x25    '@' 0x40
    '&' 0x26    '[' 0x5b
    '.' 0x2e    ']' 0x5d
    ':' 0x3a
*/

namespace boost {
namespace urls {

struct url_base_test
{
    template<class F>
    static
    void
    modify(
        string_view before,
        string_view after,
        F&& f)
    {
        url u(before);
        f(u);
        auto s = u.buffer();
        BOOST_TEST_EQ(s, after);
    }

    //--------------------------------------------
    //
    // Scheme
    //
    //--------------------------------------------

    void
    testSetScheme()
    {
        auto const remove = [](
            string_view s1, string_view s2)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST_EQ(u.remove_scheme().buffer(), s2);
            BOOST_TEST_EQ(u.scheme_id(), scheme::none);
            BOOST_TEST(u.scheme().empty());
        };

        auto const set = [](
            scheme id, string_view s1,
            string_view s2, string_view s3)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s2));
            BOOST_TEST(u.set_scheme(s1).buffer() == s3);
            BOOST_TEST_EQ(u.scheme(), s1);
            BOOST_TEST_EQ(u.scheme_id(), id);
        };

        auto const setid = [](
            scheme id, string_view s1, string_view s2)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.set_scheme_id(id).buffer() == s2);
            BOOST_TEST_EQ(u.scheme_id(), id);
        };

        remove("",      "");
        remove("x",     "x");
        remove("x:",    "");
        remove("x:/",   "/");
        remove("x:a",   "a");
        remove("x:a/",  "a/");
        remove("x:/a", "/a");
        remove("x://a", "//a");
        remove("x:///a", "///a");
        remove("x://",  "//");
        remove("x:a:",  "a%3A");
        remove("x:a:/", "a%3A/");
        remove("yabba:dabba:doo", "dabba%3Adoo");
        remove("yabba:dabba:doo:doo", "dabba%3Adoo%3Adoo");
        remove("x::::", "%3A%3A%3A");


        remove("x://a.b/1/2",      "//a.b/1/2");
        remove("x://a:b@c.d/1/?#", "//a:b@c.d/1/?#");

        set(scheme::ftp,  "ftp",  "",     "ftp:");
        set(scheme::ws,   "ws",   "/",    "ws:/");
        set(scheme::ws,   "ws",   "a",    "ws:a");
        set(scheme::ws,   "ws",   "a/",   "ws:a/");
        set(scheme::ws,   "ws",   "//",   "ws://");
        set(scheme::ws,   "ws",   "a:/",  "ws:/");
        set(scheme::http, "http", "./a:", "http:a:");

        set(scheme::ws, "ws", "//a.b/1/2",      "ws://a.b/1/2");
        set(scheme::ws, "ws", "//a:b@c.d/1/?#", "ws://a:b@c.d/1/?#");

        setid(scheme::ftp, "",    "ftp:");
        setid(scheme::ws,  "/",   "ws:/");
        setid(scheme::ws,  "a",   "ws:a");
        setid(scheme::ws,  "a/",  "ws:a/");
        setid(scheme::ws,  "//",  "ws://");
        setid(scheme::ws,  "a:/", "ws:/");

        setid(scheme::ws,
            "//a.b/1/2", "ws://a.b/1/2");

        setid(scheme::ws,
            "//a:b@c.d/1/?#",  "ws://a:b@c.d/1/?#");

        setid(scheme::none, "a:/", "/");

        BOOST_TEST_THROWS(
            url().set_scheme(""),
            system_error);

        BOOST_TEST_THROWS(
            url().set_scheme("http~"),
            system_error);

        BOOST_TEST_THROWS(
            url().set_scheme_id(scheme::unknown),
            system_error);

        // self-intersection
        modify(
            "x://?mailto",
            "mailto://?mailto",
            [](url_base& u)
            {
                u.set_scheme(
                    u.encoded_query());
            });
    }

    //--------------------------------------------
    //
    // Authority
    //
    //--------------------------------------------

    void
    testSetAuthority()
    {
        auto const remove = [](
            string_view s1, string_view s2)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST_EQ(u.remove_authority().buffer(), s2);
            BOOST_TEST(u.encoded_authority().empty());
            BOOST_TEST(! u.has_authority());
            BOOST_TEST(! u.has_userinfo());
            BOOST_TEST(! u.has_password());
            BOOST_TEST(! u.has_port());
            BOOST_TEST(u.host_address().empty());
        };

        auto const set = [](string_view s1,
            string_view s2, string_view s3)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.set_encoded_authority(s2).buffer() == s3);
            BOOST_TEST_EQ(u.encoded_authority(), s2);
            BOOST_TEST(u.has_authority());
        };

        BOOST_TEST_THROWS(
            url().set_encoded_authority("x:y"),
            system_error);

        BOOST_TEST_THROWS(
            url().set_encoded_authority("%2"),
            system_error);

        remove("", "");
        remove("/", "/");
        remove("/x", "/x");
        remove("/x/", "/x/");
        remove("/x/y", "/x/y");
        remove("x/", "x/");
        remove("x/y", "x/y");
        remove("x/y/", "x/y/");
        remove("x/y/?#", "x/y/?#");
        remove("x/y/?#", "x/y/?#");



        remove("z:", "z:");
        remove("z:/", "z:/");
        remove("z:/x", "z:/x");
        remove("z:/x/", "z:/x/");
        remove("z:/x/y", "z:/x/y");
        remove("z:x/", "z:x/");
        remove("z:x/y", "z:x/y");
        remove("z:x/y/", "z:x/y/");
        remove("z:x/y/?#", "z:x/y/?#");
        remove("z:x:/y/?#", "z:x:/y/?#");

        remove("//", "");
        remove("///", "/");
        remove("///x", "/x");
        remove("///x/", "/x/");
        remove("///x/y", "/x/y");
        remove("//x/", "/");
        remove("//x/y", "/y");
        remove("//x/y/", "/y/");
        remove("//x/y/?#", "/y/?#");

        remove("z://", "z:");
        remove("z:///", "z:/");
        remove("z:///x", "z:/x");
        remove("z:///x/", "z:/x/");
        remove("z:///x/y", "z:/x/y");
        remove("z://x/", "z:/");
        remove("z://x/y", "z:/y");
        remove("z://x/y/", "z:/y/");
        remove("z://x/y/?#", "z:/y/?#");
        remove("z://x:/y/?#", "z:/y/?#");
        remove("z://x//y/?q#f", "z:/.//y/?q#f");

        set("", "", "//");
        set("", "x@", "//x@");
        set("", ":x@", "//:x@");
        set("", "x:y@", "//x:y@");
        set("", "x", "//x");
        set("", "x.y", "//x.y");
        set("", "x:", "//x:");
        set("", ":", "//:");
        set("", ":0", "//:0");
        set("", ":443", "//:443");
        set("", ":65536", "//:65536");
        set("", "1.2.3.4", "//1.2.3.4");
        set("", "[v1.0]", "//[v1.0]");
        set("", "[::]", "//[::]");
        set("", "[::ffff:127.0.0.1]",
                "//[::ffff:127.0.0.1]");
        set("", "[::ffff:127.0.0.1]:80",
                "//[::ffff:127.0.0.1]:80");
        set("", "user:pass@example.com:80",
                "//user:pass@example.com:80");
        set("ws:",
                "user:pass@example.com:80",
                "ws://user:pass@example.com:80");

        set("///a", "", "///a");
        set("///a", "x@", "//x@/a");
        set("///a", ":x@", "//:x@/a");
        set("///a", "x:y@", "//x:y@/a");
        set("///a", "x", "//x/a");
        set("///a", "x.y", "//x.y/a");
        set("///a", "x:", "//x:/a");
        set("///a", ":", "//:/a");
        set("///a", ":0", "//:0/a");
        set("///a", ":443", "//:443/a");
        set("///a", ":65536", "//:65536/a");
        set("///a", "1.2.3.4", "//1.2.3.4/a");
        set("///a", "[v1.0]", "//[v1.0]/a");
        set("///a", "[::]", "//[::]/a");
        set("///a", "[::ffff:127.0.0.1]",
                    "//[::ffff:127.0.0.1]/a");
        set("///a", "[::ffff:127.0.0.1]:80",
                    "//[::ffff:127.0.0.1]:80/a");
        set("///a", "user:pass@example.com:80",
                    "//user:pass@example.com:80/a");
        set("ws:///a",
                    "user:pass@example.com:80",
                    "ws://user:pass@example.com:80/a");

        // self-intersection
        modify(
            "x://@?user:pass@example.com:8080",
            "x://user:pass@example.com:8080?user:pass@example.com:8080",
            [](url_base& u)
            {
                u.set_encoded_authority(
                    u.encoded_query());
            });
    }

    //--------------------------------------------
    //
    // Userinfo
    //
    //--------------------------------------------

    void
    testSetUserinfo()
    {
        auto const remove = [](
            string_view s1, string_view s2)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST_EQ(u.remove_userinfo().buffer(), s2);
            BOOST_TEST(u.encoded_userinfo().empty());
            BOOST_TEST(u.userinfo().empty());
            BOOST_TEST(! u.has_userinfo());
        };

        auto const set = [](string_view s1,
            string_view s2, string_view s3)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST_EQ(u.set_userinfo(s2).buffer(), s3);
        };

        auto const enc = [](string_view s1,
            string_view s2, string_view s3)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.set_encoded_userinfo(s2).buffer() == s3);
            BOOST_TEST_EQ(u.encoded_userinfo(), s2);
            BOOST_TEST(u.has_userinfo());
        };

        BOOST_TEST_THROWS(url().set_encoded_userinfo("%2"), std::exception);

        remove("", "");
        remove("/", "/");
        remove("//", "//");
        remove("//@", "//");
        remove("//a@", "//");
        remove("//a:@", "//");
        remove("//a:b@", "//");
        remove("//@x", "//x");
        remove("//a@x", "//x");
        remove("//a:b@x", "//x");
        remove("//a:b@x/", "//x/");

        remove("z:", "z:");
        remove("z:/", "z:/");
        remove("z://", "z://");
        remove("z://@", "z://");
        remove("z://a@", "z://");
        remove("z://a:@", "z://");
        remove("z://a:b@", "z://");
        remove("z://@x", "z://x");
        remove("z://a@x", "z://x");
        remove("z://a:b@x", "z://x");

        set("", "", "//@");
        set("/", "", "//@/");
        set("//", "", "//@");
        set("//@", "", "//@");
        set("//a@", "", "//@");
        set("//a:@", "", "//@");
        set("//a:b@", "", "//@");
        set("//@x", "", "//@x");
        set("//a@x", "", "//@x");
        set("//a:b@x", "", "//@x");
        set("//a:b@x/", "", "//@x/");

        set("w:", "", "w://@");
        set("w:/", "", "w://@/");
        set("w://", "", "w://@");
        set("w://@", "", "w://@");
        set("w://a@", "", "w://@");
        set("w://a:@", "", "w://@");
        set("w://a:b@", "", "w://@");
        set("w://@x", "", "w://@x");
        set("w://a@x", "", "w://@x");
        set("w://a:b@x", "", "w://@x");
        set("w://a:b@x/", "", "w://@x/");

        set("", ":", "//:@");
        set("/", "a", "//a@/");
        set("//", "@", "//%40@");
        set("//@", "xyz", "//xyz@");
        set("//a@", ":@", "//:%40@");
        set("//a:@", "x", "//x@");
        set("//a:b@", "p:q", "//p:q@");
        set("//@x", "z", "//z@x");
        set("//a@x", "42", "//42@x");
        set("//a:b@x", "UV", "//UV@x");
        set("//a:b@x/", "NR", "//NR@x/");

        set("w:", ":", "w://:@");
        set("w:/", "a", "w://a@/");
        set("w://", "@", "w://%40@");
        set("w://@", "xyz", "w://xyz@");
        set("w://a@", ":@", "w://:%40@");
        set("w://a:@", "x", "w://x@");
        set("w://a:b@", "p:q", "w://p:q@");
        set("w://@x", "z", "w://z@x");
        set("w://a@x", "42", "w://42@x");
        set("w://a:b@x", "UV", "w://UV@x");
        set("w://a:b@x/", "NR", "w://NR@x/");

        enc("", "", "//@");
        enc("/", "", "//@/");
        enc("//", "", "//@");
        enc("//@", "", "//@");
        enc("//a@", "", "//@");
        enc("//a:@", "", "//@");
        enc("//a:b@", "", "//@");
        enc("//@x", "", "//@x");
        enc("//a@x", "", "//@x");
        enc("//a:b@x", "", "//@x");
        enc("//a:b@x/", "", "//@x/");

        enc("w:", "", "w://@");
        enc("w:/", "", "w://@/");
        enc("w://", "", "w://@");
        enc("w://@", "", "w://@");
        enc("w://a@", "", "w://@");
        enc("w://a:@", "", "w://@");
        enc("w://a:b@", "", "w://@");
        enc("w://@x", "", "w://@x");
        enc("w://a@x", "", "w://@x");
        enc("w://a:b@x", "", "w://@x");
        enc("w://a:b@x/", "", "w://@x/");

        enc("", ":", "//:@");
        enc("", "%3a", "//%3a@");
        enc("/", "%41", "//%41@/");
        enc("//", "x", "//x@");
        enc("//@", "xyz", "//xyz@");
        enc("//a@", "%3a%40", "//%3a%40@");
        enc("//a:@", "x", "//x@");
        enc("//a:b@", "p:q", "//p:q@");
        enc("//@x", "z", "//z@x");
        enc("//a@x", "42", "//42@x");
        enc("//a:b@x", "UV", "//UV@x");
        enc("//a:b@x/", "NR", "//NR@x/");

        enc("w:", ":", "w://:@");
        enc("w:", "%3a", "w://%3a@");
        enc("w:/", "%41", "w://%41@/");
        enc("w://", "x", "w://x@");
        enc("w://@", "xyz", "w://xyz@");
        enc("w://a@", "%3a%40", "w://%3a%40@");
        enc("w://a:@", "x", "w://x@");
        enc("w://a:b@", "p:q", "w://p:q@");
        enc("w://@x", "z", "w://z@x");
        enc("w://a@x", "42", "w://42@x");
        enc("w://a:b@x", "UV", "w://UV@x");
        enc("w://a:b@x/", "NR", "w://NR@x/");

        // self-intersection
        modify(
            "x://?user:pass",
            "x://user:pass@?user:pass",
            [](url_base& u)
            {
                u.set_encoded_userinfo(u.encoded_query());
            });
        modify(
            "x://?user:pass",
            "x://user:pass@?user:pass",
            [](url_base& u)
            {
                u.set_userinfo(
                    u.encoded_query());
            });
        modify(
            "x://?user:pass",
            "x://user:pass@?user:pass",
            [](url_base& u)
            {
                u.set_userinfo(
                    u.query());
            });
    }

    void
    testSetUser()
    {
        auto const set = [](
            string_view s0,
            string_view s,
            string_view s1)
        {
            modify(s0, s1,
                [s](url_base& u)
                {
                    u.set_user(s);
                    BOOST_TEST(u.user() == s);
                    BOOST_TEST(u.has_userinfo());
                });
        };

        auto const enc = [](
            string_view s0,
            string_view s,
            string_view s1)
        {
            modify(s0, s1,
                [s](url_base& u)
                {
                    BOOST_TEST_NO_THROW(u.set_encoded_user(s));
                    BOOST_TEST_EQ(
                        decode_view(s),
                        decode_view(u.encoded_user()));
                    BOOST_TEST(u.has_userinfo());
                });
        };

        BOOST_TEST_THROWS(
            url().set_encoded_user("%2"),
            system_error);

        set("", "", "//@");
        set("/y", "", "//@/y");
        set("//", "", "//@");
        set("//y", "", "//@y");
        set("//@", "", "//@");
        set("//:@", "", "//:@");
        set("//y@", "", "//@");
        set("//y@z", "", "//@z");
        set("//y:@", "", "//:@");
        set("//y:z@", "", "//:z@");
        set("//a:b@c", "", "//:b@c");
        set("ws:", "", "ws://@");
        set("ws:/y", "", "ws://@/y");
        set("ws://", "", "ws://@");
        set("ws://y", "", "ws://@y");
        set("ws://@", "", "ws://@");
        set("ws://:@", "", "ws://:@");
        set("ws://y@", "", "ws://@");
        set("ws://y@z", "", "ws://@z");
        set("ws://y:@", "", "ws://:@");
        set("ws://y:z@", "", "ws://:z@");
        set("ws://a:b@c", "", "ws://:b@c");
        set("", "", "//@");
        set("", "x", "//x@");
        set("/y", "x", "//x@/y");
        set("//", "x", "//x@");
        set("//y", "x", "//x@y");
        set("//@", "x", "//x@");
        set("//:@", "x", "//x:@");
        set("//y@", "x", "//x@");
        set("//y@z", "x", "//x@z");
        set("//y:@", "x", "//x:@");
        set("//y:z@", "x", "//x:z@");
        set("//a:b@c", "x", "//x:b@c");
        set("ws:", "x", "ws://x@");
        set("ws:/y", "x", "ws://x@/y");
        set("ws://", "x", "ws://x@");
        set("ws://y", "x", "ws://x@y");
        set("ws://@", "x", "ws://x@");
        set("ws://:@", "x", "ws://x:@");
        set("ws://y@", "x", "ws://x@");
        set("ws://y@z", "x", "ws://x@z");
        set("ws://y:@", "x", "ws://x:@");
        set("ws://y:z@", "x", "ws://x:z@");
        set("ws://a:b@c", "x", "ws://x:b@c");
        set("ws://a:b@c", ":", "ws://%3A:b@c");
        set("ws://a:b@c", "@", "ws://%40:b@c");

        enc("", "", "//@");
        enc("", "%41", "//%41@");
        enc("/y", "%41", "//%41@/y");
        enc("//", "%41", "//%41@");
        enc("//y", "%41", "//%41@y");
        enc("//@", "%41", "//%41@");
        enc("//:@", "%41", "//%41:@");
        enc("//y@", "%41", "//%41@");
        enc("//y@z", "%41", "//%41@z");
        enc("//y:@", "%41", "//%41:@");
        enc("//y:z@", "%41", "//%41:z@");
        enc("//a:b@c", "%41", "//%41:b@c");
        enc("ws:", "%41", "ws://%41@");
        enc("ws:/y", "%41", "ws://%41@/y");
        enc("ws://", "%41", "ws://%41@");
        enc("ws://y", "%41", "ws://%41@y");
        enc("ws://@", "%41", "ws://%41@");
        enc("ws://:@", "%41", "ws://%41:@");
        enc("ws://y@", "%41", "ws://%41@");
        enc("ws://y@z", "%41", "ws://%41@z");
        enc("ws://y:@", "%41", "ws://%41:@");
        enc("ws://y:z@", "%41", "ws://%41:z@");
        enc("ws://a:b@c", "%41", "ws://%41:b@c");
        enc("x:", "user%3apass", "x://user%3apass@");
        enc("x:", "user@local", "x://user%40local@");

        // self-intersection
        modify(
            "x://u@/?johndoe",
            "x://johndoe@/?johndoe",
            [](url_base& u)
            {
                u.set_encoded_user(
                    u.encoded_query());
            });
        modify(
            "x://u@/?johndoe",
            "x://johndoe@/?johndoe",
            [](url_base& u)
            {
                u.set_user(
                    u.query());
            });
        modify(
            "x://u@/?johndoe",
            "x://johndoe@/?johndoe",
            [](url_base& u)
            {
                u.set_user(
                    u.encoded_query());
            });
    }

    void
    testSetPassword()
    {
        auto const remove = [](
            string_view s1, string_view s2)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.remove_password().buffer() == s2);
            BOOST_TEST_EQ(u.encoded_password(), "");
            BOOST_TEST_EQ(u.password(), "");
        };

        auto const set = [](
            string_view s1, string_view s2,
                string_view s3)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.set_password(s2).buffer() == s3);
        };

        auto const enc = [](
            string_view s1, string_view s2,
                string_view s3)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.set_encoded_password(
                s2).buffer() == s3);
            BOOST_TEST_EQ(u.encoded_password(), s2);
            BOOST_TEST(u.has_userinfo());
        };

        BOOST_TEST_THROWS(url().set_encoded_password(
            "%2"), system_error);

        remove("", "");
        remove("/", "/");
        remove("//", "//");
        remove("//", "//");
        remove("//@", "//@");
        remove("//y@", "//y@");
        remove("//:y@", "//@");
        remove("//y:z@", "//y@");
        remove("//y:z@a", "//y@a");

        remove("x:", "x:");
        remove("x:/", "x:/");
        remove("x://", "x://");
        remove("x://", "x://");
        remove("x://@", "x://@");
        remove("x://y@", "x://y@");
        remove("x://:y@", "x://@");
        remove("x://y:z@", "x://y@");
        remove("x://y:z@a", "x://y@a");

        set("", "", "//:@");
        set("/", "", "//:@/");
        set("//", "", "//:@");
        set("//@", "", "//:@");
        set("//y@", "", "//y:@");
        set("//:y@", "", "//:@");
        set("//y:z@", "", "//y:@");
        set("//y:z@a", "", "//y:@a");

        set("x:", "", "x://:@");
        set("x:/", "", "x://:@/");
        set("x://", "", "x://:@");
        set("x://@", "", "x://:@");
        set("x://y@", "", "x://y:@");
        set("x://:y@", "", "x://:@");
        set("x://y:z@", "", "x://y:@");
        set("x://y:z@a", "", "x://y:@a");

        set("", "x", "//:x@");
        set("/", "x", "//:x@/");
        set("//", "x", "//:x@");
        set("//x", "y", "//:y@x");
        set("//x@", "y", "//x:y@");
        set("//x:y@", "z", "//x:z@");
        set("//x:abc@", "z", "//x:z@");
        set("//x:z@", "abc", "//x:abc@");

        set("w:", "x", "w://:x@");
        set("w:/", "x", "w://:x@/");
        set("w://", "x", "w://:x@");
        set("w://x", "y", "w://:y@x");
        set("w://x@", "y", "w://x:y@");
        set("w://x:y@", "z", "w://x:z@");
        set("w://x:abc@", "z", "w://x:z@");
        set("w://x:z@", "abc", "w://x:abc@");

        set("w://x:z@", ":", "w://x::@");
        set("w://x:z@", "@", "w://x:%40@");

        enc("", "", "//:@");
        enc("", "%41", "//:%41@");
        enc("/y", "%41", "//:%41@/y");
        enc("//", "%41", "//:%41@");
        enc("//y", "%41", "//:%41@y");
        enc("//@", "%41", "//:%41@");
        enc("//:@", "%41", "//:%41@");
        enc("//y@", "%41", "//y:%41@");
        enc("//y@z", "%41", "//y:%41@z");
        enc("//y:@", "%41", "//y:%41@");
        enc("//y:z@", "%41", "//y:%41@");
        enc("//a:b@c", "%41", "//a:%41@c");

        enc("ws:", "%41", "ws://:%41@");
        enc("ws:/y", "%41", "ws://:%41@/y");
        enc("ws://", "%41", "ws://:%41@");
        enc("ws://y", "%41", "ws://:%41@y");
        enc("ws://@", "%41", "ws://:%41@");
        enc("ws://:@", "%41", "ws://:%41@");
        enc("ws://y@", "%41", "ws://y:%41@");
        enc("ws://y@z", "%41", "ws://y:%41@z");
        enc("ws://y:@", "%41", "ws://y:%41@");
        enc("ws://y:z@", "%41", "ws://y:%41@");
        enc("ws://a:b@c", "%41", "ws://a:%41@c");

        // self-intersection
        modify(
            "x://:p@/?johndoe",
            "x://:johndoe@/?johndoe",
            [](url_base& u)
            {
                u.set_encoded_password(
                    u.encoded_query());
            });
        modify(
            "x://:p@/?johndoe",
            "x://:johndoe@/?johndoe",
            [](url_base& u)
            {
                u.set_password(
                    u.query());
            });
        modify(
            "x://:p@/?johndoe",
            "x://:johndoe@/?johndoe",
            [](url_base& u)
            {
                u.set_password(
                    u.encoded_query());
            });
    }

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
    testSetHost()
    {
        auto const set_host = [](
            string_view s,
            string_view s1,
            host_type ht)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_host(s));
            BOOST_TEST_EQ(u.host_type(), ht);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            switch(ht)
            {
            case host_type::none:
                BOOST_TEST_FAIL();
                break;
            case host_type::ipv4:
                BOOST_TEST_EQ(u.encoded_host(), s);
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address(s));
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::ipv6:
                BOOST_TEST_EQ(u.encoded_host(), s);
                BOOST_TEST_EQ(bracketed(u.host_address()), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address(u.host_address()));
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::ipvfuture:
                BOOST_TEST_EQ(u.encoded_host(), s);
                BOOST_TEST_EQ(bracketed(u.host_address()), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), u.host_address());
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::name:
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), s);
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            }
        };

        auto const set_encoded_host = [](
            string_view s,
            string_view s1,
            host_type ht)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_encoded_host(s));
            BOOST_TEST_EQ(u.host_type(), ht);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            switch(ht)
            {
            case host_type::none:
                BOOST_TEST_FAIL();
                break;
            case host_type::ipv4:
                BOOST_TEST_EQ(u.encoded_host(), s);
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address(s));
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::ipv6:
                BOOST_TEST_EQ(u.encoded_host(), s);
                BOOST_TEST_EQ(bracketed(u.host_address()), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address(u.host_address()));
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::ipvfuture:
                BOOST_TEST_EQ(u.encoded_host(), s);
                BOOST_TEST_EQ(bracketed(u.host_address()), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), u.host_address());
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::name:
                BOOST_TEST_EQ(u.host_address(), pct_string_view(s).decode());
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), pct_string_view(s).decode());
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            }
        };

        auto const set_host_address = [](
            string_view s,
            string_view s1,
            host_type ht)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_host_address(s));
            BOOST_TEST_EQ(u.host_type(), ht);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            switch(ht)
            {
            case host_type::none:
                BOOST_TEST_FAIL();
                break;
            case host_type::ipv4:
                BOOST_TEST_EQ(u.encoded_host(), s);
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address(s));
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::ipv6:
                BOOST_TEST_EQ(u.encoded_host(), bracketed(s));
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address(s));
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::ipvfuture:
                BOOST_TEST_EQ(u.encoded_host(), bracketed(s));
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), s);
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::name:
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), s);
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            }
        };

        auto const set_encoded_host_address = [](
            string_view s,
            string_view s1,
            host_type ht)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_encoded_host_address(s));
            BOOST_TEST_EQ(u.host_type(), ht);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            switch(ht)
            {
            case host_type::none:
                BOOST_TEST_FAIL();
                break;
            case host_type::ipv4:
                BOOST_TEST_EQ(u.encoded_host(), s);
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address(s));
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::ipv6:
                BOOST_TEST_EQ(u.encoded_host(), bracketed(s));
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address(s));
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::ipvfuture:
                BOOST_TEST_EQ(u.encoded_host(), bracketed(s));
                BOOST_TEST_EQ(u.host_address(), s);
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), s);
                BOOST_TEST_EQ(u.host_name(), "");
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            case host_type::name:
                BOOST_TEST_EQ(u.host_address(), pct_string_view(s).decode());
                BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
                BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
                BOOST_TEST_EQ(u.host_ipvfuture(), "");
                BOOST_TEST_EQ(u.host_name(), pct_string_view(s).decode());
                BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
                break;
            }
        };

        auto const set_host_ipv4 = [](
            string_view s,
            string_view s1)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_host_ipv4(ipv4_address(s)));
            BOOST_TEST_EQ(u.host_type(), host_type::ipv4);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            BOOST_TEST_EQ(u.encoded_host(), s);
            BOOST_TEST_EQ(u.host_address(), s);
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address(s));
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
            BOOST_TEST_EQ(u.host_ipvfuture(), "");
            BOOST_TEST_EQ(u.host_name(), "");
            BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
        };

        auto const set_host_ipv6 = [](
            string_view s,
            string_view s1)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_host_ipv6(ipv6_address(s)));
            BOOST_TEST_EQ(u.host_type(), host_type::ipv6);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            BOOST_TEST_EQ(u.encoded_host(), bracketed(s));
            BOOST_TEST_EQ(u.host_address(), s);
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address(s));
            BOOST_TEST_EQ(u.host_ipvfuture(), "");
            BOOST_TEST_EQ(u.host_name(), "");
            BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
        };

        auto const set_host_ipvfuture = [](
            string_view s,
            string_view s1)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_host_ipvfuture(s))
            BOOST_TEST_EQ(u.host_type(), host_type::ipvfuture);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            BOOST_TEST_EQ(u.encoded_host(), bracketed(s));
            BOOST_TEST_EQ(u.host_address(), s);
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
            BOOST_TEST_EQ(u.host_ipvfuture(), s);
            BOOST_TEST_EQ(u.host_name(), "");
            BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
        };

        auto const set_host_name = [](
            string_view s,
            string_view s1)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_host_name(s))
            BOOST_TEST_EQ(u.host_type(), host_type::name);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), s);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            BOOST_TEST_EQ(u.host_address(), s);
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
            BOOST_TEST_EQ(u.host_ipvfuture(), "");
            BOOST_TEST_EQ(u.host_name(), s);
            BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
        };

        auto const set_encoded_host_name = [](
            string_view s,
            string_view s1)
        {
            url u;
            BOOST_TEST_NO_THROW(u.set_encoded_host_name(s))
            BOOST_TEST_EQ(u.host_type(), host_type::name);
            BOOST_TEST_EQ(u.buffer(), s1);
            BOOST_TEST_EQ(u.host(), u.encoded_host().decode());
            auto rv = parse_ipv4_address(s);
            if(! rv)
                BOOST_TEST_EQ(u.encoded_host(), s);
            else
                BOOST_TEST_EQ(u.host(), rv->to_string());
            BOOST_TEST_EQ(u.host_address(), pct_string_view(s).decode());
            BOOST_TEST_EQ(u.host_ipv4_address(), ipv4_address());
            BOOST_TEST_EQ(u.host_ipv6_address(), ipv6_address());
            BOOST_TEST_EQ(u.host_ipvfuture(), "");
            BOOST_TEST_EQ(u.host_name(), pct_string_view(s).decode());
            BOOST_TEST_EQ(u.host_name(), u.encoded_host_name().decode());
        };

        set_host("", "//", host_type::name);
        set_host("127.0.0.1", "//127.0.0.1", host_type::ipv4);
        set_host("[1::6:c0a8:1]", "//[1::6:c0a8:1]", host_type::ipv6 );
        set_host("[v42.69]", "//[v42.69]", host_type::ipvfuture );
        set_host("www.example.com", "//www.example.com", host_type::name);
        set_host("%5b%3a", "//%255b%253a", host_type::name);

        set_encoded_host("", "//", host_type::name);
        set_encoded_host("127.0.0.1", "//127.0.0.1", host_type::ipv4);
        set_encoded_host("10.2.201.1", "//10.2.201.1", host_type::ipv4);
        set_encoded_host("0.5.15.20", "//0.5.15.20", host_type::ipv4);
        set_encoded_host("100.101.110.115", "//100.101.110.115", host_type::ipv4);
        set_encoded_host("200.205.210.255", "//200.205.210.255", host_type::ipv4);
        set_encoded_host("[1::6:c0a8:1]", "//[1::6:c0a8:1]", host_type::ipv6 );
        set_encoded_host("[::ffff:192.168.102.0]", "//[::ffff:192.168.102.0]", host_type::ipv6 );
        set_encoded_host("[v42.69]", "//[v42.69]", host_type::ipvfuture );
        set_encoded_host("www.example.com", "//www.example.com", host_type::name);
        set_encoded_host("%5b%3a", "//%5b%3a", host_type::name);
        set_encoded_host("%00", "//%00", host_type::name);

        set_host_address("", "//", host_type::name);
        set_host_address("127.0.0.1", "//127.0.0.1", host_type::ipv4);
        set_host_address("1::6:c0a8:1", "//[1::6:c0a8:1]", host_type::ipv6 );
        set_host_address("v42.69", "//[v42.69]", host_type::ipvfuture );
        set_host_address("www.example.com", "//www.example.com", host_type::name);
        set_host_address("%5b%3a", "//%255b%253a", host_type::name);

        set_encoded_host_address("", "//", host_type::name);
        set_encoded_host_address("127.0.0.1", "//127.0.0.1", host_type::ipv4);
        set_encoded_host_address("127%2e0.0.1", "//127%2e0.0.1", host_type::name);
        set_encoded_host_address("1::6:c0a8:1", "//[1::6:c0a8:1]", host_type::ipv6 );
        set_encoded_host_address("v42.69", "//[v42.69]", host_type::ipvfuture );
        set_encoded_host_address("www.example.com", "//www.example.com", host_type::name);
        set_encoded_host_address("%5b%3a", "//%5b%3a", host_type::name);

        set_host_ipv4("0.0.0.0", "//0.0.0.0");
        set_host_ipv4("127.0.0.1", "//127.0.0.1");
        set_host_ipv4("255.255.255.255", "//255.255.255.255");

        set_host_ipv6("1::6:c0a8:1", "//[1::6:c0a8:1]");

        set_host_ipvfuture("v42.69", "//[v42.69]");
        BOOST_TEST_THROWS(url().set_host_ipvfuture("127.0.0.1"), system_error);

        set_host_name("www.example.com", "//www.example.com");
        set_host_name("%5b%3a", "//%255b%253a");
        set_host_name("127.0.0.1", "//127%2E0%2E0%2E1");

        set_encoded_host_name("www.example.com", "//www.example.com");
        set_encoded_host_name("%5b%3a", "//%5b%3a");
        set_encoded_host_name("127.0.0.1", "//127%2E0%2E0%2E1");
        BOOST_TEST_THROWS(url().set_encoded_host_name("%go"), system_error);

        // self-intersection
        modify(
            "x://@?www.example.com",
            "x://@www.example.com?www.example.com",
            [](url_base& u)
            {
                u.set_encoded_host(
                    u.encoded_query());
            });
        modify(
            "x://@?www.example.com",
            "x://@www.example.com?www.example.com",
            [](url_base& u)
            {
                u.set_host(
                    u.encoded_query());
            });
        modify(
            "x://@?www.example.com",
            "x://@www.example.com?www.example.com",
            [](url_base& u)
            {
                u.set_host(
                    u.query());
            });
    }

    void
    testSetPort()
    {
        auto const remove = [](
            string_view s1, string_view s2)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.remove_port().buffer() == s2);
            BOOST_TEST(! u.has_port());
            BOOST_TEST(u.port().empty());
            BOOST_TEST_EQ(u.port_number(), 0);
        };

        auto const setn = [](string_view s1,
            std::uint16_t n, string_view s2)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.set_port_number(n).buffer() == s2);
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port_number(), n);
        };

        auto const set = [](string_view s1,
            std::uint16_t n, string_view s2,
                string_view s3)
        {
            url u;
            BOOST_TEST_NO_THROW(u = url(s1));
            BOOST_TEST(u.set_port(s2).buffer() == s3);
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port_number(), n);
            BOOST_TEST_EQ(u.port(), s2);
        };

        BOOST_TEST_THROWS(
            url().set_port("x"),
            std::exception);

        remove("", "");
        remove("/", "/");
        remove("//", "//");
        remove("//:0", "//");
        remove("//:80", "//");
        remove("//:65535", "//");
        remove("//:999999", "//");
        remove("//:999999/", "///");
        remove("//x:999999/", "//x/");
        remove("//a:b@x.y:8080/path/to/file.txt?#",
               "//a:b@x.y/path/to/file.txt?#");

        remove("x:", "x:");
        remove("x:/", "x:/");
        remove("x://", "x://");
        remove("x://:0", "x://");
        remove("x://:80", "x://");
        remove("x://:65535", "x://");
        remove("x://:999999", "x://");
        remove("x://:999999/", "x:///");
        remove("x://x:999999/", "x://x/");
        remove("x://a:b@x.y:8080/path/to/file.txt?#",
               "x://a:b@x.y/path/to/file.txt?#");

        setn("", 0, "//:0");
        setn("", 443, "//:443");
        setn("", 65535, "//:65535");
        setn("/", 0, "//:0/");
        setn("//", 0, "//:0");
        setn("///", 0, "//:0/");
        setn("//x/", 0, "//x:0/");
        setn("//x/y", 0, "//x:0/y");
        setn("//a:b@/y", 0, "//a:b@:0/y");
        setn("//a:b@c/y", 0, "//a:b@c:0/y");
        setn("//a:b@x.y/path/to/file.txt?#", 8080,
             "//a:b@x.y:8080/path/to/file.txt?#");

        setn("g:", 0, "g://:0");
        setn("g:", 443, "g://:443");
        setn("g:", 65535, "g://:65535");
        setn("g:/", 0, "g://:0/");
        setn("g://", 0, "g://:0");
        setn("g:///", 0, "g://:0/");
        setn("g://x/", 0, "g://x:0/");
        setn("g://x/y", 0, "g://x:0/y");
        setn("g://a:b@/y", 0, "g://a:b@:0/y");
        setn("g://a:b@c/y", 0, "g://a:b@c:0/y");
        setn("g://a:b@x.y/path/to/file.txt?#", 8080,
            "g://a:b@x.y:8080/path/to/file.txt?#");

        set("", 0, "", "//:");
        set("/", 0, "", "//:/");
        set("//", 0, "", "//:");
        set("///", 0, "", "//:/");
        set("//x/", 0, "", "//x:/");
        set("//x/y", 0, "", "//x:/y");
        set("//a:b@/y", 0, "", "//a:b@:/y");
        set("//a:b@c/y", 0, "", "//a:b@c:/y");
        set("//a:b@x.y/path/to/file.txt?#", 0, "",
            "//a:b@x.y:/path/to/file.txt?#");

        set("g:", 0, "", "g://:");
        set("g:/", 0, "", "g://:/");
        set("g://", 0, "", "g://:");
        set("g:///", 0, "", "g://:/");
        set("g://x/", 0, "", "g://x:/");
        set("g://x/y", 0, "", "g://x:/y");
        set("g://a:b@/y", 0, "", "g://a:b@:/y");
        set("g://a:b@c/y", 0, "", "g://a:b@c:/y");
        set("g://a:b@x.y/path/to/file.txt?#", 0, "",
            "g://a:b@x.y:/path/to/file.txt?#");

        set("", 0, "0", "//:0");
        set("", 443, "443", "//:443");
        set("", 65535, "65535", "//:65535");
        set("/", 0, "0", "//:0/");
        set("//", 0, "0", "//:0");
        set("///", 0, "0", "//:0/");
        set("//x/", 0, "0", "//x:0/");
        set("//x/y", 0, "0", "//x:0/y");
        set("//a:b@/y", 0, "0", "//a:b@:0/y");
        set("//a:b@c/y", 0, "0", "//a:b@c:0/y");
        set("//a:b@x.y/path/to/file.txt?#", 8080, "8080",
            "//a:b@x.y:8080/path/to/file.txt?#");

        set("g:", 0, "0", "g://:0");
        set("g:", 443, "443", "g://:443");
        set("g:", 65535, "65535", "g://:65535");
        set("g:/", 0, "0", "g://:0/");
        set("g://", 0, "0", "g://:0");
        set("g:///", 0, "0", "g://:0/");
        set("g://x/", 0, "0", "g://x:0/");
        set("g://x/y", 0, "0", "g://x:0/y");
        set("g://a:b@/y", 0, "0", "g://a:b@:0/y");
        set("g://a:b@c/y", 0, "0", "g://a:b@c:0/y");
        set("g://a:b@x.y/path/to/file.txt?#", 8080, "8080",
            "g://a:b@x.y:8080/path/to/file.txt?#");

        // self-intersection
        modify(
            "x://@?65535",
            "x://@:65535?65535",
            [](url_base& u)
            {
                u.set_port(u.encoded_query());
            });
    }

    //--------------------------------------------

    void
    testQuery()
    {
        // has_query
        {
            {
                url u;
                BOOST_TEST(! u.has_query());
            }
            {
                url u("?");
                BOOST_TEST(u.has_query());
            }
            {
                url u("?x");
                BOOST_TEST(u.has_query());
            }
        }

        // set_encoded_query
        {
            {
                url u;
                u.set_encoded_query("");
                BOOST_TEST(u.has_query());
                BOOST_TEST_EQ(u.buffer(), "?");
                BOOST_TEST_EQ(u.encoded_query(), "");
            }
            {
                url u;
                u.set_encoded_query("x");
                BOOST_TEST(u.has_query());
                BOOST_TEST_EQ(u.buffer(), "?x");
                BOOST_TEST_EQ(u.encoded_query(), "x");
            }
            {
                url u;
                u.set_encoded_query("%41");
                BOOST_TEST(u.has_query());
                BOOST_TEST_EQ(u.buffer(), "?%41");
                BOOST_TEST_EQ(u.encoded_query(), "%41");
                BOOST_TEST_EQ(u.query(), "A");
            }
            {
                url u;
                BOOST_TEST_THROWS(
                    u.set_encoded_query("%%"),
                    system_error);
                BOOST_TEST_THROWS(
                    u.set_encoded_query("%fg"),
                    system_error);
            }
        }

        // set_query
        {
            auto good = [](
                string_view q, string_view us)
            {
                url u;
                u.set_query(q);
                BOOST_TEST(u.has_query());
                BOOST_TEST_EQ(u.buffer(), us);
                BOOST_TEST_EQ(u.query(), q);
            };
            good("", "?");
            good("x", "?x");
            good("%41", "?%2541");
            good("%%fg", "?%25%25fg");
            good("{}", "?%7B%7D");

            // issue #245
            {
                url u;
                u.set_query("");
                u.set_query("");
                BOOST_TEST_EQ(u.buffer(), "?");
            }
        }

        // has_query
        {
            url u;
            BOOST_TEST_NO_THROW(u = url("?query"));
            BOOST_TEST(u.has_query());
            u.clear();
            BOOST_TEST(! u.has_query());
            BOOST_TEST_NO_THROW(u = url("?"));
            BOOST_TEST(u.has_query());
        }

        // remove_query
        {
            url u;
            BOOST_TEST_NO_THROW(u = url("?query"));
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "query");
            BOOST_TEST_EQ(u.params().size(), 1u);
            BOOST_TEST_EQ(u.remove_query().has_query(), false);
            BOOST_TEST_EQ(u.encoded_query(), "");
            BOOST_TEST_EQ(u.query(), "");
            BOOST_TEST_EQ(u.params().size(), 0u);
            BOOST_TEST_EQ(u.encoded_params().size(), 0u);
        }

        // set_encoded_query
        {
            url u;
            BOOST_TEST(! u.has_query());
            u.set_encoded_query("k1=v1&k2=#v2");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.params().size(), 2u);
            BOOST_TEST_EQ((*u.params().begin()).key, "k1");
            BOOST_TEST_EQ((*u.params().begin()).value, "v1");
            BOOST_TEST_EQ((*std::next(u.params().begin())).key, "k2");
            BOOST_TEST_EQ((*std::next(u.params().begin())).value, "#v2");

            u.set_encoded_query("");
            BOOST_TEST(u.has_query());
            BOOST_TEST(u.encoded_query().empty());
            BOOST_TEST_EQ(u.params().size(), 1);
        }

        // set_query
        {
            url u;
            BOOST_TEST(! u.has_query());
            u.set_query("!@#$%^&*()_+=-;:'{}[]|\\?/>.<,");
            BOOST_TEST(u.has_query());
            BOOST_TEST(u.encoded_query() ==
                "!@%23$%25%5E&*()_+=-;:'%7B%7D%5B%5D%7C%5C?/%3E.%3C,");
            BOOST_TEST_EQ(u.params().size(), 2u);
            BOOST_TEST_EQ((*u.params().begin()).key, "!@#$%^");
            BOOST_TEST_EQ((*u.params().begin()).value, "");
            BOOST_TEST_EQ((*std::next(u.params().begin())).key, "*()_ ");
            BOOST_TEST_EQ((*std::next(u.params().begin())).value,
                "-;:'{}[]|\\?/>.<,");
        }

        // remove_query
        {
            {
                url u;
                u.remove_query();
                BOOST_TEST(! u.has_query());
            }
            {
                url u("?");
                u.remove_query();
                BOOST_TEST(! u.has_query());
            }
            {
                url u("?x");
                u.remove_query();
                BOOST_TEST(! u.has_query());
            }
        }

        // self-intersection
        modify(
            "#abracadabra",
            "?abracadabra#abracadabra",
            [](url_base& u)
            {
                u.set_encoded_query(
                    u.encoded_fragment());
            });
        modify(
            "#abracadabra",
            "?abracadabra#abracadabra",
            [](url_base& u)
            {
                u.set_query(
                    u.encoded_fragment());
            });
    }

    void
    testJavadocs()
    {
        //----------------------------------------
        //
        // Scheme
        //
        //----------------------------------------

        // set_scheme
        assert( url( "http://www.example.com" ).set_scheme( "https" ).scheme_id() == scheme::https );

        // set_scheme_id
        assert( url( "http://example.com/echo.cgi" ).set_scheme_id( scheme::wss ).buffer() == "wss://example.com/echo.cgi" );

        // remove_scheme
        assert( url("http://www.example.com/index.htm" ).remove_scheme().buffer() == "//www.example.com/index.htm" );

        //----------------------------------------
        //
        // Authority
        //
        //----------------------------------------

        // set_encoded_authority
        assert( url().set_encoded_authority( "My%20Computer" ).has_authority() );

        // remove_authority
        assert( url( "http://example.com/echo.cgi" ).remove_authority().buffer() == "http:/echo.cgi" );

        //----------------------------------------
        //
        // Userinfo
        //
        //----------------------------------------

        // set_userinfo
        assert( url( "http://example.com" ).set_userinfo( "user:pass" ).encoded_user() == "user" );

        // set_encoded_userinfo
        assert( url( "http://example.com" ).set_encoded_userinfo( "john%20doe" ).user() == "john doe" );

        // remove_userinfo
        assert( url( "http://user@example.com" ).remove_userinfo().has_userinfo() == false );

        //----------------------------------------

        // set_user
        assert( url().set_user("john doe").encoded_userinfo() == "john%20doe" );

        // set_encoded_user
        assert( url().set_encoded_user("john%20doe").userinfo() == "john doe" );

        // set_password
        assert( url("http://user@example.com").set_password( "pass" ).encoded_userinfo() == "user:pass" );

        // set_encoded_password
        assert( url("http://user@example.com").set_encoded_password( "pass" ).encoded_userinfo() == "user:pass" );

        // remove_password
        assert( url( "http://user:pass@example.com" ).remove_password().encoded_authority() == "user@example.com" );

        //----------------------------------------
        //
        // Host
        //
        //----------------------------------------

        // set_host
        assert( url( "http://www.example.com" ).set_host( "127.0.0.1" ).buffer() == "http://127.0.0.1" );

        // set_encoded_host
        assert( url( "http://www.example.com" ).set_host( "127.0.0.1" ).buffer() == "http://127.0.0.1" );

        // set_host_address
        assert( url( "http://www.example.com" ).set_host_address( "127.0.0.1" ).buffer() == "http://127.0.0.1" );

        // set_encoded_host_address
        assert( url( "http://www.example.com" ).set_host( "127.0.0.1" ).buffer() == "http://127.0.0.1" );

        // set_host_ipv4
        assert( url("http://www.example.com").set_host_ipv4( ipv4_address( "127.0.0.1" ) ).buffer() == "http://127.0.0.1" );

        // set_host_ipv6
        assert( url().set_host_ipv6( ipv6_address( "1::6:c0a8:1" ) ).encoded_authority() == "[1::6:c0a8:1]" );

        // set_host_ipvfuture
        assert( url().set_host_ipvfuture( "v42.bis" ).buffer() == "//[v42.bis]" );

        // set_host_name
        assert( url( "http://www.example.com/index.htm").set_host_name( "localhost" ).host_address() == "localhost" );

        // set_encoded_host_name
        assert( url( "http://www.example.com/index.htm").set_encoded_host_name( "localhost" ).host_address() == "localhost" );

        //----------------------------------------
        //
        // Port
        //
        //----------------------------------------

        // set_port
        assert( url( "http://www.example.com" ).set_port_number( 8080 ).encoded_authority() == "www.example.com:8080" );

        // set_port
        assert( url( "http://www.example.com" ).set_port( "8080" ).encoded_authority() == "www.example.com:8080" );

        // remove_port
        assert( url( "http://www.example.com:80" ).remove_port().encoded_authority() == "www.example.com" );

        //----------------------------------------
        //
        // Path
        //
        //----------------------------------------

        // set_path_absolute
        {
        url u( "path/to/file.txt" );
        assert( u.set_path_absolute( true ) );
        assert( u.buffer() == "/path/to/file.txt" );
        }

        // set_path
        {
        url u( "http://www.example.com" );

        u.set_path( "path/to/file.txt" );

        assert( u.path() == "/path/to/file.txt" );
        }

        // set_encoded_path
        {
        url u( "http://www.example.com" );

        u.set_encoded_path( "path/to/file.txt" );

        assert( u.encoded_path() == "/path/to/file.txt" );
        }

        // segments
        {
        url u( "http://example.com/path/to/file.txt" );

        segments_ref sv = u.segments();

        (void)sv;
        }

        // encoded_segments
        {
        url u( "http://example.com/path/to/file.txt" );

        segments_encoded_ref sv = u.encoded_segments();

        (void)sv;
        }

        //----------------------------------------
        //
        // Query
        //
        //----------------------------------------

        // set_query
        assert( url( "http://example.com" ).set_query( "id=42" ).query() == "id=42" );

        // set_encoded_query
        assert( url( "http://example.com" ).set_encoded_query( "id=42" ).encoded_query() == "id=42" );

        // params
        {
        params_ref pv = url( "/sql?id=42&name=jane%2Ddoe&page+size=20" ).params();

        (void)pv;
        }

        // encoded_params
        {
        params_encoded_ref pv = url( "/sql?id=42&name=jane%2Ddoe&page+size=20" ).encoded_params();

        (void)pv;
        }

        // remove_query
        assert( url( "http://www.example.com?id=42" ).remove_query().buffer() == "http://www.example.com" );

        //----------------------------------------
        //
        // Fragment
        //
        //----------------------------------------

        // remove_fragment
        assert( url( "?first=john&last=doe#anchor" ).remove_fragment().buffer() == "?first=john&last=doe" );

        // set_fragment
        assert( url("?first=john&last=doe" ).set_encoded_fragment( "john doe" ).encoded_fragment() == "john%20doe" );

        // set_encoded_fragment
        assert( url("?first=john&last=doe" ).set_encoded_fragment( "john%2Ddoe" ).fragment() == "john-doe" );

        //----------------------------------------
        //
        // Compound Fields
        //
        //----------------------------------------

        // remove_origin
        assert( url( "http://www.example.com/index.htm" ).remove_origin().buffer() == "/index.htm" );
    }

    void
    run()
    {
        testSetScheme();
        testSetAuthority();
        testSetUserinfo();
        testSetUser();
        testSetPassword();
        testSetHost();
        testSetPort();
        testQuery();
        testJavadocs();
    }
};

TEST_SUITE(
    url_base_test,
    "boost.url.url_base");

} // urls
} // boost
