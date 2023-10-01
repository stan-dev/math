//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/url_view.hpp>

#include <boost/url/parse.hpp>
#include <boost/url/url.hpp>
#include <boost/core/ignore_unused.hpp>

#include "test_rule.hpp"

#include <sstream>

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

namespace boost {
namespace urls {

class url_view_test
{
public:
    void
    testSpecialMembers()
    {
        // url_view()
        // ~url_view()
        {
            url_view u;
            BOOST_TEST(u.empty());
            BOOST_TEST_EQ(u.size(), 0u);
        }

        // url_view(url_view const&)
        {
            url_view u1("x://y/z?#");
            url_view u2(u1);
            BOOST_TEST_EQ(u2.data(), u1.data());
            BOOST_TEST_EQ(u2.size(), u1.size());
        }

        // operator=(url_view const&)
        {
            url_view u1("x://y/z?#");
            url_view u2;
            u2 = u1;
            BOOST_TEST_EQ(u2.data(), u1.data());
        }

#if defined(__clang__) && defined(__has_warning)
# if __has_warning("-Wself-assign-overloaded")
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wself-assign-overloaded"
# endif
#endif
        // operator=(url_view const&)
        {
            url_view u1("x://y/z?#");
            u1 = u1;
            BOOST_TEST_EQ(u1.data(), u1.data());
        }

#if defined(__clang__) && defined(__has_warning)
# if __has_warning("-Wself-assign-overloaded")
#  pragma clang diagnostic pop
# endif
#endif

        // operator=(url_view_base const&)
        {
            url u1("x://y/z?#");
            url_view u2 = u1;
            u2 = u1;
            BOOST_TEST_EQ(u1.data(), u2.data());
        }

        // url_view(string_view)
        {
            BOOST_TEST_NO_THROW(url_view(
                "http://example.com/path/to/file.txt?#"));
            BOOST_TEST_THROWS(url_view("{}"),
                std::exception);
        }
    }

    void
    testObservers()
    {
        // size()
        {
            url_view u;
            BOOST_TEST_EQ(u.size(), 0u);
            u = url_view("/");
            BOOST_TEST_EQ(u.size(), 1u);
        }

        // empty()
        {
            url_view u;
            BOOST_TEST(u.empty());
            u = url_view("/");
            BOOST_TEST(! u.empty());
        }

        // data()
        {
            string_view s = "/index.htm";
            url_view u(s);
            BOOST_TEST_NE(u.data(), nullptr);
            BOOST_TEST_EQ(u.data(), s.data());
        }

        // string()
        {
            string_view s = "/index.htm";
            url_view u = parse_relative_ref(s).value();
            BOOST_TEST_EQ(u.buffer(), s);
            BOOST_TEST_EQ(u.buffer().data(), s.data());
        }

        // persist()
        {
        std::shared_ptr<url_view const> sp;
        {
            std::string s( "http://example.com" );
            url_view u( s );                        // u references characters in s

            assert( u.data() == s.data() );         // same buffer

            sp = u.persist();

            assert( sp->data() != s.data() );       // different buffer
            assert( sp->buffer() == s);        // same contents

            // s is destroyed and thus u
            // becomes invalid, but sp remains valid.
        }
        }
    }

    void
    testScheme()
    {
        auto const check = [](
            string_view s,
            char const* m,
            scheme id)
        {
            result<url_view> r =
                parse_uri_reference(s);
            if(! BOOST_TEST(r))
                return;
            url_view u = r.value();
            if(m)
            {
                BOOST_TEST(u.scheme() ==
                           string_view(m));
                BOOST_TEST(
                    u.scheme_id() == id);
            }
            else
            {
                BOOST_TEST(! u.has_scheme());
                BOOST_TEST(u.scheme_id() ==
                           scheme::none);
            }
        };

        auto const bad = [](
            string_view s)
        {
            result<url_view> r =
                parse_uri_reference(s);
            BOOST_TEST(r.has_error());
        };

        check("http://", "http", scheme::http);
        check("ou812://", "ou812", scheme::unknown);
        check("/x", nullptr, scheme::unknown);

        check("http://", "http", scheme::http);
        check("HTTP://", "HTTP", scheme::http);
        check("HtTp://", "HtTp", scheme::http);
        check("a1steak://", "a1steak", scheme::unknown);

        bad("1x:");
        bad(" ");
        bad(" http");
        bad("http ");
    }

    void
    testAuthority()
    {
        auto const no =
            [](string_view s)
        {
            BOOST_TEST_NO_THROW([s]
            {
                url_view u(s);
                BOOST_TEST(! u.has_authority());
            }());
        };
        auto const yes =
            [](string_view s, string_view m)
        {
            //BOOST_TEST_NO_THROW(
            //[&]
            {
                url_view u(s);
                BOOST_TEST(u.has_authority());
                BOOST_TEST_EQ(
                    u.encoded_authority(), m);
                BOOST_TEST_EQ(
                    u.authority().buffer(), m);
            }
            //());
        };

        no("http:xyz/");
        no("http:/x");
        no("http:/x");
        no("http:%2f%2f");
        no("http:/%40");

        yes("http://", "");
        yes("http://a", "a");
        yes("http://a@", "a@");
        yes("http://:@", ":@");
        yes("http://@", "@");
        yes("http://@x", "@x");

        {
            url_view u("http:/path");
            BOOST_TEST_EQ(u.encoded_host(), "");
        }

        // Docs
        assert( url_view( "http://www.example.com/index.htm" ).has_authority() == true );

        assert( url_view( "//" ).has_authority() == true );

        assert( url_view( "/file.txt" ).has_authority() == false );
    }

    void
    testUserinfo()
    {
        auto const no =
            [](string_view s)
        {
            BOOST_TEST_NO_THROW(
            [s]{
                url_view u(s);
                BOOST_TEST(! u.has_userinfo());
            }());
        };
        auto const yes =
            []( string_view s,
                string_view m1,
                string_view m2)
        {
            BOOST_TEST_NO_THROW(
            [&]{
                url_view u(s);
                BOOST_TEST(u.has_userinfo());
                BOOST_TEST_EQ(
                    u.encoded_userinfo(), m1);
                BOOST_TEST_EQ(
                    u.userinfo(), m2);
                BOOST_TEST_EQ(
                    u.authority().encoded_userinfo(), m1);
                BOOST_TEST_EQ(
                    u.authority().userinfo(), m2);
            }());
        };

        no("http:");
        no("http:xyz/");
        no("http:/x");
        no("http:/x");
        no("http:%2f%2f");
        no("http:/%40");
        no("http://");
        no("http://a");

        yes("http://a@", "a", "a");
        yes("http://:@", ":", ":");
        yes("http://@", "", "");
        yes("http://@x", "", "");
        yes("http://%61@x", "%61", "a");
        yes("http://:%61@x", ":%61", ":a");
        yes("http://%61%3a%62@x", "%61%3a%62", "a:b");

        {
            url_view u("x://@");
            BOOST_TEST(u.has_userinfo());
            BOOST_TEST_EQ(u.encoded_userinfo(), "");
            BOOST_TEST_EQ(u.userinfo(), "");
            BOOST_TEST_EQ(u.encoded_user(), "");
            BOOST_TEST_EQ(u.user(), "");
            BOOST_TEST_EQ(u.has_password(), false);
            BOOST_TEST_EQ(u.encoded_password(), "");
            BOOST_TEST_EQ(u.password(), "");
            BOOST_TEST(u.authority().has_userinfo());
            BOOST_TEST_EQ(u.authority().encoded_userinfo(), "");
            BOOST_TEST_EQ(u.authority().userinfo(), "");
            BOOST_TEST_EQ(u.authority().encoded_user(), "");
            BOOST_TEST_EQ(u.authority().user(), "");
            BOOST_TEST_EQ(u.authority().has_password(), false);
            BOOST_TEST_EQ(u.authority().encoded_password(), "");
            BOOST_TEST_EQ(u.authority().password(), "");
        }
        {
            url_view u("x://:@");
            BOOST_TEST(u.has_userinfo());
            BOOST_TEST_EQ(u.encoded_userinfo(), ":");
            BOOST_TEST_EQ(u.userinfo(), ":");
            BOOST_TEST_EQ(u.encoded_user(), "");
            BOOST_TEST_EQ(u.user(), "");
            BOOST_TEST_EQ(u.has_password(), true);
            BOOST_TEST_EQ(u.encoded_password(), "");
            BOOST_TEST_EQ(u.password(), "");
            BOOST_TEST(u.authority().has_userinfo());
            BOOST_TEST_EQ(u.authority().encoded_userinfo(), ":");
            BOOST_TEST_EQ(u.authority().userinfo(), ":");
            BOOST_TEST_EQ(u.authority().encoded_user(), "");
            BOOST_TEST_EQ(u.authority().user(), "");
            BOOST_TEST_EQ(u.authority().has_password(), true);
            BOOST_TEST_EQ(u.authority().encoded_password(), "");
            BOOST_TEST_EQ(u.authority().password(), "");
        }
        {
            url_view u("x://a%41:@");
            BOOST_TEST(u.has_userinfo());
            BOOST_TEST_EQ(u.encoded_userinfo(), "a%41:");
            BOOST_TEST_EQ(u.encoded_user(), "a%41");
            BOOST_TEST_EQ(u.user(), "aA");
            BOOST_TEST_EQ(u.has_password(), true);
            BOOST_TEST_EQ(u.encoded_password(), "");
            BOOST_TEST_EQ(u.password(), "");
            BOOST_TEST(u.authority().has_userinfo());
            BOOST_TEST_EQ(u.authority().encoded_userinfo(), "a%41:");
            BOOST_TEST_EQ(u.authority().encoded_user(), "a%41");
            BOOST_TEST_EQ(u.authority().user(), "aA");
            BOOST_TEST_EQ(u.authority().has_password(), true);
            BOOST_TEST_EQ(u.authority().encoded_password(), "");
            BOOST_TEST_EQ(u.authority().password(), "");
        }
        {
            url_view u("x://:b%42@");
            BOOST_TEST(u.has_userinfo());
            BOOST_TEST_EQ(u.encoded_userinfo(), ":b%42");
            BOOST_TEST_EQ(u.encoded_user(), "");
            BOOST_TEST_EQ(u.user(), "");
            BOOST_TEST_EQ(u.has_password(), true);
            BOOST_TEST_EQ(u.encoded_password(), "b%42");
            BOOST_TEST_EQ(u.password(), "bB");
            BOOST_TEST(u.authority().has_userinfo());
            BOOST_TEST_EQ(u.authority().encoded_userinfo(), ":b%42");
            BOOST_TEST_EQ(u.authority().encoded_user(), "");
            BOOST_TEST_EQ(u.authority().user(), "");
            BOOST_TEST_EQ(u.authority().has_password(), true);
            BOOST_TEST_EQ(u.authority().encoded_password(), "b%42");
            BOOST_TEST_EQ(u.authority().password(), "bB");
        }
        {
            url_view u("x://a:b@");
            BOOST_TEST(u.has_userinfo());
            BOOST_TEST_EQ(u.encoded_userinfo(), "a:b");
            BOOST_TEST_EQ(u.encoded_user(), "a");
            BOOST_TEST_EQ(u.has_password(), true);
            BOOST_TEST_EQ(u.encoded_password(), "b");
            BOOST_TEST(u.authority().has_userinfo());
            BOOST_TEST_EQ(u.authority().encoded_userinfo(), "a:b");
            BOOST_TEST_EQ(u.authority().encoded_user(), "a");
            BOOST_TEST_EQ(u.authority().has_password(), true);
            BOOST_TEST_EQ(u.authority().encoded_password(), "b");
        }
        {
            url_view u("x://%3a:%3a@");
            BOOST_TEST(u.has_userinfo());
            BOOST_TEST_EQ(u.encoded_userinfo(), "%3a:%3a");
            BOOST_TEST_EQ(u.userinfo(), ":::");
            BOOST_TEST_EQ(u.encoded_user(), "%3a");
            BOOST_TEST_EQ(u.user(), ":");
            BOOST_TEST_EQ(u.has_password(), true);
            BOOST_TEST_EQ(u.encoded_password(), "%3a");
            BOOST_TEST_EQ(u.password(), ":");
            BOOST_TEST(u.authority().has_userinfo());
            BOOST_TEST_EQ(u.authority().encoded_userinfo(), "%3a:%3a");
            BOOST_TEST_EQ(u.authority().userinfo(), ":::");
            BOOST_TEST_EQ(u.authority().encoded_user(), "%3a");
            BOOST_TEST_EQ(u.authority().user(), ":");
            BOOST_TEST_EQ(u.authority().has_password(), true);
            BOOST_TEST_EQ(u.authority().encoded_password(), "%3a");
            BOOST_TEST_EQ(u.authority().password(), ":");
        }
        {
            url_view u("x://%2525@");
            BOOST_TEST(u.has_userinfo());
            BOOST_TEST_EQ(u.encoded_userinfo(), "%2525");
            BOOST_TEST_EQ(u.userinfo(), "%25");
            BOOST_TEST_EQ(u.encoded_user(), "%2525");
            BOOST_TEST_EQ(u.user(), "%25");
            BOOST_TEST_EQ(u.has_password(), false);
            BOOST_TEST_EQ(u.encoded_password(), "");
            BOOST_TEST_EQ(u.password(), "");
            BOOST_TEST(u.authority().has_userinfo());
            BOOST_TEST_EQ(u.authority().encoded_userinfo(), "%2525");
            BOOST_TEST_EQ(u.authority().userinfo(), "%25");
            BOOST_TEST_EQ(u.authority().encoded_user(), "%2525");
            BOOST_TEST_EQ(u.authority().user(), "%25");
            BOOST_TEST_EQ(u.authority().has_password(), false);
            BOOST_TEST_EQ(u.authority().encoded_password(), "");
            BOOST_TEST_EQ(u.authority().password(), "");
        }
    }

    //--------------------------------------------
    //
    // Host
    //
    //--------------------------------------------

    void
    testPort()
    {
        {
            url_view u("http://");
            BOOST_TEST(! u.has_port());
            BOOST_TEST_EQ(u.port(), "");
            BOOST_TEST_EQ(u.port_number(), 0);
            BOOST_TEST(! u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "");
            BOOST_TEST_EQ(u.authority().port_number(), 0);
        }
        {
            url_view u("http://www");
            BOOST_TEST(! u.has_port());
            BOOST_TEST_EQ(u.port(), "");
            BOOST_TEST_EQ(u.port_number(), 0);
            BOOST_TEST(! u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "");
            BOOST_TEST_EQ(u.authority().port_number(), 0);
        }
        {
            url_view u("http://:");
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port(), "");
            BOOST_TEST_EQ(u.port_number(), 0);
            BOOST_TEST(u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "");
            BOOST_TEST_EQ(u.authority().port_number(), 0);
        }
        {
            url_view u("http://:0");
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port(), "0");
            BOOST_TEST_EQ(u.port_number(), 0);
            BOOST_TEST(u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "0");
            BOOST_TEST_EQ(u.authority().port_number(), 0);
        }
        {
            url_view u("http://:42");
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port(), "42");
            BOOST_TEST_EQ(u.port_number(), 42);
            BOOST_TEST(u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "42");
            BOOST_TEST_EQ(u.authority().port_number(), 42);
        }
        {
            url_view u("http://:00000");
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port(), "00000");
            BOOST_TEST_EQ(u.port_number(), 0);
            BOOST_TEST(u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "00000");
            BOOST_TEST_EQ(u.authority().port_number(), 0);
        }
        {
            url_view u("http://:000001");
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port(), "000001");
            BOOST_TEST_EQ(u.port_number(), 1);
            BOOST_TEST(u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "000001");
            BOOST_TEST_EQ(u.authority().port_number(), 1);
        }
        {
            url_view u("http://:65535");
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port(), "65535");
            BOOST_TEST_EQ(u.port_number(), 65535);
            BOOST_TEST(u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "65535");
            BOOST_TEST_EQ(u.authority().port_number(), 65535);
        }
        {
            url_view u("http://:65536");
            BOOST_TEST(u.has_port());
            BOOST_TEST_EQ(u.port(), "65536");
            BOOST_TEST_EQ(u.port_number(), 0);
            BOOST_TEST(u.authority().has_port());
            BOOST_TEST_EQ(u.authority().port(), "65536");
            BOOST_TEST_EQ(u.authority().port_number(), 0);
        }
    }

    void
    testHostAndPort()
    {
        {
            url_view u("http://x:1");
            BOOST_TEST(u.encoded_host_and_port() ==
                "x:1");
            BOOST_TEST(u.authority().encoded_host_and_port() ==
                "x:1");
        }
        {
            url_view u("http://x%3a:1");
            BOOST_TEST(u.encoded_host_and_port() ==
                "x%3a:1");
            BOOST_TEST(u.authority().encoded_host_and_port() ==
                "x%3a:1");
        }
        {
            url_view u("http://:1");
            BOOST_TEST(u.encoded_host_and_port() ==
                ":1");
            BOOST_TEST(u.authority().encoded_host_and_port() ==
                ":1");
        }
        {
            url_view u("http://:000001");
            BOOST_TEST(u.encoded_host_and_port() ==
                ":000001");
            BOOST_TEST(u.authority().encoded_host_and_port() ==
                ":000001");
        }
        {
            url_view u("http://xyz:99999");
            BOOST_TEST(u.encoded_host_and_port() ==
                "xyz:99999");
            BOOST_TEST(u.authority().encoded_host_and_port() ==
                "xyz:99999");
        }
    }

    void
    testOrigin()
    {
        BOOST_TEST(url_view(
            "x://p:q@a.b.c/f.z?a=b#frag"
                ).encoded_origin() == "x://p:q@a.b.c");
        BOOST_TEST(url_view(
            "/file.txt").encoded_origin() == "");
        BOOST_TEST(url_view("x:/path/file/txt"
            ).encoded_origin() == "");
    }

    void
    testPath()
    {
        {
            url_view u;
            BOOST_TEST_NO_THROW(u =
                url_view("/path/to/file.htm"));
            BOOST_TEST(u.encoded_path() ==
                "/path/to/file.htm");
            auto const p = u.encoded_segments();
            BOOST_TEST(! p.empty());
            BOOST_TEST_EQ(p.size(), 3u);
            auto it = p.begin();
            BOOST_TEST_EQ(*it, "path");
            ++it;
            BOOST_TEST_EQ(*it, "to");
            ++it;
            BOOST_TEST_EQ(*it, "file.htm");
            ++it;
            BOOST_TEST_EQ(it, p.end());
        }
    }

    void
    testQuery()
    {
        {
            url_view u("http://");
            BOOST_TEST(! u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "");
            BOOST_TEST_EQ(u.query(), "");
        }
        {
            url_view u("http://?");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "");
            BOOST_TEST_EQ(u.query(), "");
        }
        {
            url_view u("http://?k");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "k");
            BOOST_TEST_EQ(u.query(), "k");
        }
        {
            url_view u("http://?k=");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "k=");
            BOOST_TEST_EQ(u.query(), "k=");
        }
        {
            url_view u("http://?#");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "");
            BOOST_TEST_EQ(u.query(), "");
        }
        {
            url_view u("http://?%3f");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "%3f");
            BOOST_TEST_EQ(u.query(), "?");
        }
        {
            url_view u("http://?%25");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "%25");
            BOOST_TEST_EQ(u.query(), "%");
        }
        {
            url_view u("http://?&");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "&");
            BOOST_TEST_EQ(u.query(), "&");
        }
        {
            url_view u("http://?%26");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "%26");
            BOOST_TEST_EQ(u.query(), "&");
        }
        {
            url_view u("http://?a%3db%26");
            BOOST_TEST(u.has_query());
            BOOST_TEST_EQ(u.encoded_query(), "a%3db%26");
            BOOST_TEST_EQ(u.query(), "a=b&");
        }

        // VFALCO TODO test params()
    }

    void
    testFragment()
    {
        auto const check = [](
            string_view s,
            char const* encoded,
            string_view plain)
        {
            result<url_view> r =
                parse_uri_reference(s);
            if(! BOOST_TEST(r))
                return;
            url_view u = r.value();
            if(encoded)
            {
                BOOST_TEST(
                    u.has_fragment());
                BOOST_TEST(
                    u.encoded_fragment() ==
                        string_view(encoded));
                BOOST_TEST_EQ(
                    u.fragment(), plain);
            }
            else
            {
                BOOST_TEST(
                    ! u.has_fragment());
            }
        };

        auto const bad = [](string_view s)
        {
            result<url_view> r =
                parse_uri_reference(s);
            BOOST_TEST(r.has_error());
        };

        check("", nullptr, {});
        check("#", "", "");
        check("/#", "", "");
        check("/#A", "A", "A");
        check("/#%41", "%41", "A");
        check("/?#%41", "%41", "A");
        check(
            "#/?:@!$&'()*+,;=",
            "/?:@!$&'()*+,;=",
            "/?:@!$&'()*+,;=");

        bad("#%%");

        // coverage
        {
            url_view u;
            BOOST_TEST(
                u.encoded_fragment() == "");
            BOOST_TEST_EQ(u.fragment(), "");
        }
    }

    //--------------------------------------------

    void
    testParse()
    {
        // parse_absolute_uri
        {
            result<url_view> r;

            r = parse_absolute_uri(
                "http://user:pass@example.com:443/path/to/file.txt?q");
            BOOST_TEST(r.has_value());
            BOOST_TEST_NO_THROW(r.value());

            r = parse_absolute_uri("");
            BOOST_TEST(r.has_error());
            BOOST_TEST_THROWS(
                r.value(), std::exception);
        }

        // parse_uri
        {
            result<url_view> r;

            r = parse_uri(
                "http://user:pass@example.com:443/path/to/file.txt?q#f");
            BOOST_TEST(r.has_value());
            BOOST_TEST_NO_THROW(r.value());

            r = parse_uri("");
            BOOST_TEST(r.has_error());
            BOOST_TEST_THROWS(
                r.value(), std::exception);
        }

        // parse_relative_ref
        {
            result<url_view> r;

            r = parse_relative_ref(
                "//example.com/path/to/file.txt?q#f");
            BOOST_TEST(r.has_value());
            BOOST_TEST_NO_THROW(r.value());

            r = parse_relative_ref("http:file.txt");
            BOOST_TEST(r.has_error());
            BOOST_TEST_THROWS(
                r.value(), std::exception);
        }

        // parse_uri_reference
        {
            result<url_view> r;

            r = parse_uri_reference(
                "http://user:pass@example.com:443/path/to/file.txt?q#f");
            BOOST_TEST(r.has_value());
            BOOST_TEST_NO_THROW(r.value());

            r = parse_uri_reference(
                "//example.com/path/to/file.txt?q#f");
            BOOST_TEST(r.has_value());

            r = parse_uri_reference("");
            BOOST_TEST(r.has_value());
            BOOST_TEST_NO_THROW(r.value());

            r = parse_uri_reference("1000://");
            BOOST_TEST(r.has_error());
            BOOST_TEST_THROWS(r.value(),
                std::exception);
        }

    }

    void
    testOutput()
    {
        url_view u( "http://example.com" );
        std::stringstream ss;
        ss << u;
        BOOST_TEST(
            ss.str() == "http://example.com");
    }

    void
    testCases()
    {
        BOOST_TEST_NO_THROW(url_view(
            "javascript:alert(1)"));
    }

    void
    testRelativePart()
    {
        auto const ok = [](
            string_view s)
        {
            BOOST_TEST_NO_THROW(
                parse_relative_ref(s).value());
        };

        // "//" authority path-abempty
        {
            ok("//");
            ok("///");
            ok("////");
            ok("///x");
            ok("///:");
            ok("///x/");
            ok("///%3a/");
            ok("///%20");
            ok("///%20");
            ok("///%25");
            ok("///%25%2e");

            ok("//x");
            ok("//x/");
            ok("//x//");
            ok("//x/x");
            ok("//x/:");
            ok("//x/x/");
            ok("//x/%3a/");
            ok("//x/%20");
            ok("//x/%20");
            ok("//x/%25");
            ok("//x/%25%2e");

            ok("");
            ok("/");
            ok("//");
            ok("//user:pass@");
            ok("//boost.org");
            ok("//1.2.3.4:8080");
            ok("//1.2.3.4:8080/");
            ok("//1.2.3.4:8080/x");
            ok("//1.2.3.4:8080/x/");
            ok("//1.2.3.4:8080////");
            ok("/x");
            ok("/x/");
            ok("/x/y");
            ok("/x/y//");
            ok("x");
            ok("x/");
            ok("x//");
            ok("x/y/z");
            ok("x//y///z///");

            //bad(":/"); // colon not ok in relative-part
        }

        // path-absolute
        {
            ok("/");
            ok("/x");
            ok("/x/");
            ok("/:/");
            ok("/x//");
            ok("/%20");
            ok("/:%20");
            ok("/%20");
            ok("/%25");
            ok("/%25%2e");
        }

        // path-noscheme
        {
            ok(".");
            ok("x");
            ok("%20");
            ok("%2f");
            ok("a/");
            ok("a//");
            ok("a/x");
            ok("a/x/");
            ok("a/x//");
            ok("a///");
        }

        // path-abempty
        {
            ok("");
            ok("/");
            ok("//");
            ok("/x");
            ok("/:");
            ok("/x/");
            ok("/%3a/");
            ok("/%20");
            ok("/%20");
            ok("/%25");
            ok("/%25%2e");
        }
    }

    void
    testParseOriginForm()
    {
        BOOST_TEST(parse_origin_form("/").has_value());
        BOOST_TEST(parse_origin_form("/x").has_value());
        BOOST_TEST(parse_origin_form("//").has_value());
        BOOST_TEST(parse_origin_form("/x/").has_value());
        BOOST_TEST(parse_origin_form("/x/y").has_value());
        BOOST_TEST(parse_origin_form("/?").has_value());
        BOOST_TEST(parse_origin_form("/?a").has_value());
        BOOST_TEST(parse_origin_form("/?a=").has_value());
        BOOST_TEST(parse_origin_form("/?a=b").has_value());
        BOOST_TEST(parse_origin_form("/x/y?a=b&c=d").has_value());

        BOOST_TEST(parse_origin_form("").has_error());
        BOOST_TEST(parse_origin_form(" ").has_error());
        BOOST_TEST(parse_origin_form("*").has_error());
        BOOST_TEST(parse_origin_form("?").has_error());
    }

    void
    testJavadocs()
    {
        // {class}
        {
    url_view u( "https://www.example.com/index.htm?text=none#a1" );

    ignore_unused(u);
        }
        {
    result< url_view > rv = parse_uri_reference( "https://www.example.com/index.htm?text=none#a1" );

    ignore_unused(rv);
        }

        // url_view()
        {
        url_view u;

        ignore_unused(u);
        }

        // url_view(string_view)
        {
        url_view u( "http://www.example.com/index.htm" );

        ignore_unused(u);
        }
    }

    void
    run()
    {
        testSpecialMembers();
        testObservers();
        testScheme();
        testAuthority();
        testUserinfo();
        testPort();
        testHostAndPort();
        testOrigin();
        testPath();
        testQuery();
        testFragment();
        testParse();
        testOutput();
        testCases();
        testRelativePart();

        testParseOriginForm();

        testJavadocs();
    }
};

TEST_SUITE(
    url_view_test,
    "boost.url.url_view");

} // urls
} // boost
