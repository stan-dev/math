//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/format.hpp>

#include <boost/url/parse.hpp>
#include <boost/url/static_url.hpp>

#include "test_suite.hpp"

#ifdef BOOST_TEST_CSTR_EQ
#undef BOOST_TEST_CSTR_EQ
#define BOOST_TEST_CSTR_EQ(expr1,expr2) \
    BOOST_TEST_EQ( boost::urls::detail::to_sv(expr1), boost::urls::detail::to_sv(expr2) )
#endif

namespace boost {
namespace urls {

struct X
{};

namespace detail {
template <>
struct formatter<X>;
}

struct format_test
{
    void
    testFormat()
    {
        {
            url u = urls::format("http:");
            BOOST_TEST_CSTR_EQ( u.buffer(), "http:" );
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
        }
        {
            url u = urls::format("{}:", "http");
            BOOST_TEST_CSTR_EQ( u.buffer(), "http:" );
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST_NOT( u.has_authority() );
        }
        {
            url u = urls::format("{}:", "http");
            BOOST_TEST_CSTR_EQ( u.buffer(), "http:" );
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST_NOT( u.has_authority() );
        }
        {
            url u = urls::format("{}://", "http");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.encoded_host().empty() );
        }
        auto segs_equal = [](
            segments_encoded_view segs,
            std::initializer_list<core::string_view> il)
        {
            if (BOOST_TEST_EQ(segs.size(), il.size()))
            {
                auto it0 = segs.begin();
                auto it1 = il.begin();
                auto end0 = segs.end();
                auto end1 = il.end();
                while (
                    it0 != end0 &&
                    it1 != end1)
                {
                    BOOST_TEST_CSTR_EQ(*it0++, *it1++);
                }
            }
        };

        auto params_equal = [](
            params_encoded_view params,
            std::initializer_list<std::pair<core::string_view, core::string_view>> il)
        {
            if (BOOST_TEST_EQ(params.size(), il.size()))
            {
                auto it0 = params.begin();
                auto it1 = il.begin();
                auto end0 = params.end();
                auto end1 = il.end();
                while (
                    it0 != end0 &&
                    it1 != end1)
                {
                    auto p = *it0;
                    BOOST_TEST_CSTR_EQ(p.key, it1->first);
                    BOOST_TEST_CSTR_EQ(p.value, it1->second);
                    ++it0;
                    ++it1;
                }
            }
        };

        {
            url u = urls::format("{}:///", "http");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http:///");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/" );
            BOOST_TEST( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {""} );

        }
        {
            url u = urls::format("{}://{}", "http", "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b");
            BOOST_TEST_CSTR_EQ(u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.encoded_path().empty() );
        }
        {
            url u = urls::format("{}://[{}]", "http", "fe80::1ff:fe23:4567:890a");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://[fe80::1ff:fe23:4567:890a]");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "[fe80::1ff:fe23:4567:890a]" );
            BOOST_TEST_EQ(u.host_type(), host_type::ipv6);
            BOOST_TEST( u.encoded_path().empty() );
        }
        {
            url u = urls::format("{}:?q", "http");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http:?q");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_query(), "q" );
        }
        {
            url u = urls::format("{}:/", "http");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http:/");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/" );
            segs_equal( u.encoded_segments(), {""} );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}:path:to:joe", "mailto");
            BOOST_TEST_CSTR_EQ(u.buffer(), "mailto:path:to:joe");
            BOOST_TEST_CSTR_EQ( u.scheme(), "mailto" );
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "path:to:joe" );
            segs_equal( u.encoded_segments(), {"path:to:joe"} );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}:{}", "mailto", "joe");
            BOOST_TEST_CSTR_EQ(u.buffer(), "mailto:joe");
            BOOST_TEST_CSTR_EQ( u.scheme(), "mailto" );
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "joe" );
            segs_equal( u.encoded_segments(), {"joe"} );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}:{}/a/{}/b", "http", 'a', 'b');
            BOOST_TEST_CSTR_EQ(u.buffer(), "http:a/a/b/b");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "a/a/b/b" );
            segs_equal( u.encoded_segments(), {"a", "a", "b", "b"} );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("http://www.a.org");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://www.a.org");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "www.a.org" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://www.a.org", "http");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://www.a.org");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "www.a.org" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://{}", "http", "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://u@www.a.org", "http");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://u@www.a.org");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_user(), "u" );
            BOOST_TEST_NOT( u.has_password() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "www.a.org" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://{}@www.a.org", "http", "user");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://user@www.a.org");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_user(), "user" );
            BOOST_TEST_NOT( u.has_password() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "www.a.org" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://{}:{}@www.a.org", "http", "u", "p");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://u:p@www.a.org");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_user(), "u" );
            BOOST_TEST( u.has_password() );
            BOOST_TEST_CSTR_EQ( u.encoded_password(), "p" );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "www.a.org" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://{}:{}@{}", "http", "u", "p", "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://u:p@a.b");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_user(), "u" );
            BOOST_TEST( u.has_password() );
            BOOST_TEST_CSTR_EQ( u.encoded_password(), "p" );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://{}:{}@{}:80", "http", 'u', 'p', "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://u:p@a.b:80");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_user(), "u" );
            BOOST_TEST( u.has_password() );
            BOOST_TEST_CSTR_EQ( u.encoded_password(), "p" );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.has_port() );
            BOOST_TEST_CSTR_EQ( u.port(), "80" );
            BOOST_TEST_EQ( u.port_number(), 80 );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://{}:{}@{}:", "http", 'u', 'p', "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://u:p@a.b:");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_user(), "u" );
            BOOST_TEST( u.has_password() );
            BOOST_TEST_CSTR_EQ( u.encoded_password(), "p" );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.has_port() );
            BOOST_TEST_CSTR_EQ( u.port(), "" );
            BOOST_TEST_EQ( u.port_number(), 0 );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://{}:{}@{}:{}", "http", 'u', 'p', "a.b", 80);
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://u:p@a.b:80");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_user(), "u" );
            BOOST_TEST( u.has_password() );
            BOOST_TEST_CSTR_EQ( u.encoded_password(), "p" );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.has_port() );
            BOOST_TEST_CSTR_EQ( u.port(), "80" );
            BOOST_TEST_EQ( u.port_number(), 80 );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}://{}?k=v", "http", "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b?k=v");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_user().empty() );
            BOOST_TEST_NOT( u.has_password() );
            BOOST_TEST( u.encoded_password().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST_NOT( u.has_port() );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST( u.has_query() );
            BOOST_TEST_EQ( u.encoded_query(), "k=v" );
            params_equal( u.encoded_params(), {{"k", "v"}});
        }
        {
            url u = urls::format("{}://{}?{}", "http", "a.b", "k=v");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b?k=v");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_user().empty() );
            BOOST_TEST_NOT( u.has_password() );
            BOOST_TEST( u.encoded_password().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST_NOT( u.has_port() );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST( u.has_query() );
            BOOST_TEST_EQ( u.encoded_query(), "k=v" );
            params_equal( u.encoded_params(), {{"k", "v"}});
        }
        {
            url u = urls::format("{}://{}:{}@{}:{}/{}/{}/{}?{}", "http", 'u', 'p', "a.b", 80, 'a', 'b', 'c', "k=v");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://u:p@a.b:80/a/b/c?k=v");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_user(), "u" );
            BOOST_TEST( u.has_password() );
            BOOST_TEST_CSTR_EQ( u.encoded_password(), "p" );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.has_port() );
            BOOST_TEST_CSTR_EQ( u.port(), "80" );
            BOOST_TEST_EQ( u.port_number(), 80 );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a/b/c" );
            segs_equal( u.encoded_segments(), {"a", "b", "c"});
            BOOST_TEST( u.has_query() );
            BOOST_TEST_EQ( u.encoded_query(), "k=v" );
            params_equal( u.encoded_params(), {{"k", "v"}});
        }

        {
            url u = urls::format("/a/b/c/d");
            BOOST_TEST_CSTR_EQ(u.buffer(), "/a/b/c/d");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_user().empty() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_NOT( u.has_port() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a/b/c/d" );
            segs_equal( u.encoded_segments(), {"a", "b", "c", "d"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("/{}/b/{}/d", 'a', 'c');
            BOOST_TEST_CSTR_EQ(u.buffer(), "/a/b/c/d");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a/b/c/d" );
            segs_equal( u.encoded_segments(), {"a", "b", "c", "d"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("/{}/b/{}/d?q", 'a', 'c');
            BOOST_TEST_CSTR_EQ(u.buffer(), "/a/b/c/d?q");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a/b/c/d" );
            segs_equal( u.encoded_segments(), {"a", "b", "c", "d"});
            BOOST_TEST( u.has_query() );
            params_equal( u.encoded_params(), {{"q", ""}});
        }
        {
            url u = urls::format("/{}/b/{}/d?{}", 'a', 'c', "k=v");
            BOOST_TEST_CSTR_EQ(u.buffer(), "/a/b/c/d?k=v");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a/b/c/d" );
            segs_equal( u.encoded_segments(), {"a", "b", "c", "d"});
            BOOST_TEST( u.has_query() );
            params_equal( u.encoded_params(), {{"k", "v"}});
        }

        {
            url u = urls::format("");
            BOOST_TEST_CSTR_EQ(u.buffer(), "");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("/");
            BOOST_TEST_CSTR_EQ(u.buffer(), "/");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/" );
            BOOST_TEST( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {""});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("a");
            BOOST_TEST_CSTR_EQ(u.buffer(), "a");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "a" );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"a"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("//www.a.com");
            BOOST_TEST_CSTR_EQ(u.buffer(), "//www.a.com");
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.host(), "www.a.com" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("//{}", "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "//a.b");
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.host(), "a.b" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("abc");
            BOOST_TEST_CSTR_EQ(u.buffer(), "abc");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "abc" );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"abc"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}", "v");
            BOOST_TEST_CSTR_EQ(u.buffer(), "v");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "v" );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"v"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}", 0);
            BOOST_TEST_CSTR_EQ(u.buffer(), "0");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "0" );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"0"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("{}", 0U);
            BOOST_TEST_CSTR_EQ(u.buffer(), "0");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "0" );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"0"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("/a");
            BOOST_TEST_CSTR_EQ(u.buffer(), "/a");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a" );
            BOOST_TEST( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"a"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("/{}", 'a');
            BOOST_TEST_CSTR_EQ(u.buffer(), "/a");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a" );
            BOOST_TEST( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"a"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("/a/b/c/d");
            BOOST_TEST_CSTR_EQ(u.buffer(), "/a/b/c/d");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a/b/c/d" );
            BOOST_TEST( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"a", "b", "c", "d"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("/{}/b/{}/d", 'a', 'c');
            BOOST_TEST_CSTR_EQ(u.buffer(), "/a/b/c/d");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/a/b/c/d" );
            BOOST_TEST( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"a", "b", "c", "d"});
            BOOST_TEST_NOT( u.has_query() );
        }
        {
            url u = urls::format("a?q#f");
            BOOST_TEST_CSTR_EQ(u.buffer(), "a?q#f");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "a" );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"a"});
            BOOST_TEST( u.has_query() );
            params_equal( u.encoded_params(), {{"q", ""}});
            BOOST_TEST( u.has_fragment() );
            BOOST_TEST_CSTR_EQ( u.encoded_fragment(), "f" );
        }
        {
            url u = urls::format("a?{}#{}", 'q', 'f');
            BOOST_TEST_CSTR_EQ(u.buffer(), "a?q#f");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "a" );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"a"});
            BOOST_TEST( u.has_query() );
            params_equal( u.encoded_params(), {{"q", ""}});
            BOOST_TEST( u.has_fragment() );
            BOOST_TEST_CSTR_EQ( u.encoded_fragment(), "f" );
        }

        {
            url u = urls::format("http://www.a.com");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://www.a.com");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "www.a.com" );
            BOOST_TEST( u.encoded_path().empty() );
        }
        {
            url u = urls::format("{}://{}", "http", "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.encoded_path().empty() );
        }
        {
            url u = urls::format("{}://{}?q#f", "http", "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b?q#f");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {});
            BOOST_TEST( u.has_query() );
            params_equal( u.encoded_params(), {{"q", ""}});
            BOOST_TEST( u.has_fragment() );
            BOOST_TEST_CSTR_EQ( u.encoded_fragment(), "f" );
        }
        {
            url u = urls::format("{}://{}?{}#{}", "http", "a.b", 'q', 'f');
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b?q#f");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {});
            BOOST_TEST( u.has_query() );
            params_equal( u.encoded_params(), {{"q", ""}});
            BOOST_TEST( u.has_fragment() );
            BOOST_TEST_CSTR_EQ( u.encoded_fragment(), "f" );
        }

        {
            url u = urls::format("{}://{}", "http", "a.b");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {});
            BOOST_TEST_NOT( u.has_query() );
            BOOST_TEST_NOT( u.has_fragment() );
        }
        {
            url u = urls::format("{}", 'p');
            BOOST_TEST_CSTR_EQ(u.buffer(), "p");
            BOOST_TEST_NOT( u.has_authority() );
            BOOST_TEST( u.encoded_host().empty() );
            BOOST_TEST_NOT( u.is_path_absolute() );
            segs_equal( u.encoded_segments(), {"p"} );
            BOOST_TEST_NOT( u.has_query() );
            BOOST_TEST_NOT( u.has_fragment() );
        }

        {
            std::string p = "path/to";
            url u =
                urls::format("http{}://{}.{}.com:{}/{}/file.txt?k={}#frag-{}",
                    's',
                    "www",
                    "h/o/s/t",
                    80,
                    p,
                    "v",
                    X{});
            BOOST_TEST_CSTR_EQ(u.buffer(), "https://www.h%2Fo%2Fs%2Ft.com:80/path/to/file.txt?k=v#frag-X");
            BOOST_TEST_CSTR_EQ( u.scheme(), "https" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_NOT( u.has_userinfo() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "www.h%2Fo%2Fs%2Ft.com" );
            BOOST_TEST( u.has_port() );
            BOOST_TEST_CSTR_EQ( u.port(), "80" );
            BOOST_TEST_EQ( u.port_number(), 80 );
            BOOST_TEST( u.is_path_absolute() );
            BOOST_TEST_CSTR_EQ( u.encoded_path(), "/path/to/file.txt" );
            segs_equal( u.encoded_segments(), {"path", "to", "file.txt"});
            BOOST_TEST( u.has_query() );
            params_equal( u.encoded_params(), {{"k", "v"}});
            BOOST_TEST( u.has_fragment() );
            BOOST_TEST_CSTR_EQ( u.encoded_fragment(), "frag-X" );
        }

        BOOST_TEST_CSTR_EQ(urls::format("{}://{}?{}#{}", "http", "a.b", 'q', 'f').buffer(), "http://a.b?q#f");
        BOOST_TEST_CSTR_EQ(urls::format("{}://{}?{}#{}", "http", "a.b", 'q').buffer(), "http://a.b?q#");

        BOOST_TEST_CSTR_EQ(urls::format("{}", 'c').buffer(), "c");
        BOOST_TEST_CSTR_EQ(urls::format("//{}", ':').buffer(), "//%3A");
        BOOST_TEST_CSTR_EQ(urls::format("mailto:{}", "joe").buffer(), "mailto:joe");
        BOOST_TEST_NO_THROW(urls::format("user/{}", static_cast<unsigned int>(-1)));
        BOOST_TEST_NO_THROW(urls::format("user/{}", static_cast<unsigned long long int>(-1)));
        BOOST_TEST_CSTR_EQ(urls::format("user/{}", 1).buffer(), "user/1");
        BOOST_TEST_CSTR_EQ(urls::format("user/{}", 5678).buffer(), "user/5678");
        BOOST_TEST_CSTR_EQ(urls::format("user/{}", static_cast<long long int>(-1)).buffer(), "user/-1");
        BOOST_TEST_CSTR_EQ(urls::format("{}", X{}).buffer(), "X");

        // const& arg
        {
            const unsigned short int port{80};
            url u = urls::format("{}://{}:{}", "http", "a.b", port);
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://a.b:80");
            BOOST_TEST_CSTR_EQ( u.scheme(), "http" );
            BOOST_TEST( u.has_authority() );
            BOOST_TEST_CSTR_EQ( u.encoded_host(), "a.b" );
            BOOST_TEST( u.has_port() );
            BOOST_TEST_CSTR_EQ( u.port(), "80" );
            BOOST_TEST_EQ( u.port_number(), 80 );
            BOOST_TEST( u.encoded_path().empty() );
            BOOST_TEST_NOT( u.has_query() );
        }

        // first segment contains ':'
        BOOST_TEST_CSTR_EQ(urls::format("{}:{}", "http", "joe:").buffer(), "http:joe:");
        {
            url u = urls::format("{}", "joe:");
            BOOST_TEST_CSTR_EQ(u.buffer(), "joe%3A");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "joe%3A");
            BOOST_TEST_NOT(u.is_path_absolute());
            segs_equal(u.encoded_segments(), {"joe%3A"});
        }
        BOOST_TEST_CSTR_EQ(urls::format("{}", "::joe:").buffer(), "%3A%3Ajoe%3A");
        {
            url u = urls::format("{}", "::joe:/b:");
            BOOST_TEST_CSTR_EQ(u.buffer(), "%3A%3Ajoe%3A/b:");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "%3A%3Ajoe%3A/b:");
            BOOST_TEST_NOT(u.is_path_absolute());
            segs_equal(u.encoded_segments(), {"%3A%3Ajoe%3A", "b:"});
        }

        // path starts with "//"
        {
            url u = urls::format("{}://{}", "http", "joe");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http://joe");
            BOOST_TEST(u.has_authority());
            BOOST_TEST_CSTR_EQ(u.host(), "joe");
        }
        {
            url u = urls::format("{}:{}", "http", "//joe");
            BOOST_TEST_CSTR_EQ(u.buffer(), "http:/.//joe");
            BOOST_TEST_NOT(u.has_authority());
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/.//joe");
            BOOST_TEST(u.is_path_absolute());
            segs_equal(u.encoded_segments(), {"", "joe"});
        }
        {
            url u = urls::format("//{}", "joe");
            BOOST_TEST_CSTR_EQ(u.buffer(), "//joe");
            BOOST_TEST(u.has_authority());
            BOOST_TEST_CSTR_EQ(u.host(), "joe");
        }
        {
            url u = urls::format("{}", "//joe");
            BOOST_TEST_CSTR_EQ(u.buffer(), "/.//joe");
            BOOST_TEST_NOT(u.has_authority());
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/.//joe");
            BOOST_TEST(u.is_path_absolute());
            segs_equal(u.encoded_segments(), {"", "joe"});
        }

        // invalid format strings
        BOOST_TEST_THROWS(urls::format("{:"), system::system_error);
        BOOST_TEST_THROWS(urls::format("{}://www.a.com", "1nvalid scheme"), system::system_error);
        BOOST_TEST_THROWS(urls::format("{}://{}:{}@{}:a", "http", 'u', 'p', "a.b"), system::system_error);
        BOOST_TEST_THROWS(urls::format("{}://[", "http"), system::system_error);
        BOOST_TEST_THROWS(urls::format("{://"), system::system_error);
        BOOST_TEST_THROWS(urls::format("http:%"), system::system_error);
        BOOST_TEST_THROWS(urls::format("{}:\\", "A"), system::system_error);

        // static_url
        {
            static_url<30> u;
            urls::format_to(u, "{}://{}", "https", "www.boost.org");
            BOOST_TEST_CSTR_EQ(u.buffer(), "https://www.boost.org");
            BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            BOOST_TEST(u.has_authority());
            BOOST_TEST_CSTR_EQ(u.host(), "www.boost.org");
            BOOST_TEST(u.encoded_path().empty());
            BOOST_TEST_NOT(u.is_path_absolute());
        }
        {
            static_url<30> u;
            urls::format_to(u, "{s}://{h}", {{"s", "https"}, {"h", "www.boost.org"}});
            BOOST_TEST_CSTR_EQ(u.buffer(), "https://www.boost.org");
        }
        {
            static_url<10> u;
            BOOST_TEST_THROWS(urls::format_to(u, "{}://{}", "https", "www.boost.org"), system::system_error);
        }


        // escaped '{' always throws because '{'s are not allowed in URLs
        BOOST_TEST_THROWS(urls::format("{scheme}:{path}/{{}", "mailto", 'a'), system::system_error);

        // "{}" with no format arg is ignored
        BOOST_TEST_CSTR_EQ(urls::format("/{}/{}/{}", 'a', 'b').buffer(), "/a/b/");

        // format specs
        {
            BOOST_TEST_CSTR_EQ(urls::format("{:}", 'a').buffer(), "a");
            BOOST_TEST_CSTR_EQ(urls::format("{:c}", 'a').buffer(), "a");
            BOOST_TEST_CSTR_EQ(urls::format("{:^1s}", 'a').buffer(), "a");
            BOOST_TEST_CSTR_EQ(urls::format("{:^3s}", 'a').buffer(), "%20a%20");
            BOOST_TEST_CSTR_EQ(urls::format("{:.^5s}", 'a').buffer(), "..a..");
            BOOST_TEST_CSTR_EQ(urls::format("{:.<5s}", 'a').buffer(), "a....");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>5s}", 'a').buffer(), "....a");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{}s}", 'a', 5).buffer(), "....a");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{1}s}", 'a', 5).buffer(), "....a");
            BOOST_TEST_CSTR_EQ(urls::format("{0:.>{2}s}/{1}", 'a', 'b', 5).buffer(), "....a/b");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{b}s}", 'a', arg("b", 5)).buffer(), "....a");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{b}s}", 'a', arg("b", '5')).buffer(), "....a");
            BOOST_TEST_THROWS(urls::format("{:cx}", 'a'), system::system_error);

            BOOST_TEST_CSTR_EQ(urls::format("{:d}", 99).buffer(), "99");
            BOOST_TEST_CSTR_EQ(urls::format("{:#d}", 99).buffer(), "99");
            BOOST_TEST_CSTR_EQ(urls::format("{:^1d}", 99).buffer(), "99");
            BOOST_TEST_CSTR_EQ(urls::format("{:+d}", 99).buffer(), "+99");
            BOOST_TEST_CSTR_EQ(urls::format("{: d}", std::size_t(99)).buffer(), "%2099");
            BOOST_TEST_CSTR_EQ(urls::format("{:+d}", std::size_t(99)).buffer(), "+99");
            BOOST_TEST_CSTR_EQ(urls::format("{: d}", 99).buffer(), "%2099");
            BOOST_TEST_CSTR_EQ(urls::format("{:^6d}", 99).buffer(), "%20%2099%20%20");
            BOOST_TEST_CSTR_EQ(urls::format("{:.^6d}", 99).buffer(), "..99..");
            BOOST_TEST_CSTR_EQ(urls::format("{:.<6d}", 99).buffer(), "99....");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>6d}", 99).buffer(), "....99");
            BOOST_TEST_CSTR_EQ(urls::format("{:.^6d}", std::size_t(99)).buffer(), "..99..");
            BOOST_TEST_CSTR_EQ(urls::format("{:.<6d}", std::size_t(99)).buffer(), "99....");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>6d}", std::size_t(99)).buffer(), "....99");
            BOOST_TEST_CSTR_EQ(urls::format("{:>06d}", 99).buffer(), "000099");
            BOOST_TEST_CSTR_EQ(urls::format("{:>06d}", std::size_t(99)).buffer(), "000099");
            BOOST_TEST_CSTR_EQ(urls::format("{:>+06d}", 99).buffer(), "+00099");
            BOOST_TEST_CSTR_EQ(urls::format("{:> 06d}", 99).buffer(), "%2000099");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{}d}", 99, 6).buffer(), "....99");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{1}d}", 99, 6).buffer(), "....99");
            BOOST_TEST_CSTR_EQ(urls::format("{0:.>{2}d}/{1}", 99, 'b', 6).buffer(), "....99/b");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{b}d}", 99, arg("b", 6)).buffer(), "....99");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{}d}", std::size_t(99), 6).buffer(), "....99");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{1}d}", std::size_t(99), 6).buffer(), "....99");
            BOOST_TEST_CSTR_EQ(urls::format("{0:.>{2}d}/{1}", std::size_t(99), 'b', 6).buffer(), "....99/b");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{b}d}", std::size_t(99), arg("b", 6)).buffer(), "....99");
            BOOST_TEST_CSTR_EQ(urls::format("{:.>{b}d}", 99, arg("b", '6')).buffer(), "....99");
            BOOST_TEST_THROWS(urls::format("{:dx}", 99), system::system_error);

        }

        // positional arguments
        {
            BOOST_TEST_CSTR_EQ(urls::format("{}/{}/{}", 'a', 'b', 'c').buffer(), "a/b/c");
            BOOST_TEST_CSTR_EQ(urls::format("{0}/{1}/{2}", 'a', 'b', 'c').buffer(), "a/b/c");
            BOOST_TEST_CSTR_EQ(urls::format("{2}/{1}/{0}", 'a', 'b', 'c').buffer(), "c/b/a");
            BOOST_TEST_CSTR_EQ(urls::format("//{0}{1}{0}.com", "abra", "cad").buffer(), "//abracadabra.com");
            BOOST_TEST_CSTR_EQ(urls::format("{0}/{}/{}/{1}/{}", 'a', 'b', 'c').buffer(), "a/a/b/b/c");
            BOOST_TEST_CSTR_EQ(urls::format("https://www.boost.org/en/{0}-{1}", "fast", "library").buffer(), "https://www.boost.org/en/fast-library");
            BOOST_TEST_CSTR_EQ(urls::format("https://www.boost.org/en/{1}-{0}", "rapida", "biblioteca").buffer(), "https://www.boost.org/en/biblioteca-rapida");
        }

        // named arguments
        {
            // examples from openAPI
            // https://swagger.io/specification/#paths-object
            BOOST_TEST_CSTR_EQ(urls::format("/pets/{petId}", arg("petId", 30)).buffer(), "/pets/30");
            BOOST_TEST_CSTR_EQ(urls::format("/pets/{name}", arg("name", "ted")).buffer(), "/pets/ted");
            BOOST_TEST_CSTR_EQ(urls::format("/users/{userid}/address", arg("userid", 30)).buffer(), "/users/30/address");
            BOOST_TEST_CSTR_EQ(urls::format("~12.0~1repositories~1{username}/get", arg("username", "boostorg")).buffer(), "~12.0~1repositories~1boostorg/get");
            BOOST_TEST_CSTR_EQ(
                urls::format(
                    "https://{username}.gigantic-server.com:{port}/{basePath}/{path}",
                    arg("username", "joe"),
                    arg("port", 80),
                    arg("basePath", "v2"),
                    arg("path", "index.html")).buffer(),
                "https://joe.gigantic-server.com:80/v2/index.html");
            BOOST_TEST_CSTR_EQ(
                urls::format(
                    "https://{username}.gigantic-server.com:{port}/{basePath}/{path}",
                    arg("basePath", "v2"),
                    arg("path", "index.html"),
                    arg("port", 80),
                    arg("username", "joe")).buffer(),
                "https://joe.gigantic-server.com:80/v2/index.html");
            BOOST_TEST_CSTR_EQ(urls::format("/pets/{petId}/{unknown}", arg("petId", 30)).buffer(), "/pets/30/");
            // initializer-list overload
            BOOST_TEST_CSTR_EQ(
                urls::format(
                    "https://{username}.gigantic-server.com:{port}/{basePath}/{path}",
                    {{"basePath", "v2"}, {"path", "index.html"}, {"port", 80}, {"username", "joe"}}).buffer(),
                "https://joe.gigantic-server.com:80/v2/index.html");
        }

    }

    void
    run()
    {
        // I have spent a lot of time on this and have no
        // idea how to fix this bug in GCC 4.8 and GCC 5.0
        // without help from the pros.
#if !BOOST_WORKAROUND( BOOST_GCC_VERSION, < 60000 )
        testFormat();
#endif
    }
};

TEST_SUITE(format_test, "boost.url.format");

// Test extending supported types in implementation
namespace detail {
template <>
struct formatter<X>
{
public:
    char const*
    parse(format_parse_context& ctx) const
    {
        return formatter<ignore_format>::parse_empty_spec(
            ctx.begin(), ctx.end());
    }

    std::size_t
    measure(
        X,
        measure_context& ctx,
        grammar::lut_chars const& cs) const
    {
        return ctx.out() + measure_one('X', cs);
    }

    char*
    format(X, format_context& ctx, grammar::lut_chars const& cs) const
    {
        char* o = ctx.out();
        encode_one(o, 'X', cs);
        return o;
    }
};
}

} // urls
} // boost
