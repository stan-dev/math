//
// Copyright (c) 2023 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/url.hpp>
#include <boost/url/parse.hpp>

#include "test_suite.hpp"

#ifdef BOOST_TEST_CSTR_EQ
#undef BOOST_TEST_CSTR_EQ
#define BOOST_TEST_CSTR_EQ(expr1,expr2) \
    BOOST_TEST_EQ( boost::urls::detail::to_sv(expr1), boost::urls::detail::to_sv(expr2) )
#endif

namespace boost {
namespace urls {

struct ada_test
{
    static
    void
    setHostShouldReturnFalseSometimes()
    {
        system::result<url> r = parse_uri("mailto:a@b.com");
        BOOST_TEST_CSTR_EQ(r->scheme(), "mailto");
        BOOST_TEST_CSTR_EQ(r->path(), "a@b.com");
        // This is an invalid operation in ada, which is reasonable since
        // the path is relative, so it can't have a host.
        // Boost.URL will encode the host, add any required prefixes, and
        // make the path absolute.
        BOOST_TEST_NO_THROW(r->set_encoded_host("something"));
        BOOST_TEST_CSTR_EQ(r->scheme(), "mailto");
        BOOST_TEST_CSTR_EQ(r->host(), "something");
        BOOST_TEST_CSTR_EQ(r->path(), "/a@b.com");
        BOOST_TEST_CSTR_EQ(r->buffer(), "mailto://something/a@b.com");
    }

    static
    void
    setHostShouldReturnTrueSometimes()
    {
        system::result<url> r = parse_uri("https://www.google.com");
        BOOST_TEST_NO_THROW(r->set_encoded_host("something"));
    }

    static
    void
    setHostnameShouldReturnFalseSometimes()
    {
        system::result<url> r = parse_uri("mailto:a@b.com");
        // This is an invalid operation in ada, while Boost.URL will
        // encode the host and adjust the URL.
        BOOST_TEST_NO_THROW(r->set_host_address("something"));
        BOOST_TEST_CSTR_EQ(r->scheme(), "mailto");
        BOOST_TEST_CSTR_EQ(r->host(), "something");
        BOOST_TEST_CSTR_EQ(r->path(), "/a@b.com");
        BOOST_TEST_CSTR_EQ(r->buffer(), "mailto://something/a@b.com");
    }

    static
    void
    setHostnameShouldReturnTrueSometimes()
    {
        system::result<url> r = parse_uri("https://www.google.com");
        BOOST_TEST_NO_THROW(r->set_host_address("something"));
    }

    static
    void
    readme1()
    {
        system::result<url> u = parse_uri("https://www.google.com");
        BOOST_TEST(u);
    }

    static
    void
    readme2()
    {
        system::result<url> u = parse_uri("https://www.google.com");
        if (!BOOST_TEST(u))
            return;
        u->set_encoded_user("username");
        u->set_encoded_password("password");
        // Ada would add a "/" path at the end of the authority
        // Boost.URL keeps the original path
        BOOST_TEST_CSTR_EQ(u->buffer(), "https://username:password@www.google.com");
    }

    static
    void
    readme3()
    {
        system::result<url> u = parse_uri("https://www.google.com");
        if (!BOOST_TEST(u))
            return;
        u->set_scheme("wss");
        BOOST_TEST_CSTR_EQ(u->scheme(), "wss");
        // Ada would add a "/" path at the end of the authority
        // Boost.URL keeps the original path
        BOOST_TEST_CSTR_EQ(u->buffer(), "wss://www.google.com");
    }

    static
    void
    readme4()
    {
        system::result<url> u = parse_uri("https://www.google.com");
        if (!BOOST_TEST(u))
            return;
        u->set_encoded_host("github.com");
        BOOST_TEST_CSTR_EQ(u->encoded_host(), "github.com");
    }

    static
    void
    readme5()
    {
        system::result<url> u = parse_uri("https://www.google.com");
        if (!BOOST_TEST(u))
            return;
        u->set_port("8080");
        BOOST_TEST_CSTR_EQ(u->port(), "8080");
    }

    static
    void
    readme6()
    {
        system::result<url> u = parse_uri("https://www.google.com");
        if (!BOOST_TEST(u))
            return;
        u->set_encoded_path("/my-super-long-path");
        BOOST_TEST_CSTR_EQ(u->encoded_path(), "/my-super-long-path");
    }

    static
    void
    readme7()
    {
        system::result<url> u = parse_uri("https://www.google.com");
        if (!BOOST_TEST(u))
            return;
        u->set_encoded_query("target=self");
        // Ada would return the "?" prefix with the query
        // Boost.URL only returns the query
        BOOST_TEST_CSTR_EQ(u->encoded_query(), "target=self");
    }

    static
    void
    readme8()
    {
        system::result<url> u = parse_uri("https://www.google.com");
        if (!BOOST_TEST(u))
            return;
        u->set_encoded_fragment("is-this-the-real-life");
        // Ada would return the "#" prefix with the fragment
        // Boost.URL only returns the fragment
        BOOST_TEST_CSTR_EQ(u->encoded_fragment(), "is-this-the-real-life");
    }

    static
    void
    nodejs1()
    {
        auto base = parse_uri("http://other.com/");
        if (!BOOST_TEST(base))
            return;
        url u;
        auto r = resolve(*base, url_view("http://GOOgoo.com"), u);
        BOOST_TEST(r);
    }

    static
    void
    nodejs2()
    {
        // Ada supports URLs with whitespaces so that these whitespaces
        // are not considered part of the path
        // Boost.URL doesn't support whitespaces in the URL
        system::result<url> u = parse_uri("data:space    ?test");
        BOOST_TEST_NOT(u);
        u = parse_uri("data:space%20%20%20%20?test");
        if (!BOOST_TEST(u))
            return;
        BOOST_TEST_CSTR_EQ(u->encoded_query(), "test");
        u->set_encoded_query("");
        BOOST_TEST_CSTR_EQ(u->encoded_query(), "");
        BOOST_TEST_NOT(u->encoded_path() == "space");
        BOOST_TEST_CSTR_EQ(u->encoded_path(), "space%20%20%20%20");
        BOOST_TEST_CSTR_EQ(u->path(), "space    ");
        // Ada would remove the query is it is set to ""
        // Boost.URL keeps the empty query
        BOOST_TEST_NOT(u->buffer() == "data:space%20%20%20%20");
        BOOST_TEST_CSTR_EQ(u->buffer(), "data:space%20%20%20%20?");
        u->remove_query();
        BOOST_TEST_CSTR_EQ(u->buffer(), "data:space%20%20%20%20");
    }

    static
    void
    nodejs3()
    {
        // Ada supports URLs with whitespaces so that these whitespaces
        // are not considered part of the path
        // Boost.URL doesn't support whitespaces in the URL
        system::result<url> u = parse_uri("data:space    ?test#test");
        BOOST_TEST_NOT(u);
        u = parse_uri("data:space%20%20%20%20?test#test");
        if (!BOOST_TEST(u))
            return;
        BOOST_TEST_CSTR_EQ(u->encoded_query(), "test");
        BOOST_TEST_CSTR_EQ(u->encoded_fragment(), "test");
        u->set_encoded_query("");
        BOOST_TEST_CSTR_EQ(u->encoded_query(), "");
        BOOST_TEST_CSTR_EQ(u->path(), "space    ");
        BOOST_TEST_CSTR_EQ(u->buffer(), "data:space%20%20%20%20?#test");
    }

    static
    void
    nodejs4()
    {
        system::result<url> u = parse_uri("file:///var/log/system.log");
        if (!BOOST_TEST(u))
            return;
        // This test uses set_href in ada.
        // boost.url does not support this operation.
        u = parse_uri("http://0300.168.0xF0");
        if (!BOOST_TEST(u))
            return;
        // Ada returns ":" with the scheme. Boost.URL doesn't.
        BOOST_TEST_NOT(u->scheme() == "http:");
        BOOST_TEST_CSTR_EQ(u->scheme(), "http");
        // Ada converts "0300.168.0xF0" to "192.168.0.240"
        // Boost.URL doesn't, although "0300.168.0xF0" is a valid reg-name
        BOOST_TEST_NOT(u->buffer() == "http://192.168.0.240/");
        BOOST_TEST_CSTR_EQ(u->buffer(), "http://0300.168.0xF0");
        BOOST_TEST_CSTR_EQ(u->encoded_host_address(), "0300.168.0xF0");
    }

    static
    void
    adaBasicTests()
    {
        setHostShouldReturnFalseSometimes();
        setHostShouldReturnTrueSometimes();
        setHostnameShouldReturnFalseSometimes();
        setHostnameShouldReturnTrueSometimes();
        readme1();
        readme2();
        readme3();
        readme4();
        readme5();
        readme6();
        readme7();
        readme8();
        nodejs1();
        nodejs2();
        nodejs3();
        nodejs4();
    }

    static
    void
    urlTestData()
    {
        // Based on http://trac.webkit.org/browser/trunk/LayoutTests/fast/url/script-tests/segments.js
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://example\t.\norg");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://user:pass@foo:21/bar;par?b#c");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://user:pass@foo:21/bar;par?b#c"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "user");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "21");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/bar;par");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "b");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "c");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https://test:@test");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "https://test@test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https://:@test");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "https://test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-special://test:@test/x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-special://test@test/x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-special");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/x");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-special://:@test/x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-special://test/x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-special");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/x");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:foo.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/foo.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/foo.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("\t   :foo.com   \n");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference(" foo.com  ");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("a:\t foo.com");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("http://f:21/ b ? d # e ");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("lolscheme:x x#x x");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://f:/c");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://f/c"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "f");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/c");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://f:0/c");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://f:0/c"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "f");
            // Boost.URL does not compress port values
            BOOST_TEST_EQ(u.port_number(), 0);
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/c");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://f:00000000000000/c");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://f:0/c"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "f");
            // Boost.URL does not compress port values
            BOOST_TEST_EQ(u.port_number(), 0);
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/c");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://f:00000000000000000000080/c");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://f/c"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "f");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/c");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http://f:b/c");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("http://f: /c");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://f:\n/c");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http://f:fifty-two/c");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http://f:999999/c");
            // Boost.URL does not fail on port overflow
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("non-special://f:999999/c");
            // Boost.URL does not fail on port overflow
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("http://f: 21 / b ? d # e ");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not fail on empty references
            system::result<url> ref = parse_uri_reference("");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("  \t");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference(":foo.com/");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference(":foo.com\\");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference(":");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference(":a");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference(":/");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference(":\\");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference(":#");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/bar#"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/bar#/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("#\\");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#;?");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/bar#;?"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), ";?");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("?");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/bar?"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference(":23");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/:23");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/:23"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/:23");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("\\x");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("\\\\x\\hello");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("::");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("::23");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("foo://");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "foo://"
            BOOST_TEST_CSTR_EQ(u.scheme(), "foo");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://a:b@c:29/d");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://a:b@c:29/d"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "a");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "c");
            BOOST_TEST_CSTR_EQ(u.port(), "29");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/d");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http::@c:29");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/:@c:29"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/:@c:29");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "&" in username
            // Expected username in Ada: "&a"
            system::result<url> ref = parse_uri_reference("http://&a:foo(b]c@d:2/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not passwords with unencoded "@"
            // Expected password in Ada: "%3A%40c"
            system::result<url> ref = parse_uri_reference("http://::@c@d:2");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://foo.com:b@d/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://foo.com:b@d/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "foo.com");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "d");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://foo.com/\\@");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http:\\\\foo.com\\");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http:\\\\a\\b:c\\d@foo.com\\");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("foo:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "foo:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("foo:/bar.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "foo:/bar.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/bar.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("foo://///////");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "foo://///////"
            BOOST_TEST_CSTR_EQ(u.scheme(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "///////");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("foo://///////bar.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "foo://///////bar.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "///////bar.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("foo:////://///");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "foo:////://///"
            BOOST_TEST_CSTR_EQ(u.scheme(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//://///");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("c:/foo");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "c:/foo"
            BOOST_TEST_CSTR_EQ(u.scheme(), "c");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//foo/bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://foo/bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/bar");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://foo/path;a??e#f#g");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://foo/path;a??e#f#g"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/path;a");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "?e");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "f#g");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://foo/abcd?efgh?ijkl");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://foo/abcd?efgh?ijkl"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/abcd");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "efgh?ijkl");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://foo/abcd#foo?bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://foo/abcd#foo?bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/abcd");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "foo?bar");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "[" in paths
            system::result<url> ref = parse_uri_reference("[61:24:74]:98");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "[" in paths
            system::result<url> ref = parse_uri_reference("http:[61:27]/:foo");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http://[1::2]:3:4");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http://2001::1");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http://2001::1]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http://2001::1]:80");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://[2001::1]");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://[2001::1]"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://[::127.0.0.1]");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://[::7f00:1]"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "[::7f00:1]"
            BOOST_TEST_NOT(u.encoded_host() == "[::7f00:1]");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http://[::127.0.0.1.]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://[0:0:0:0:0:0:13.1.68.3]");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://[::d01:4403]"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "[::d01:4403]"
            BOOST_TEST_NOT(u.encoded_host() == "[::d01:4403]");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://[2001::1]:80");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://[2001::1]"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftp:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("madeupscheme:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "madeupscheme:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "madeupscheme");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("file://example:1/");
            // Boost.URL does validate individual schemes
            // Ada and Node.js fail on "file" URLs with ports
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("file://example:test/");
            // Boost.URL does validate individual schemes
            // Ada and Node.js fail on "file" URLs with ports
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("file://example%/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("file://[example]/");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftps:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ftps:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ftps");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("gopher:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "gopher:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "gopher");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ws:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("wss:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("data:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "data:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "data");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("javascript:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "javascript:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "javascript");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("mailto:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "mailto:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "mailto");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftp:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("madeupscheme:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "madeupscheme:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "madeupscheme");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftps:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ftps:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ftps");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("gopher:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "gopher:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "gopher");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ws:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("wss:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("data:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "data:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "data");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("javascript:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "javascript:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "javascript");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("mailto:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "mailto:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "mailto");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/a/b/c");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/a/b/c"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/a/b/c");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("/a/ /c");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/a%2fc");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/a%2fc"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/a%2Fc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/a%2Fc");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/a/%2f/c");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/a/%2f/c"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/a/%2F/c");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/a/%2F/c");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("#β");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("data:text/html,test#test");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "data:text/html,test#test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "data");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "text/html,test");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "test");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("tel:1234567890");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "tel:1234567890"
            BOOST_TEST_CSTR_EQ(u.scheme(), "tel");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "1234567890");
        }();
        // Based on https://felixfbecker.github.io/whatwg-url-custom-host-repro/
        []{
            system::result<url> base = parse_uri("http://example.org/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ssh://example.com/foo/bar.git");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ssh://example.com/foo/bar.git"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ssh");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar.git");
        }();
        // Based on http://trac.webkit.org/browser/trunk/LayoutTests/fast/url/file.html
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("file:c:\\foo\\bar.html");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("  File:c|////foo\\bar.html");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("C|/foo/bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("/C|\\foo\\bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("//C|/foo/bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//server/file");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://server/file"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "server");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/file");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("\\\\server\\file");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("/\\server/file");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///foo/bar.txt");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///foo/bar.txt"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar.txt");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///home/me");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///home/me"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/home/me");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("///");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("///test");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://test");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://localhost");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://localhost/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://localhost/test");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("test");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///tmp/mock/test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/tmp/mock/test");
        }();
        []{
            system::result<url> base = parse_uri("file:///tmp/mock/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:test");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///tmp/mock/test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/tmp/mock/test");
        }();
        // Based on http://trac.webkit.org/browser/trunk/LayoutTests/fast/url/script-tests/path.js
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/././foo");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/./.foo");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/.foo"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/.foo");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/.");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/./");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/bar/..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/bar/../");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/..bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/..bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/..bar");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/bar/../ton");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/ton"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/ton");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/bar/../ton/../../a");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/a"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/a");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/../../..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['foo', '..', '..', '..']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/../../../ton");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/ton"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['foo', '..', '..', '..', 'ton']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/ton");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/%2e");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://example.com/foo/%2e%2");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/%2e./%2e%2e/.%2e/%2e.bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/%2e.bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['foo', '%2e.', '%2e%2e', '.%2e', '%2e.bar']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%2E.bar");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com////../..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com//"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/bar//../..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo/bar//..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo/bar/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/%20foo");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/%20foo"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%20foo");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://example.com/foo%");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://example.com/foo%2");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://example.com/foo%2zbar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://example.com/foo%2Â©zbar");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo%41%7a");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo%41%7a"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo%41%7A");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/fooAz");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://example.com/foo\t\x91" "%91");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/foo%00%51");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/foo%00%51"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo%00%51");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo%00Q");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/(%28:%3A%29)");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/(%28:%3A%29)"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/(%28:%3A%29)");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/((::))");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/%3A%3a%3C%3c");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/%3A%3a%3C%3c"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%3A%3A%3C%3C");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/::%3C%3C");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://example.com/foo\tbar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://example.com\\\\foo\\\\bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/%7Ffp3%3Eju%3Dduvgw%3Dd");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/%7Ffp3%3Eju%3Dduvgw%3Dd"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%7Ffp3%3Eju%3Dduvgw%3Dd");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%7Ffp3%3Eju=duvgw=d");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.com/@asdf%40");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com/@asdf%40"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/@asdf%40");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/@asdf@");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://example.com/你好你好");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://example.com/‥/foo");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://example.com/\ufeff/foo");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://example.com/\u202e/foo/\u202d/bar");
            // BOOST_TEST_NOT(ref);
        }();
        // Based on http://trac.webkit.org/browser/trunk/LayoutTests/fast/url/script-tests/relative.js
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://www.google.com/foo?bar=baz#");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.google.com/foo?bar=baz#"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.google.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "bar=baz");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://www.google.com/foo?bar=baz# »");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("data:test# »");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://www.google.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.google.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.google.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://192.0x00A80001");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://192.168.0.1"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "192.168.0.1"
            BOOST_TEST_NOT(u.encoded_host() == "192.168.0.1");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://www/foo%2Ehtml");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www/foo%2Ehtml"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo%2Ehtml");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo.html");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://www/foo/%2E/html");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www/foo/html"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/html");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://user:pass@/");
            // Boost.URL and Node.js do not fail
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://%25DOMAIN:foobar@foodomain.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://%25DOMAIN:foobar@foodomain.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "%25DOMAIN");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foodomain.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http:\\\\www.google.com\\foo");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://foo:80/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://foo:81/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://foo:81/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "81");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("httpa://foo:80/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "httpa://foo:80/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "httpa");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "80");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://foo:-80/");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https://foo:443/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "https://foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https://foo:80/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "https://foo:80/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "80");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftp://foo:21/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ftp://foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ftp");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftp://foo:80/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ftp://foo:80/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ftp");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "80");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("gopher://foo:70/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "gopher://foo:70/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "gopher");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "70");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("gopher://foo:443/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "gopher://foo:443/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "gopher");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "443");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ws://foo:80/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ws://foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ws");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ws://foo:81/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ws://foo:81/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ws");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "81");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ws://foo:443/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ws://foo:443/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ws");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "443");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ws://foo:815/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ws://foo:815/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ws");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "815");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("wss://foo:80/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "wss://foo:80/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "wss");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "80");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("wss://foo:81/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "wss://foo:81/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "wss");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "81");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("wss://foo:443/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "wss://foo/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "wss");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("wss://foo:815/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "wss://foo:815/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "wss");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foo");
            BOOST_TEST_CSTR_EQ(u.port(), "815");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftp:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("madeupscheme:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "madeupscheme:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "madeupscheme");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftps:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ftps:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ftps");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("gopher:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "gopher:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "gopher");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ws:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("wss:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("data:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "data:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "data");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("javascript:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "javascript:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "javascript");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("mailto:/example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "mailto:/example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "mailto");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftp:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("madeupscheme:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "madeupscheme:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "madeupscheme");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ftps:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ftps:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ftps");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("gopher:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "gopher:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "gopher");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ws:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("wss:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "example.com"
            BOOST_TEST_NOT(u.encoded_host() == "example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("data:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "data:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "data");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("javascript:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "javascript:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "javascript");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("mailto:example.com/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "mailto:example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "mailto");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "example.com/");
        }();
        // Based on http://trac.webkit.org/browser/trunk/LayoutTests/fast/url/segments-userinfo-vs-host.html
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "www.example.com"
            BOOST_TEST_NOT(u.encoded_host() == "www.example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:/@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "www.example.com"
            BOOST_TEST_NOT(u.encoded_host() == "www.example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:a:b@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "www.example.com"
            BOOST_TEST_NOT(u.encoded_host() == "www.example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:/a:b@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "www.example.com"
            BOOST_TEST_NOT(u.encoded_host() == "www.example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://a:b@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://a:b@www.example.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "a");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://@pple.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://pple.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "pple.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http::b@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "www.example.com"
            BOOST_TEST_NOT(u.encoded_host() == "www.example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:/:b@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "www.example.com"
            BOOST_TEST_NOT(u.encoded_host() == "www.example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://:b@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://:b@www.example.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http:/:@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://user@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http:@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http:/@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("https:@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http:a:b@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http:/a:b@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://a:b@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http::@/www.example.com");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:a:@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "www.example.com"
            BOOST_TEST_NOT(u.encoded_host() == "www.example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:/a:@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL does not convert paths into hostnames
            // Expected hostname in Ada: "www.example.com"
            BOOST_TEST_NOT(u.encoded_host() == "www.example.com");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://a:@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://a@www.example.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "a");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://www.@pple.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.@pple.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "www.");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "pple.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http:@:www.example.com");
            // Boost.URL does not fail on paths that resemble hosts
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http:/@:www.example.com");
            // Boost.URL does not fail on paths that resemble hosts
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://@:www.example.com");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://:@www.example.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        // Others
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/test.txt");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com/test.txt"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test.txt");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference(".");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("test.txt");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com/test.txt"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test.txt");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("./test.txt");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com/test.txt"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test.txt");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("../test.txt");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com/test.txt"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', 'test.txt']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test.txt");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("../aaa/test.txt");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com/aaa/test.txt"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', 'aaa', 'test.txt']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/aaa/test.txt");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("../../test.txt");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example.com/test.txt"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example.com");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '..', 'test.txt']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test.txt");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("中/test.txt");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://www.example2.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example2.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example2.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//www.example2.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://www.example2.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "www.example2.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:...");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///..."
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL normalization does not treat "file:" differently
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/...");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://www.example.com/test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:a");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///a"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL normalization does not treat "file:" differently
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/a");
        }();
        // Based on http://trac.webkit.org/browser/trunk/LayoutTests/fast/url/host.html
        // Basic canonicalization, uppercase should be converted to lowercase
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://ExAmPlE.CoM");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("http://example example.com");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("http://Goo%20 goo%7C|.com");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[:]");
            BOOST_TEST_NOT(base && ref);
        }();
        // U+3000 is mapped to U+0020 (space) which is disallowed
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://GOO\xa0" "\u3000goo.com");
            // BOOST_TEST_NOT(ref);
        }();
        // Other types of space (no-break, zero-width, zero-width-no-break) are name-prepped away to nothing. U+200B, U+2060, and U+FEFF, are ignored
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://GOO\u200b\u2060\ufeffgoo.com");
            // BOOST_TEST_NOT(ref);
        }();
        // Leading and trailing C0 control or space
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("\x00" "\x1b" "\x04" "\x12" " http://example.com/\x1f" " \r ");
            // BOOST_TEST_NOT(ref);
        }();
        // Ideographic full stop (full-width period for Chinese, etc.) should be treated as a dot. U+3002 is mapped to U+002E (dot)
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://www.foo。bar.com");
            // BOOST_TEST_NOT(ref);
        }();
        // Invalid unicode characters should fail... U+FDD0 is disallowed; %ef%b7%90 is U+FDD0
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://\ufdd0zyx.com");
            // BOOST_TEST_NOT(ref);
        }();
        // This is the same as previous but escaped
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("http://%ef%b7%90zyx.com");
            BOOST_TEST(ref);
        }();
        // U+FFFD
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("https://�");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("https://%EF%BF%BD");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("https://x/�?�#�");
            // BOOST_TEST_NOT(ref);
        }();
        // Domain is ASCII, but a label is invalid IDNA
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid IDNA
            system::result<url> ref = parse_uri_reference("http://a.b.c.xn--pokxncvks");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid IDNA
            system::result<url> ref = parse_uri_reference("http://10.0.0.xn--pokxncvks");
            BOOST_TEST(ref);
        }();
        // IDNA labels should be matched case-insensitively
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid IDNA
            system::result<url> ref = parse_uri_reference("http://a.b.c.XN--pokxncvks");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid IDNA
            system::result<url> ref = parse_uri_reference("http://a.b.c.Xn--pokxncvks");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid IDNA
            system::result<url> ref = parse_uri_reference("http://10.0.0.XN--pokxncvks");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid IDNA
            system::result<url> ref = parse_uri_reference("http://10.0.0.xN--pokxncvks");
            BOOST_TEST(ref);
        }();
        // Test name prepping, fullwidth input should be converted to ASCII and NOT IDN-ized. This is 'Go' in fullwidth UTF-8/UTF-16.
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://Ｇｏ.com");
            // BOOST_TEST_NOT(ref);
        }();
        // URL spec forbids the following. https://www.w3.org/Bugs/Public/show_bug.cgi?id=24257
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://％４１.com");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("http://%ef%bc%85%ef%bc%94%ef%bc%91.com");
            BOOST_TEST(ref);
        }();
        // ...%00 in fullwidth should fail (also as escaped UTF-8 input)
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://％００.com");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("http://%ef%bc%85%ef%bc%90%ef%bc%90.com");
            BOOST_TEST(ref);
        }();
        // Basic IDN support, UTF-8 and UTF-16 input should be converted to IDN
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://你好你好");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("https://faß.ExAmPlE/");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("sc://faß.ExAmPlE/");
            // BOOST_TEST_NOT(ref);
        }();
        // Invalid escaped characters should fail and the percents should be escaped. https://www.w3.org/Bugs/Public/show_bug.cgi?id=24191
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://%zz%66%a.com");
            BOOST_TEST_NOT(ref);
        }();
        // If we get an invalid character that has been escaped.
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://%25");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://hello%00");
            BOOST_TEST(ref);
        }();
        // Escaped numbers should be treated like IP addresses if they are.
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://%30%78%63%30%2e%30%32%35%30.01");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://192.168.0.1"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "192.168.0.1"
            BOOST_TEST_NOT(u.encoded_host() == "192.168.0.1");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://%30%78%63%30%2e%30%32%35%30.01%2e");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://192.168.0.1"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "192.168.0.1"
            BOOST_TEST_NOT(u.encoded_host() == "192.168.0.1");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://192.168.0.257");
            BOOST_TEST(ref);
        }();
        // Invalid escaping in hosts causes failure
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://%3g%78%63%30%2e%30%32%35%30%2E.01");
            BOOST_TEST_NOT(ref);
        }();
        // A space in a host causes failure
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("http://192.168.0.1 hello");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("https://x x:12");
            BOOST_TEST_NOT(ref);
        }();
        // Fullwidth and escaped UTF-8 fullwidth should still be treated as IP
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://０Ｘｃ０．０２５０．０１");
            // BOOST_TEST_NOT(ref);
        }();
        // Domains with empty labels
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://./");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://./"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), ".");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://../");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://../"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "..");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        // Non-special domains with empty labels
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("h://.");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "h://."
            BOOST_TEST_CSTR_EQ(u.scheme(), "h");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), ".");
        }();
        // Broken IPv6
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://[www.google.com]/");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[google.com]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[::1.2.3.4x]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[::1.2.3.]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[::1.2.]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[::.1.2]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[::1.]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[::.1]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://[::%31]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://%5B::1]");
            BOOST_TEST_NOT(base && ref);
        }();
        // Misc Unicode
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://foo:💩@example.com/bar");
            // BOOST_TEST_NOT(ref);
        }();
        // resolving a fragment against any scheme succeeds
        []{
            system::result<url> base = parse_uri("test:test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "test:test#"
            BOOST_TEST_CSTR_EQ(u.scheme(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "test");
        }();
        []{
            system::result<url> base = parse_uri("mailto:x@x.com");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "mailto:x@x.com#x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "mailto");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "x@x.com");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "x");
        }();
        []{
            system::result<url> base = parse_uri("data:,");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "data:,#x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "data");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), ",");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "x");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "about:blank#x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "about");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "blank");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "x");
        }();
        []{
            system::result<url> base = parse_uri("test:test?test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "test:test?test#"
            BOOST_TEST_CSTR_EQ(u.scheme(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "test");
        }();
        // multiple @ in authority state
        []{
            system::result<url> base = parse_uri("http://doesnotmatter/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support multiple "@" in the authority
            // Ada and Node.js would consider the last "@" as the separator
            system::result<url> ref = parse_uri_reference("https://@test@test@example:800/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://doesnotmatter/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support multiple "@" in the authority
            // Ada and Node.js would consider the last "@" as the separator
            system::result<url> ref = parse_uri_reference("https://@@@example");
            BOOST_TEST_NOT(ref);
        }();
        // non-az-09 characters
        []{
            system::result<url> base = parse_uri("http://doesnotmatter/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid authority char "`"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://`{}:`{}@h/`{}?`{}");
            BOOST_TEST_NOT(ref);
        }();
        // byte is ' and url is special
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://host/?'");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://host/?%27"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "'");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("notspecial://host/?'");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "notspecial://host/?'"
            BOOST_TEST_CSTR_EQ(u.scheme(), "notspecial");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "'");
        }();
        // Credentials in base
        []{
            system::result<url> base = parse_uri("http://user@example.org/smth");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/some/path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://user@example.org/some/path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "user");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/some/path");
        }();
        []{
            system::result<url> base = parse_uri("http://user:pass@example.org:21/smth");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not fail on empty references
            system::result<url> ref = parse_uri_reference("");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://user:pass@example.org:21/smth");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/some/path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://user:pass@example.org:21/some/path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "user");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.port(), "21");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/some/path");
        }();
        // a set of tests designed by zcorpan for relative URLs with unknown schemes
        []{
            system::result<url> base = parse_uri("sc:/pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:/pa/i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/pa/i");
        }();
        []{
            system::result<url> base = parse_uri("sc://ho/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc://ho/i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "ho");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/i");
        }();
        []{
            system::result<url> base = parse_uri("sc:///pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:///pa/i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/pa/i");
        }();
        []{
            system::result<url> base = parse_uri("sc:/pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("../i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:/i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/i");
        }();
        []{
            system::result<url> base = parse_uri("sc://ho/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("../i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc://ho/i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "ho");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', 'i']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/i");
        }();
        []{
            system::result<url> base = parse_uri("sc:///pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("../i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:///i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/i");
        }();
        []{
            system::result<url> base = parse_uri("sc:/pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:/i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/i");
        }();
        []{
            system::result<url> base = parse_uri("sc://ho/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc://ho/i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "ho");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/i");
        }();
        []{
            system::result<url> base = parse_uri("sc:///pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:///i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/i");
        }();
        []{
            system::result<url> base = parse_uri("sc:/pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("?i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:/pa/pa?i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/pa/pa");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "i");
        }();
        []{
            system::result<url> base = parse_uri("sc://ho/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("?i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc://ho/pa?i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "ho");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/pa");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "i");
        }();
        []{
            system::result<url> base = parse_uri("sc:///pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("?i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:///pa/pa?i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/pa/pa");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "i");
        }();
        []{
            system::result<url> base = parse_uri("sc:sd");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:sd#i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "sd");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "i");
        }();
        []{
            system::result<url> base = parse_uri("sc:sd/sd");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:sd/sd#i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "sd/sd");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "i");
        }();
        []{
            system::result<url> base = parse_uri("sc:/pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:/pa/pa#i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/pa/pa");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "i");
        }();
        []{
            system::result<url> base = parse_uri("sc://ho/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc://ho/pa#i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "ho");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/pa");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "i");
        }();
        []{
            system::result<url> base = parse_uri("sc:///pa/pa");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#i");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:///pa/pa#i"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/pa/pa");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "i");
        }();
        // make sure that relative URL logic works on known typically non-relative schemes too
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("about:/../");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "about:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "about");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("data:/../");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "data:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "data");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("javascript:/../");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "javascript:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "javascript");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("mailto:/../");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "mailto:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "mailto");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        // unknown schemes and their hosts
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("sc://ñ.test/");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("sc://%/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("sc://@/");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support multiple "@" in the authority
            // Ada and Node.js would consider the last "@" as the separator
            system::result<url> ref = parse_uri_reference("sc://te@s:t@/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("sc://:/");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("sc://:12/");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            // Boost.URL does not support unicode
            // system::result<url> base = parse_uri("sc://ñ");
            // BOOST_TEST_NOT(base);
        }();
        // unknown schemes and backslashes
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("sc:\\../");
            BOOST_TEST_NOT(ref);
        }();
        // unknown scheme with path looking like a password
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("sc::a@example.net");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc::a@example.net"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), ":a@example.net");
        }();
        // unknown scheme with bogus percent-encoding
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("wow:%NBD");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("wow:%1G");
            BOOST_TEST_NOT(ref);
        }();
        // unknown scheme with non-URL characters
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("wow:\uffff");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://example.com/\ud800\U000107fe\udfff\ufdd0\ufdcf\ufdefﷰ\ufffe\uffff?\ud800\U000107fe\udfff\ufdd0\ufdcf\ufdefﷰ\ufffe\uffff");
            // BOOST_TEST_NOT(ref);
        }();
        // Forbidden host code points
        []{
            system::result<url> base = parse_uri("about:blank");
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("sc://a\x00" "b/");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("sc://a b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char "<"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("sc://a<b");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char ">"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("sc://a>b");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("sc://a[b/");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("sc://a\\b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("sc://a]b/");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char "^"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("sc://a^b");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char "|"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("sc://a|b/");
            BOOST_TEST_NOT(ref);
        }();
        // Forbidden host codepoints: tabs and newlines are removed during preprocessing
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("foo://ho\tst/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("foo://ho\nst/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("foo://ho\rst/");
            BOOST_TEST_NOT(ref);
        }();
        // Forbidden domain code-points
        []{
            system::result<url> base = parse_uri("about:blank");
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("http://a\x00" "b/");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x01
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x01" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x02
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x02" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x03
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x03" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x04
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x04" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x05
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x05" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x06
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x06" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://a\x07" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://a\x08" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://a\x0b" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://a\x0c" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x0e
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x0e" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x0f
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x0f" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x10
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x10" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x11
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x11" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x12
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x12" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x13
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x13" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x14
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x14" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x15
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x15" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x16
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x16" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x17
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x17" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x18
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x18" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x19
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x19" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x1a
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x1a" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x1b
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x1b" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x1c
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x1c" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x1d
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x1d" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x1e
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x1e" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x1f
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x1f" "b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("http://a b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://a%b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char "<"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a<b");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char ">"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a>b");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://a[b/");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://a]b/");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char "^"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a^b");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char "|"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a|b/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char \x7f
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://a\x7f" "b/");
            BOOST_TEST_NOT(ref);
        }();
        // Forbidden domain codepoints: tabs and newlines are removed during preprocessing
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://ho\tst/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://ho\nst/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://ho\rst/");
            BOOST_TEST_NOT(ref);
        }();
        // Encoded forbidden domain codepoints in special URLs
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%00st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%01st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%02st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%03st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%04st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%05st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%06st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%07st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%08st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%09st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%0Ast/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%0Bst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%0Cst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%0Dst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%0Est/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%0Fst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%10st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%11st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%12st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%13st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%14st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%15st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%16st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%17st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%18st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%19st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%1Ast/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%1Bst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%1Cst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%1Dst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%1Est/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%1Fst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%20st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%23st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%25st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%2Fst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%3Ast/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%3Cst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%3Est/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%3Fst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%40st/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%5Bst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%5Cst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%5Dst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%7Cst/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("http://ho%7Fst/");
            BOOST_TEST(ref);
        }();
        // Allowed host/domain code points
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid authority char """
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://!\"$&\'()*+,-.;=_`{}~/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("sc://\x01" "\x02" "\x03" "\x04" "\x05" "\x06" "\x07" "\x08" "\x0b" "\x0c" "\x0e" "\x0f" "\x10" "\x11" "\x12" "\x13" "\x14" "\x15" "\x16" "\x17" "\x18" "\x19" "\x1a" "\x1b" "\x1c" "\x1d" "\x1e" "\x1f" "\x7f" "!\"$%&\'()*+,-.;=_`{}~/");
            BOOST_TEST_NOT(ref);
        }();
        // Hosts and percent-encoding
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("ftp://example.com%80/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("ftp://example.com%A0/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("https://example.com%80/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("https://example.com%A0/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("ftp://%e2%98%83");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("https://%e2%98%83");
            BOOST_TEST(ref);
        }();
        // tests from jsdom/whatwg-url designed for code coverage
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://127.0.0.1:10100/relative_import.html");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://127.0.0.1:10100/relative_import.html"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.port(), "10100");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/relative_import.html");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://facebook.com/?foo=%7B%22abc%22");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://facebook.com/?foo=%7B%22abc%22"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "facebook.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "foo=%7B%22abc%22");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https://localhost:3000/jqueryui@1.2.3");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "https://localhost:3000/jqueryui@1.2.3"
            BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "localhost");
            BOOST_TEST_CSTR_EQ(u.port(), "3000");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/jqueryui@1.2.3");
        }();
        // tab/LF/CR
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("h\tt\nt\rp://h\to\ns\rt:9\t0\n0\r0/p\ta\nt\rh?q\tu\ne\rry#f\tr\na\rg");
            BOOST_TEST_NOT(ref);
        }();
        // Stringification of URL.searchParams
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("?a=b&c=d");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/bar?a=b&c=d"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "a=b&c=d");
        }();
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("??a=b&c=d");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/bar??a=b&c=d"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "?a=b&c=d");
        }();
        // Scheme only
        []{
            system::result<url> base = parse_uri("http://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/foo/bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar");
        }();
        []{
            system::result<url> base = parse_uri("https://example.org/foo/bar");
            system::result<url> ref = parse_uri_reference("http:");
            // Boost.URL does validate individual schemes
            // Ada and Node.js fail on "http:" with no authority
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("https://example.org/foo/bar");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("sc:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
        }();
        // Percent encoding of fragments
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("http://foo.bar/baz?qux#foo\x08" "bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid fragment char """
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://foo.bar/baz?qux#foo\"bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid fragment char "<"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://foo.bar/baz?qux#foo<bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid fragment char ">"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://foo.bar/baz?qux#foo>bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid fragment char "`"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://foo.bar/baz?qux#foo`bar");
            BOOST_TEST_NOT(ref);
        }();
        // IPv4 parsing (via https://github.com/nodejs/node/pull/10317)
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://1.2.3.4/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://1.2.3.4/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://1.2.3.4./");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://1.2.3.4/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://192.168.257");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://192.168.1.1"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "192.168.1.1"
            BOOST_TEST_NOT(u.encoded_host() == "192.168.1.1");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://192.168.257.");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://192.168.1.1"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "192.168.1.1"
            BOOST_TEST_NOT(u.encoded_host() == "192.168.1.1");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://192.168.257.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://192.168.257.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "192.168.257.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://256");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://0.0.1.0"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "0.0.1.0"
            BOOST_TEST_NOT(u.encoded_host() == "0.0.1.0");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://256.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://256.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "256.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://999999999");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://59.154.201.255"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "59.154.201.255"
            BOOST_TEST_NOT(u.encoded_host() == "59.154.201.255");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://999999999.");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://59.154.201.255"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "59.154.201.255"
            BOOST_TEST_NOT(u.encoded_host() == "59.154.201.255");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://999999999.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://999999999.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "999999999.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            system::result<url> ref = parse_uri_reference("http://10000000000");
            // Boost.URL does not fail on port overflow
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://10000000000.com");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://10000000000.com"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "10000000000.com");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://4294967295");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://255.255.255.255"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "255.255.255.255"
            BOOST_TEST_NOT(u.encoded_host() == "255.255.255.255");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL does support invalid IDNA
            system::result<url> ref = parse_uri_reference("http://4294967296");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://0xffffffff");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://255.255.255.255"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "255.255.255.255"
            BOOST_TEST_NOT(u.encoded_host() == "255.255.255.255");
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://0xffffffff1");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://256.256.256.256");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("https://0x.0x.0");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "https://0.0.0.0"
            BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "0.0.0.0"
            BOOST_TEST_NOT(u.encoded_host() == "0.0.0.0");
        }();
        // More IPv4 parsing (via https://github.com/jsdom/whatwg-url/issues/92)
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("https://0x100000000/test");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("https://256.0.0.1/test");
            BOOST_TEST(ref);
        }();
        // file URLs containing percent-encoded Windows drive letters (shouldn't work)
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///C%3A/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///C%3A/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL normalizes octets when resolving the URL
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C%3A/");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///C%7C/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///C%7C/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C%7C/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("file://%43%3A");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("file://%43%7C");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char "|"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("file://%43|");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("file://C%7C");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("file://%43%7C/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("https://%43%7C/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support invalid authority char "|"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("asdf://%43|/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does support invalid encoded chars
            system::result<url> ref = parse_uri_reference("asdf://%43%7C/");
            BOOST_TEST(ref);
        }();
        // file URLs relative to other file URLs (via https://github.com/jsdom/whatwg-url/pull/60)
        []{
            system::result<url> base = parse_uri("file:///C:/Users/Domenic/Dropbox/GitHub/tmpvar/jsdom/test/level2/html/files/anchor.html");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("pix/submit.gif");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///C:/Users/Domenic/Dropbox/GitHub/tmpvar/jsdom/test/level2/html/files/pix/submit.gif"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/Users/Domenic/Dropbox/GitHub/tmpvar/jsdom/test/level2/html/files/pix/submit.gif");
        }();
        []{
            system::result<url> base = parse_uri("file:///C:/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///C:/");
        }();
        []{
            system::result<url> base = parse_uri("file:///");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        // More file URL tests by zcorpan and annevk
        []{
            system::result<url> base = parse_uri("file:///C:/a/b");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///C:/");
        }();
        []{
            system::result<url> base = parse_uri("file://h/C:/a/b");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file://h/C:/");
        }();
        []{
            system::result<url> base = parse_uri("file://h/a/b");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://h/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "h");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("file:///C:/a/b");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//d:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///d:");
        }();
        []{
            system::result<url> base = parse_uri("file:///C:/a/b");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//d:/..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///d:/");
        }();
        []{
            system::result<url> base = parse_uri("file:///ab:/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("file:///1:/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("..");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("file:///test?test#test");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not fail on empty references
            system::result<url> ref = parse_uri_reference("");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///test?test#test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///test?test"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "test");
        }();
        []{
            system::result<url> base = parse_uri("file:///test?test#test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("?x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///test?x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "x");
        }();
        []{
            system::result<url> base = parse_uri("file:///test?test#test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:?x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///test?x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "x");
        }();
        []{
            system::result<url> base = parse_uri("file:///test?test#test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///test?test#x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "x");
        }();
        []{
            system::result<url> base = parse_uri("file:///test?test#test");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:#x");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///test?test#x"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "test");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "x");
        }();
        // File URLs and many (back)slashes
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("file:\\\\//");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("file:\\\\\\\\");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("file:\\\\\\\\?fox");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("file:\\\\\\\\#guppy");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://spider///");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://spider///"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "spider");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "///");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("file:\\\\localhost//");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///localhost//cat");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///localhost//cat"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/localhost//cat");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("file://\\/localhost//cat");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://localhost//a//../..//");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://///"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "///");
        }();
        []{
            system::result<url> base = parse_uri("file:///elephant");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/////mouse");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://///mouse"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "///mouse");
        }();
        []{
            system::result<url> base = parse_uri("file://lion/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("\\//pig");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://lion/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("\\/localhost//pig");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://lion/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//localhost//pig");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:////pig");
        }();
        []{
            system::result<url> base = parse_uri("file://lion/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/..//localhost//pig");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file://lion//localhost//pig");
        }();
        []{
            system::result<url> base = parse_uri("file://ape/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///");
        }();
        // File URLs with non-empty hosts
        []{
            system::result<url> base = parse_uri("file://tea/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/rooibos");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://tea/rooibos"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "tea");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/rooibos");
        }();
        []{
            system::result<url> base = parse_uri("file://tea/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/?chai");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://tea/?chai"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "tea");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "chai");
        }();
        // Windows drive letter handling with the 'file:' base URL
        []{
            system::result<url> base = parse_uri("file://host/dir/file");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("C|");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/D:/dir1/dir2/file");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("C|");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/dir/file");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("C|#");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/dir/file");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("C|?");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/dir/file");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("C|/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/dir/file");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("C|\n/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/dir/file");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("C|\\");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/dir/file");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("C");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://host/dir/C"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/dir/C");
        }();
        []{
            system::result<url> base = parse_uri("file://host/dir/file");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("C|a");
            BOOST_TEST_NOT(ref);
        }();
        // Windows drive letter quirk in the file slash state
        []{
            system::result<url> base = parse_uri("file:///c:/baz/qux");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/c:/foo/bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///c:/foo/bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/c:/foo/bar");
        }();
        []{
            system::result<url> base = parse_uri("file:///c:/baz/qux");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("/c|/foo/bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file:///c:/baz/qux");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("file:\\c:\\foo\\bar");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/c:/foo/bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://host/c:/foo/bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/c:/foo/bar");
        }();
        // Do not drop the host in the presence of a drive letter
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://example.net/C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://example.net/C:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://1.2.3.4/C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://1.2.3.4/C:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://[1::8]/C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://[1::8]/C:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        // Copy the host from the base URL in the following cases
        []{
            system::result<url> base = parse_uri("file://host/");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("C|/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("file://host/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://host/C:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        []{
            system::result<url> base = parse_uri("file://host/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://host/C:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        []{
            system::result<url> base = parse_uri("file://host/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:/C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://host/C:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        // Copy the empty host from the input in the following cases
        []{
            system::result<url> base = parse_uri("file://host/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///C:/");
        }();
        []{
            system::result<url> base = parse_uri("file://host/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///C:/");
        }();
        []{
            system::result<url> base = parse_uri("file://host/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("///C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///C:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        []{
            system::result<url> base = parse_uri("file://host/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///C:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///C:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/C:/");
        }();
        // Windows drive letter quirk (no host)
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support "|" in paths
            system::result<url> ref = parse_uri_reference("file:/C|/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid authority char "|"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("file://C|/");
            BOOST_TEST_NOT(ref);
        }();
        // file URLs without base URL by Rimas Misevičius
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:?q=v");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///?q=v");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:#frag");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:///#frag");
        }();
        // file: drive letter cases from https://crbug.com/1078698
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///Y:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///Y:"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/Y:");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///Y:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///Y:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/Y:/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///./Y");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///Y"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/Y");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///./Y:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///Y:"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/Y:");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("\\\\\\.\\Y:");
            BOOST_TEST_NOT(ref);
        }();
        // file: drive letter cases from https://crbug.com/1078698 but lowercased
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///y:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///y:"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/y:");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///y:/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///y:/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/y:/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///./y");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///y"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/y");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///./y:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///y:"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/y:");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("\\\\\\.\\y:");
            BOOST_TEST_NOT(ref);
        }();
        // Additional file URL tests for (https://github.com/whatwg/url/issues/405)
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://localhost//a//../..//foo");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://///foo"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "///foo");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file://localhost////foo");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://////foo"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "////foo");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:////foo");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:////foo"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//foo");
        }();
        []{
            system::result<url> base = parse_uri("file:///");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///one/two");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///one/two"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/one/two");
        }();
        []{
            system::result<url> base = parse_uri("file:///");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:////one/two");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:////one/two"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//one/two");
        }();
        []{
            system::result<url> base = parse_uri("file:///");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("//one/two");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file://one/two"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "one");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/two");
        }();
        []{
            system::result<url> base = parse_uri("file:///");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("///one/two");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///one/two"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/one/two");
        }();
        []{
            system::result<url> base = parse_uri("file:///");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("////one/two");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:////one/two"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//one/two");
        }();
        []{
            system::result<url> base = parse_uri("file:////");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:///.//");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:////"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//");
        }();
        // File URL tests for https://github.com/whatwg/url/issues/549
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:.//p");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:////p");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("file:/.//p");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no exceptional logic for "file:"
            BOOST_TEST_NOT(u.buffer() == "file:////p");
        }();
        // IPv6 tests
        []{
            system::result<url> base = parse_uri("http://example.net/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://[1:0::]");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://[1::]"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "[1::]"
            BOOST_TEST_NOT(u.encoded_host() == "[1::]");
        }();
        []{
            system::result<url> base = parse_uri("http://example.net/");
            system::result<url> ref = parse_uri_reference("http://[0:1:2:3:4:5:6:7:8]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("https://[0::0::0]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("https://[0:.0]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("https://[0:0:]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("https://[0:1:2:3:4:5:6:7.0.0.0.1]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("https://[0:1.00.0.0.0]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("https://[0:1.290.0.0.0]");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("https://[0:1.23.23]");
            BOOST_TEST_NOT(base && ref);
        }();
        // Empty host
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://?");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://#");
            // Boost.URL does not fail on empty hostname
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        // Port overflow (2^32 + 81)
        []{
            system::result<url> base = parse_uri("http://example.org/");
            system::result<url> ref = parse_uri_reference("http://f:4294967377/c");
            // Boost.URL does not fail on port overflow
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        // Port overflow (2^64 + 81)
        []{
            system::result<url> base = parse_uri("http://example.org/");
            system::result<url> ref = parse_uri_reference("http://f:18446744073709551697/c");
            // Boost.URL does not fail on port overflow
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        // Port overflow (2^128 + 81)
        []{
            system::result<url> base = parse_uri("http://example.org/");
            system::result<url> ref = parse_uri_reference("http://f:340282366920938463463374607431768211537/c");
            // Boost.URL does not fail on port overflow
            BOOST_TEST(base);
            BOOST_TEST(ref);
        }();
        // Non-special-URL path tests
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("sc://ñ");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("sc://ñ?x");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("sc://ñ#x");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            // Boost.URL does not support unicode
            // system::result<url> base = parse_uri("sc://ñ");
            // BOOST_TEST_NOT(base);
        }();
        []{
            // Boost.URL does not support unicode
            // system::result<url> base = parse_uri("sc://ñ");
            // BOOST_TEST_NOT(base);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("sc://?");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc://?"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("sc://#");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc://#"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
        }();
        []{
            system::result<url> base = parse_uri("sc://x/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("///");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:///"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("sc://x/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("////");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:////"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//");
        }();
        []{
            system::result<url> base = parse_uri("sc://x/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("////x/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "sc:////x/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "sc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//x/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("tftp://foobar.com/someconfig;mode=netascii");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "tftp://foobar.com/someconfig;mode=netascii"
            BOOST_TEST_CSTR_EQ(u.scheme(), "tftp");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foobar.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/someconfig;mode=netascii");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("telnet://user:pass@foobar.com:23/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "telnet://user:pass@foobar.com:23/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "telnet");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "user");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "foobar.com");
            BOOST_TEST_CSTR_EQ(u.port(), "23");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ut2004://10.10.10.10:7777/Index.ut2");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ut2004://10.10.10.10:7777/Index.ut2"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ut2004");
            BOOST_TEST_CSTR_EQ(u.port(), "7777");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/Index.ut2");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("redis://foo:bar@somehost:6379/0?baz=bam&qux=baz");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "redis://foo:bar@somehost:6379/0?baz=bam&qux=baz"
            BOOST_TEST_CSTR_EQ(u.scheme(), "redis");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "somehost");
            BOOST_TEST_CSTR_EQ(u.port(), "6379");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/0");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "baz=bam&qux=baz");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("rsync://foo@host:911/sup");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "rsync://foo@host:911/sup"
            BOOST_TEST_CSTR_EQ(u.scheme(), "rsync");
            BOOST_TEST_CSTR_EQ(u.encoded_user(), "foo");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.port(), "911");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/sup");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("git://github.com/foo/bar.git");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "git://github.com/foo/bar.git"
            BOOST_TEST_CSTR_EQ(u.scheme(), "git");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "github.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar.git");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("irc://myserver.com:6999/channel?passwd");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "irc://myserver.com:6999/channel?passwd"
            BOOST_TEST_CSTR_EQ(u.scheme(), "irc");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "myserver.com");
            BOOST_TEST_CSTR_EQ(u.port(), "6999");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/channel");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "passwd");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("dns://fw.example.org:9999/foo.bar.org?type=TXT");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "dns://fw.example.org:9999/foo.bar.org?type=TXT"
            BOOST_TEST_CSTR_EQ(u.scheme(), "dns");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "fw.example.org");
            BOOST_TEST_CSTR_EQ(u.port(), "9999");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo.bar.org");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "type=TXT");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("ldap://localhost:389/ou=People,o=JNDITutorial");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "ldap://localhost:389/ou=People,o=JNDITutorial"
            BOOST_TEST_CSTR_EQ(u.scheme(), "ldap");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "localhost");
            BOOST_TEST_CSTR_EQ(u.port(), "389");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/ou=People,o=JNDITutorial");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("git+https://github.com/foo/bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "git+https://github.com/foo/bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "git+https");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "github.com");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo/bar");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("urn:ietf:rfc:2648");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "urn:ietf:rfc:2648"
            BOOST_TEST_CSTR_EQ(u.scheme(), "urn");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "ietf:rfc:2648");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("tag:joe@example.org,2001:foo/bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "tag:joe@example.org,2001:foo/bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "tag");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "joe@example.org,2001:foo/bar");
        }();
        // Serialize /. in path
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-spec:/.//");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            // Boost.URL will not let "//" path become the authority separator
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-spec:/..//");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '', '']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-spec:/a/..//");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-spec:/.//path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            // Boost.URL will not let "//" path become the authority separator
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//path");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-spec:/..//path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '', 'path']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//path");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-spec:/a/..//path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//path");
        }();
        []{
            system::result<url> base = parse_uri("non-spec:/p");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/.//path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            // Boost.URL will not let "//" path become the authority separator
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//path");
        }();
        []{
            system::result<url> base = parse_uri("non-spec:/p");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("/..//path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '', 'path']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//path");
        }();
        []{
            system::result<url> base = parse_uri("non-spec:/p");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("..//path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '', 'path']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//path");
        }();
        []{
            system::result<url> base = parse_uri("non-spec:/p");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("a/..//path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "//path");
        }();
        []{
            system::result<url> base = parse_uri("non-spec:/..//p");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not fail on empty references
            system::result<url> ref = parse_uri_reference("");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("non-spec:/..//p");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/.//path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            // Boost.URL normalization does not remove ".." segments above root
            // Resolved path has segments ['..', '', 'path']
            // Removing ".." segments would make Errata 4547 innocuous
            // see https://www.rfc-editor.org/errata/eid4547
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//path");
        }();
        // Do not serialize /. in path
        []{
            system::result<url> base = parse_uri("non-spec:/.//p");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("../path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-spec:/path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-spec");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/path");
        }();
        // percent encoded hosts in non-special-URLs
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("non-special://%E2%80%A0/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-special://H%4fSt/path");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-special://H%4fSt/path"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-special");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "host");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/path");
        }();
        // IPv6 in non-special-URLs
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-special://[1:2:0:0:5:0:0:0]/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-special://[1:2:0:0:5::]/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-special");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "[1:2:0:0:5::]"
            BOOST_TEST_NOT(u.encoded_host() == "[1:2:0:0:5::]");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-special://[1:2:0:0:0:0:0:3]/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-special://[1:2::3]/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-special");
            // Boost.URL does not modify IP addresses
            // Expected hostname in Ada: "[1:2::3]"
            BOOST_TEST_NOT(u.encoded_host() == "[1:2::3]");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("non-special://[1:2::3]:80/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "non-special://[1:2::3]:80/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "non-special");
            BOOST_TEST_CSTR_EQ(u.port(), "80");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("non-special://[:80/");
            BOOST_TEST_NOT(base && ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("blob:https://example.com:443/");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "blob:https://example.com:443/"
            BOOST_TEST_CSTR_EQ(u.scheme(), "blob");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "https://example.com:443/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("blob:d3958f5c-0777-0845-9dcf-2cb28783acaf");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "blob:d3958f5c-0777-0845-9dcf-2cb28783acaf"
            BOOST_TEST_CSTR_EQ(u.scheme(), "blob");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "d3958f5c-0777-0845-9dcf-2cb28783acaf");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("blob:");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "blob:"
            BOOST_TEST_CSTR_EQ(u.scheme(), "blob");
        }();
        // Invalid IPv4 radix digits
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://0x7f.0.0.0x7g");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://0x7f.0.0.0x7g"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "0x7f.0.0.0x7g");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://0X7F.0.0.0X7G");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://0x7f.0.0.0x7g"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "0x7f.0.0.0x7g");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        // Invalid IPv4 portion of IPv6 address
        []{
            system::result<url> base = parse_uri("about:blank");
            system::result<url> ref = parse_uri_reference("http://[::127.0.0.0.1]");
            BOOST_TEST_NOT(base && ref);
        }();
        // Uncompressed IPv6 addresses with 0
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://[0:1:0:1:0:1:0:1]");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://[0:1:0:1:0:1:0:1]"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://[1:0:1:0:1:0:1:0]");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://[1:0:1:0:1:0:1:0]"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            // Boost.URL does not create absolute paths
            // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
        }();
        // Percent-encoded query and fragment
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid query char """
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://example.org/test?\"");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.org/test?#");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/test?#"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid query char "<"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://example.org/test?<");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid query char ">"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("http://example.org/test?>");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://example.org/test?⌣");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("http://example.org/test?%23%23");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "http://example.org/test?%23%23"
            BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
            BOOST_TEST_CSTR_EQ(u.encoded_query(), "%23%23");
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://example.org/test?%GH");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("http://example.org/test?a#%EF");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("http://example.org/test?a#%GH");
            BOOST_TEST_NOT(ref);
        }();
        // URLs that require a non-about:blank base. (Also serve as invalid base tests.)
        // Bases that don't fail to parse but fail to be bases
        // Other base URL tests, that must succeed
        []{
            system::result<url> base = parse_uri("a:/");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("test-a-colon-slash.html");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "a:/test-a-colon-slash.html"
            BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test-a-colon-slash.html");
        }();
        []{
            system::result<url> base = parse_uri("a://");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("test-a-colon-slash-slash.html");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "a:///test-a-colon-slash-slash.html"
            BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test-a-colon-slash-slash.html");
        }();
        []{
            system::result<url> base = parse_uri("a:/b");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("test-a-colon-slash-b.html");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "a:/test-a-colon-slash-b.html"
            BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test-a-colon-slash-b.html");
        }();
        []{
            system::result<url> base = parse_uri("a://b");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("test-a-colon-slash-slash-b.html");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "a://b/test-a-colon-slash-slash-b.html"
            BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "b");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test-a-colon-slash-slash-b.html");
        }();
        // Null code point in fragment
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("http://example.org/test?a#b\x00" "c");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("non-spec://example.org/test?a#b\x00" "c");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("non-spec:/test?a#b\x00" "c");
            // BOOST_TEST_NOT(ref);
        }();
        // First scheme char - not allowed: https://github.com/whatwg/url/issues/464
        []{
            system::result<url> base = parse_uri("file:///some/dir/bar.html");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("10.0.0.7:8080/foo.html");
            // Boost.URL requires segment-nz-nc
            // A URL should have a scheme or the first path segment
            // cannot contain a colon (":")
            BOOST_TEST_NOT(ref);
        }();
        // Subsequent scheme chars - not allowed
        []{
            system::result<url> base = parse_uri("file:///some/dir/bar.html");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("a!@$*=/foo.html");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "file:///some/dir/a!@$*=/foo.html"
            BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/some/dir/a!@$*=/foo.html");
        }();
        // First and subsequent scheme chars - allowed
        []{
            system::result<url> base = parse_uri("http://example.com/dir/file");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("a1234567890-+.:foo/bar");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "a1234567890-+.:foo/bar"
            BOOST_TEST_CSTR_EQ(u.scheme(), "a1234567890-+.");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "foo/bar");
        }();
        // IDNA ignored code points in file URLs hosts
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("file://a\xad" "b/p");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("file://a%C2%ADb/p");
            BOOST_TEST(ref);
        }();
        // IDNA hostnames which get mapped to 'localhost'
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("file://loC𝐀𝐋𝐇𝐨𝐬𝐭/usr/bin");
            // BOOST_TEST_NOT(ref);
        }();
        // Empty host after the domain to ASCII
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("file://\xad" "/p");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("file://%C2%AD/p");
            BOOST_TEST(ref);
        }();
        // https://bugzilla.mozilla.org/show_bug.cgi?id=1647058
        []{
            system::result<url> base = parse_uri("https://example.org/##link");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("#link");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "https://example.org/#link"
            BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.org");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
            BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "link");
        }();
        // UTF-8 percent-encode of C0 control percent-encode set and supersets
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("non-special:cannot-be-a-base-url-\x00" "\x01" "\x1f" "\x1e" "~\x7f" "\x80" "");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid path char "{"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("https://www.example.com/path{\x7f" "path.html?query'\x7f" "=query#fragment<\x7f" "fragment");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://example.org");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid query char "\x7f"
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("https://user:pass[\x7f" "@foo/bar");
            BOOST_TEST_NOT(ref);
        }();
        // Tests for the distinct percent-encode sets
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("foo:// !\"$%&\'()*+,-.;<=>@[\\]^_`{|}~@host/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("wss:// !\"$%&\'()*+,-.;<=>@[]^_`{|}~@host/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("foo://joe: !\"$%&\'()*+,-.:;<=>@[\\]^_`{|}~@host/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support whitespace in input
            system::result<url> ref = parse_uri_reference("wss://joe: !\"$%&\'()*+,-.:;<=>@[]^_`{|}~@host/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support incomplete octets
            system::result<url> ref = parse_uri_reference("foo://!\"$%&\'()*+,-.;=_`{}~/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid authority char """
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("wss://!\"$&\'()*+,-.;=_`{}~/");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("foo://host/ !\"$%&\'()*+,-./:;<=>@[\\]^_`{|}~");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("wss://host/ !\"$%&\'()*+,-./:;<=>@[\\]^_`{|}~");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("foo://host/dir/? !\"$%&\'()*+,-./:;<=>?@[\\]^_`{|}~");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("wss://host/dir/? !\"$%&\'()*+,-./:;<=>?@[\\]^_`{|}~");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("foo://host/dir/# !\"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support escape chars in input
            system::result<url> ref = parse_uri_reference("wss://host/dir/# !\"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~");
            BOOST_TEST_NOT(ref);
        }();
        // Ensure that input schemes are not ignored when resolving non-special URLs
        []{
            system::result<url> base = parse_uri("abc://host/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("abc:rootless");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no special logic for schemes
            BOOST_TEST_NOT(u.buffer() == "abc:rootless");
        }();
        []{
            system::result<url> base = parse_uri("abc:/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("abc:rootless");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no special logic for schemes
            BOOST_TEST_NOT(u.buffer() == "abc:rootless");
        }();
        []{
            system::result<url> base = parse_uri("abc:path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("abc:rootless");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // href: "abc:rootless"
            BOOST_TEST_CSTR_EQ(u.scheme(), "abc");
            BOOST_TEST_CSTR_EQ(u.encoded_path(), "rootless");
        }();
        []{
            system::result<url> base = parse_uri("abc://host/path");
            if (!BOOST_TEST(base))
                return;
            system::result<url> ref = parse_uri_reference("abc:/rooted");
            if (!BOOST_TEST(ref))
                return;
            url u;
            system::result<void> r = resolve(*base, *ref, u);
            BOOST_TEST(r);
            u.normalize_authority();
            u.normalize_path();
            // Boost.URL has no special logic for schemes
            BOOST_TEST_NOT(u.buffer() == "abc:/rooted");
        }();
        // Empty query and fragment with blank should throw an error
        []{
            system::result<url> r = parse_uri_reference("#");
            url u = *r;
        }();
        []{
            system::result<url> r = parse_uri_reference("?");
            url u = *r;
        }();
        // Last component looks like a number, but not valid IPv4
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://1.2.3.4.5");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://1.2.3.4.5.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://0..0x300/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://0..0x300./");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://256.256.256.256.256");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("http://other.com/");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://256.256.256.256.256.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://1.2.3.08");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://1.2.3.08.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://1.2.3.09");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://09.2.3.4");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://09.2.3.4.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://01.2.3.4.5");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://01.2.3.4.5.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://0x100.2.3.4");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://0x100.2.3.4.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://0x1.2.3.4.5");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://0x1.2.3.4.5.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://foo.1.2.3.4");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://foo.1.2.3.4.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://foo.2.3.4");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://foo.2.3.4.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://foo.09");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://foo.09.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://foo.0x4");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL can parse an overflown IPv4address as a valid reg-name
            // All chars in the invalid IPv4address are unreserved
            system::result<url> ref = parse_uri_reference("http://foo.0x4.");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not fail when a component looks like a number
            // All unreserved chars are accepted
            system::result<url> ref = parse_uri_reference("http://foo.09..");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not fail when a component looks like a number
            // All unreserved chars are accepted
            system::result<url> ref = parse_uri_reference("http://0999999999999999999/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not fail when a component looks like a number
            // All unreserved chars are accepted
            system::result<url> ref = parse_uri_reference("http://foo.0x");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not fail when a component looks like a number
            // All unreserved chars are accepted
            system::result<url> ref = parse_uri_reference("http://foo.0XFfFfFfFfFfFfFfFfFfAcE123");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("http://💩.123/");
            // BOOST_TEST_NOT(ref);
        }();
        // U+0000 and U+FFFF in various places
        []{
            system::result<url> base = parse_uri("about:blank");
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("https://\x00" "y");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("https://x/\x00" "y");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("https://x/?\x00" "y");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("https://x/?#\x00" "y");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("https://\uffffy");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("https://x/\uffffy");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("https://x/?\uffffy");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("https://x/?#\uffffy");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("non-special:\x00" "y");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("non-special:x/\x00" "y");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("non-special:x/?\x00" "y");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Invalid input: contains "\0"
            // system::result<url> ref = parse_uri_reference("non-special:x/?#\x00" "y");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("non-special:\uffffy");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("non-special:x/\uffffy");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("non-special:x/?\uffffy");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("non-special:x/?#\uffffy");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not fail on empty references
            system::result<url> ref = parse_uri_reference("");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does not support invalid path char """
            // Node.js and Ada would attempt to encode these chars
            system::result<url> ref = parse_uri_reference("https://example.com/\"quoted\"");
            BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            if (!BOOST_TEST(base))
                return;
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("https://a%C2%ADb/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does not support unicode input
            // system::result<url> ref = parse_uri_reference("https://\xad" "/");
            // BOOST_TEST_NOT(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid encoded unicode
            system::result<url> ref = parse_uri_reference("https://%C2%AD/");
            BOOST_TEST(ref);
        }();
        []{
            system::result<url> base = parse_uri("about:blank");
            // Boost.URL does support invalid IDNA
            system::result<url> ref = parse_uri_reference("https://xn--/");
            BOOST_TEST(ref);
        }();
    }

    static
    void
    adaSettersTests()
    {
        //  Tests for setters of https://url.spec.whatwg.org/#urlutils-members
        //
        // This file contains a JSON object.
        // Other than 'comment', each key is an attribute of the `URL` interface
        // defined in WHATWG’s URL Standard.
        // The values are arrays of test case objects for that attribute.
        //
        // To run a test case for the attribute `attr`:
        //
        // * Create a new `URL` object with the value for the 'href' key
        //   the constructor single parameter. (Without a base URL.)
        //   This must not throw.
        // * Set the attribute `attr` to (invoke its setter with)
        //   with the value of for 'new_value' key.
        // * The value for the 'expected' key is another object.
        //   For each `key` / `value` pair of that object,
        //   get the attribute `key` (invoke its getter).
        //   The returned string must be equal to `value`.
        //
        // Note: the 'href' setter is already covered by urltestdata.json.
        // protocol
        {
            // The empty string is not a valid scheme. Setter leaves the URL unchanged.
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // BOOST_TEST_NO_THROW(u.set_scheme(""));
                // BOOST_TEST_CSTR_EQ(u, "a://example.net");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "a");
                u.remove_scheme();
                BOOST_TEST_CSTR_EQ(u, "//example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "");
            }
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("b"));
                BOOST_TEST_CSTR_EQ(u, "b://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "b");
            }
            {
                system::result<url> r = parse_uri_reference("javascript:alert(1)");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("defuse"));
                BOOST_TEST_CSTR_EQ(u, "defuse:alert(1)");
                BOOST_TEST_CSTR_EQ(u.scheme(), "defuse");
            }
            // Upper-case ASCII is lower-cased
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("B"));
                // BOOST_TEST_CSTR_EQ(u, "b://example.net");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "b");
                // Although schemes are case-insensitive, the canonical form
                // is lowercase and documents that specify schemes must do so
                // with lowercase letters
                BOOST_TEST_CSTR_EQ(u, "B://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "B");
                u.normalize_scheme();
                BOOST_TEST_CSTR_EQ(u, "b://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "b");
            }
            // Non-ASCII is rejected
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("é"), system::system_error);
                BOOST_TEST_CSTR_EQ(u, "a://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            }
            // No leading digit
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("0b"), system::system_error);
                BOOST_TEST_CSTR_EQ(u, "a://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            }
            // No leading punctuation
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("+b"), system::system_error);
                BOOST_TEST_CSTR_EQ(u, "a://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            }
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("bC0+-."));
                // BOOST_TEST_CSTR_EQ(u, "bc0+-.://example.net");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "bc0+-.");
                BOOST_TEST_CSTR_EQ(u, "bC0+-.://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "bC0+-.");
                u.normalize_scheme();
                BOOST_TEST_CSTR_EQ(u, "bc0+-.://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "bc0+-.");
            }
            // Only some punctuation is acceptable
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("b,c"), system::system_error);
                BOOST_TEST_CSTR_EQ(u, "a://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            }
            // Non-ASCII is rejected
            {
                system::result<url> r = parse_uri_reference("a://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("bé"), system::system_error);
                BOOST_TEST_CSTR_EQ(u, "a://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "a");
            }
            // Ada can’t switch from URL containing username/password/port to file
            // Boost.URL handles the general syntax defined in RFC3986
            // so there's no concept of a special scheme
            {
                system::result<url> r = parse_uri_reference("http://test@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("file"));
                // BOOST_TEST_CSTR_EQ(u, "http://test@example.net/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "http");
                BOOST_TEST_CSTR_EQ(u, "file://test@example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net:1234");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("file"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net:1234/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "https");
                BOOST_TEST_CSTR_EQ(u, "file://example.net:1234");
                BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            }
            {
                system::result<url> r = parse_uri_reference("wss://x:x@example.net:1234");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("file"));
                // BOOST_TEST_CSTR_EQ(u, "wss://x:x@example.net:1234/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "wss");
                BOOST_TEST_CSTR_EQ(u, "file://x:x@example.net:1234");
                BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            }
            // Ada can’t switch from file URL with no host
            // Boost.URL handles the general syntax defined in RFC3986
            // so there's no concept of a special scheme
            {
                system::result<url> r = parse_uri_reference("file://localhost/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("http"));
                // BOOST_TEST_CSTR_EQ(u, "file:///");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "file");
                BOOST_TEST_CSTR_EQ(u, "http://localhost/");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            }
            {
                system::result<url> r = parse_uri_reference("file:///test");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("https"));
                // BOOST_TEST_CSTR_EQ(u, "file:///test");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "file");
                BOOST_TEST_CSTR_EQ(u, "https:///test");
                BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            }
            {
                system::result<url> r = parse_uri_reference("file:");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("wss"));
                // BOOST_TEST_CSTR_EQ(u, "file:///");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "file");
                BOOST_TEST_CSTR_EQ(u, "wss:");
                BOOST_TEST_CSTR_EQ(u.scheme(), "wss");
            }
            // Ada can’t switch from special scheme to non-special
            // Boost.URL handles the general syntax defined in RFC3986
            // so there's no concept of a special scheme
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("b"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "http");
                BOOST_TEST_CSTR_EQ(u, "b://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "b");
            }
            {
                system::result<url> r = parse_uri_reference("file://hi/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("s"));
                // BOOST_TEST_CSTR_EQ(u, "file://hi/path");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "file");
                BOOST_TEST_CSTR_EQ(u, "s://hi/path");
                BOOST_TEST_CSTR_EQ(u.scheme(), "s");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("s"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "https");
                BOOST_TEST_CSTR_EQ(u, "s://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "s");
            }
            {
                system::result<url> r = parse_uri_reference("ftp://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("test"));
                // BOOST_TEST_CSTR_EQ(u, "ftp://example.net/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "ftp");
                BOOST_TEST_CSTR_EQ(u, "test://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "test");
            }
            // Cannot-be-a-base URL doesn’t have a host, but URL in a special scheme must.
            {
                system::result<url> r = parse_uri_reference("mailto:me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("http"));
                // BOOST_TEST_CSTR_EQ(u, "mailto:me@example.net");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "mailto");
                BOOST_TEST_CSTR_EQ(u, "http:me@example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            }
            // Ada can’t switch from non-special scheme to special
            // Boost.URL handles the general syntax defined in RFC3986
            // so there's no concept of a special scheme
            {
                system::result<url> r = parse_uri_reference("ssh://me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("http"));
                // BOOST_TEST_CSTR_EQ(u, "ssh://me@example.net");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "ssh");
                BOOST_TEST_CSTR_EQ(u, "http://me@example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            }
            {
                system::result<url> r = parse_uri_reference("ssh://me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("https"));
                // BOOST_TEST_CSTR_EQ(u, "ssh://me@example.net");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "ssh");
                BOOST_TEST_CSTR_EQ(u, "https://me@example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            }
            {
                system::result<url> r = parse_uri_reference("ssh://me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("file"));
                // BOOST_TEST_CSTR_EQ(u, "ssh://me@example.net");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "ssh");
                BOOST_TEST_CSTR_EQ(u, "file://me@example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            }
            {
                system::result<url> r = parse_uri_reference("ssh://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("file"));
                // BOOST_TEST_CSTR_EQ(u, "ssh://example.net");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "ssh");
                BOOST_TEST_CSTR_EQ(u, "file://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "file");
            }
            {
                system::result<url> r = parse_uri_reference("nonsense:///test");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("https"));
                // BOOST_TEST_CSTR_EQ(u, "nonsense:///test");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "nonsense");
                BOOST_TEST_CSTR_EQ(u, "https:///test");
                BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            }
            // Ada ignores substrings after the first ':'
            // Boost.URL throws on invalid schemes
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("https:foo : bar"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "https");
                BOOST_TEST_CSTR_EQ(u, "http://example.net");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            }
            // Ada ignores substrings after the first ':'
            // Boost.URL throws on invalid schemes
            {
                // Ada accepts literal `<>`s in path
                // Boost.URL only accepts segment-nz = 1*pchar
                system::result<url> r = parse_uri_reference("data:text/html,<p>Test");
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_scheme("view-source+data:foo : bar"));
                // BOOST_TEST_CSTR_EQ(u, "view-source+data:text/html,<p>Test");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "view-source+data");
            }
            // Ada sets port to null if it is the default for new scheme.
            // Boost.URL keeps the URL the way the user defined it
            {
                system::result<url> r = parse_uri_reference("http://foo.com:443/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("https"));
                // BOOST_TEST_CSTR_EQ(u, "https://foo.com/");
                BOOST_TEST_CSTR_EQ(u, "https://foo.com:443/");
                BOOST_TEST_CSTR_EQ(u.scheme(), "https");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u.port(), "443");
            }
            // Ada strips tab and newline
            // Boost.URL does not accept invalid chars
            {
                system::result<url> r = parse_uri_reference("http://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("h\r\ntt\tps"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "https://test/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "https");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "http://test/");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            {
                system::result<url> r = parse_uri_reference("http://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("https\r"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "https://test/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "https");
                BOOST_TEST_CSTR_EQ(u, "http://test/");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            }
            // Non-tab/newline C0 controls result in no-op
            // Boost.URL core::string_view constructor parses until it finds null
            {
                system::result<url> r = parse_uri_reference("http://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_scheme("https\x00" ""));
                // BOOST_TEST_CSTR_EQ(u, "http://test/");
                // BOOST_TEST_CSTR_EQ(u.scheme(), "http");
                BOOST_TEST_CSTR_EQ(u, "https://test/");
                BOOST_TEST_CSTR_EQ(u.scheme(), "https");
            }
            {
                system::result<url> r = parse_uri_reference("http://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("https\x0c" ""), system::system_error);
                BOOST_TEST_CSTR_EQ(u, "http://test/");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            }
            {
                system::result<url> r = parse_uri_reference("http://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("https\x0e" ""), system::system_error);
                BOOST_TEST_CSTR_EQ(u, "http://test/");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            }
            {
                system::result<url> r = parse_uri_reference("http://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_scheme("https "), system::system_error);
                BOOST_TEST_CSTR_EQ(u, "http://test/");
                BOOST_TEST_CSTR_EQ(u.scheme(), "http");
            }
        }
        // username
        {
            // In Ada, no host means no username
            // Boost.URL creates the host when the username is set
            {
                system::result<url> r = parse_uri_reference("file:///home/you/index.html");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("me"));
                // BOOST_TEST_CSTR_EQ(u, "file:///home/you/index.html");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
                BOOST_TEST_CSTR_EQ(u, "file://me@/home/you/index.html");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "me");
            }
            // In Ada, no host means no username
            // Boost.URL creates the host when the username is set
            {
                system::result<url> r = parse_uri_reference("unix:/run/foo.socket");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("me"));
                // BOOST_TEST_CSTR_EQ(u, "unix:/run/foo.socket");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
                BOOST_TEST_CSTR_EQ(u, "unix://me@/run/foo.socket");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "me");
            }
            // In Ada, cannot-be-a-base (no host) means no username
            // Boost.URL still replaces the username
            {
                system::result<url> r = parse_uri_reference("mailto:you@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("me"));
                // BOOST_TEST_CSTR_EQ(u, "mailto:you@example.net");
                BOOST_TEST_CSTR_EQ(u, "mailto://me@/you@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "me");
            }
            {
                system::result<url> r = parse_uri_reference("javascript:alert(1)");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("wario"));
                // BOOST_TEST_CSTR_EQ(u, "javascript:alert(1)");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
                BOOST_TEST_CSTR_EQ(u, "javascript://wario@/alert(1)");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "wario");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("me"));
                // BOOST_TEST_CSTR_EQ(u, "http://me@example.net/");
                BOOST_TEST_CSTR_EQ(u, "http://me@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "me");
            }
            {
                system::result<url> r = parse_uri_reference("http://:secret@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("me"));
                // BOOST_TEST_CSTR_EQ(u, "http://me:secret@example.net/");
                BOOST_TEST_CSTR_EQ(u, "http://me:secret@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "me");
            }
            {
                // Boost.URL differentiates between an empty user an no user
                system::result<url> r = parse_uri_reference("http://me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user(""));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
                BOOST_TEST_CSTR_EQ(u, "http://@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
            }
            {
                system::result<url> r = parse_uri_reference("http://me:secret@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user(""));
                BOOST_TEST_CSTR_EQ(u, "http://:secret@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
            }
            // UTF-8 percent encoding with the userinfo encode set.
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // Boost.URL creates a string view which ends the buffer at the
                // first null char
                BOOST_TEST_NO_THROW(u.set_encoded_user("\x00" "\x01" "\t\n\r\x1f" " !\"#$%&\'()*+,-./09:;<=>?@AZ[\\]^_`az{|}~\x7f" "\x80" "\x81" "Éé"));
                // BOOST_TEST_CSTR_EQ(u, "http://%00%01%09%0A%0D%1F%20!%22%23$%&'()*+,-.%2F09%3A%3B%3C%3D%3E%3F%40AZ%5B%5C%5D%5E_%60az%7B%7C%7D~%7F%C2%80%C2%81%C3%89%C3%A9@example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "%00%01%09%0A%0D%1F%20!%22%23$%&'()*+,-.%2F09%3A%3B%3C%3D%3E%3F%40AZ%5B%5C%5D%5E_%60az%7B%7C%7D~%7F%C2%80%C2%81%C3%89%C3%A9");
                BOOST_TEST_CSTR_EQ(u, "http://@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
            }
            // Bytes already percent-encoded are left as-is.
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("%c3%89té"));
                // Boost.URL does not create unrequested absolute paths
                // BOOST_TEST_CSTR_EQ(u, "http://%c3%89t%C3%A9@example.net/");
                BOOST_TEST_CSTR_EQ(u, "http://%c3%89t%C3%A9@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "%c3%89t%C3%A9");
            }
            {
                system::result<url> r = parse_uri_reference("sc:///");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("x"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://x@/");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "x");
            }
            {
                system::result<url> r = parse_uri_reference("javascript://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("wario"));
                BOOST_TEST_CSTR_EQ(u, "javascript://wario@x/");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "wario");
            }
            {
                system::result<url> r = parse_uri_reference("file://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_user("test"));
                // BOOST_TEST_CSTR_EQ(u, "file://test/");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "");
                BOOST_TEST_CSTR_EQ(u, "file://test@test/");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "test");
            }
        }
        // password
        {
            // In Ada, no host means no password
            // Boost.URL creates an empty host
            {
                system::result<url> r = parse_uri_reference("file:///home/me/index.html");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("secret"));
                // BOOST_TEST_CSTR_EQ(u, "file:///home/me/index.html");
                // BOOST_TEST_CSTR_EQ(u.encoded_password(), "");
                BOOST_TEST_CSTR_EQ(u, "file://:secret@/home/me/index.html");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "secret");
            }
            // In Ada, no host means no password
            // Boost.URL creates an empty host
            {
                system::result<url> r = parse_uri_reference("unix:/run/foo.socket");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("secret"));
                // BOOST_TEST_CSTR_EQ(u, "unix:/run/foo.socket");
                // BOOST_TEST_CSTR_EQ(u.encoded_password(), "");
                BOOST_TEST_CSTR_EQ(u, "unix://:secret@/run/foo.socket");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "secret");
            }
            // In Ada, cannot-be-a-base means no password
            // Boost.URL creates an empty host
            {
                system::result<url> r = parse_uri_reference("mailto:me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("secret"));
                // BOOST_TEST_CSTR_EQ(u, "mailto:me@example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_password(), "");
                BOOST_TEST_CSTR_EQ(u, "mailto://:secret@/me@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "secret");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("secret"));
                // BOOST_TEST_CSTR_EQ(u, "http://:secret@example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_password(), "secret");
                BOOST_TEST_CSTR_EQ(u, "http://:secret@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "secret");
            }
            {
                system::result<url> r = parse_uri_reference("http://me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("secret"));
                // BOOST_TEST_CSTR_EQ(u, "http://me:secret@example.net/");
                BOOST_TEST_CSTR_EQ(u, "http://me:secret@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "secret");
            }
            {
                system::result<url> r = parse_uri_reference("http://:secret@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password(""));
                // Boost.URL distinguishes between empty password and no password
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                BOOST_TEST_CSTR_EQ(u, "http://:@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "");
            }
            {
                system::result<url> r = parse_uri_reference("http://me:secret@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password(""));
                // BOOST_TEST_CSTR_EQ(u, "http://me@example.net/");
                BOOST_TEST_CSTR_EQ(u, "http://me:@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "");
            }
            // UTF-8 percent encoding with the userinfo encode set.
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("\x00" "\x01" "\t\n\r\x1f" " !\"#$%&\'()*+,-./09:;<=>?@AZ[\\]^_`az{|}~\x7f" "\x80" "\x81" "Éé"));
                // BOOST_TEST_CSTR_EQ(u, "http://:%00%01%09%0A%0D%1F%20!%22%23$%&'()*+,-.%2F09%3A%3B%3C%3D%3E%3F%40AZ%5B%5C%5D%5E_%60az%7B%7C%7D~%7F%C2%80%C2%81%C3%89%C3%A9@example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_password(), "%00%01%09%0A%0D%1F%20!%22%23$%&'()*+,-.%2F09%3A%3B%3C%3D%3E%3F%40AZ%5B%5C%5D%5E_%60az%7B%7C%7D~%7F%C2%80%C2%81%C3%89%C3%A9");
                BOOST_TEST_CSTR_EQ(u, "http://:@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "");
            }
            // Bytes already percent-encoded are left as-is.
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("%c3%89té"));
                BOOST_TEST_CSTR_EQ(u, "http://:%c3%89t%C3%A9@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "%c3%89t%C3%A9");
            }
            {
                system::result<url> r = parse_uri_reference("sc:///");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("x"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_password(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://:x@/");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "x");
            }
            {
                system::result<url> r = parse_uri_reference("javascript://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("bowser"));
                BOOST_TEST_CSTR_EQ(u, "javascript://:bowser@x/");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "bowser");
            }
            {
                system::result<url> r = parse_uri_reference("file://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_password("test"));
                // BOOST_TEST_CSTR_EQ(u, "file://test/");
                // BOOST_TEST_CSTR_EQ(u.encoded_password(), "");
                BOOST_TEST_CSTR_EQ(u, "file://:test@test/");
                BOOST_TEST_CSTR_EQ(u.encoded_password(), "test");
            }
        }
        // host
        {
            // Non-special scheme
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("\x00" "");
                BOOST_TEST_NO_THROW(u.set_encoded_host("\x00" ""));
                // BOOST_TEST_CSTR_EQ(u, "sc://x/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "x");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "x");
                BOOST_TEST_CSTR_EQ(u, "sc:///");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("\t");
                BOOST_TEST_NO_THROW(u.set_encoded_host("\t"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%09/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%09");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%09");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("\n");
                BOOST_TEST_NO_THROW(u.set_encoded_host("\n"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%0A/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%0A");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%0A");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("\r");
                BOOST_TEST_NO_THROW(u.set_encoded_host("\r"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%0D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%0D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%0D");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port(" ");
                BOOST_TEST_NO_THROW(u.set_encoded_host(" "));
                // BOOST_TEST_CSTR_EQ(u, "sc://x/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "x");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "x");
                BOOST_TEST_CSTR_EQ(u, "sc://%20/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%20");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%20");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("#");
                BOOST_TEST_NO_THROW(u.set_encoded_host("#"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%23/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%23");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%23");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("/");
                BOOST_TEST_NO_THROW(u.set_encoded_host("/"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%2F/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%2F");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%2F");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("?");
                BOOST_TEST_NO_THROW(u.set_encoded_host("?"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%3F/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%3F");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%3F");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("@");
                BOOST_TEST_NO_THROW(u.set_encoded_host("@"));
                // BOOST_TEST_CSTR_EQ(u, "sc://x/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "x");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "x");
                BOOST_TEST_CSTR_EQ(u, "sc://%40/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%40");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%40");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("ß");
                BOOST_TEST_NO_THROW(u.set_encoded_host("ß"));
                BOOST_TEST_CSTR_EQ(u, "sc://%C3%9F/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%C3%9F");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%C3%9F");
            }
            // IDNA Nontransitional_Processing
            {
                system::result<url> r = parse_uri_reference("https://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("ß");
                BOOST_TEST_NO_THROW(u.set_encoded_host("ß"));
                // BOOST_TEST_CSTR_EQ(u, "https://xn--zca/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "xn--zca");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "xn--zca");
                BOOST_TEST_CSTR_EQ(u, "https://%C3%9F/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%C3%9F");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%C3%9F");
            }
            // Cannot-be-a-base means no host
            {
                system::result<url> r = parse_uri_reference("mailto:me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                // BOOST_TEST_CSTR_EQ(u, "mailto:me@example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u, "mailto://example.com/me@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
            }
            // Cannot-be-a-base means no host
            {
                system::result<url> r = parse_uri_reference("data:text/plain,Stuff");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.net");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.net"));
                // BOOST_TEST_CSTR_EQ(u, "data:text/plain,Stuff");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u, "data://example.net/text/plain,Stuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:8080");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_NO_THROW(u.set_port("8080"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/");
                BOOST_TEST_CSTR_EQ(u, "http://example.com:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "8080");
            }
            // Port number is unchanged if not specified in the new value
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/");
                BOOST_TEST_CSTR_EQ(u, "http://example.com:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "8080");
            }
            // Port number is unchanged if not specified
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_NO_THROW(u.set_port(""));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.com:");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // The empty host is not valid for special schemes
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("");
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
            }
            // The empty host is OK for non-special schemes
            {
                system::result<url> r = parse_uri_reference("view-source+http://example.net/foo");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("");
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                BOOST_TEST_CSTR_EQ(u, "view-source+http:///foo");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
            }
            // Path-only URLs can gain a host
            {
                system::result<url> r = parse_uri_reference("a:/foo");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.net");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.net"));
                BOOST_TEST_CSTR_EQ(u, "a://example.net/foo");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
            }
            // IPv4 address syntax is normalized in Ada
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("0x7F000001:8080");
                BOOST_TEST_NO_THROW(u.set_encoded_host("0x7F000001"));
                BOOST_TEST_NO_THROW(u.set_port("8080"));
                // BOOST_TEST_CSTR_EQ(u, "http://127.0.0.1:8080/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "127.0.0.1:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "127.0.0.1");
                BOOST_TEST_CSTR_EQ(u, "http://0x7F000001:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "0x7F000001:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "0x7F000001");
                BOOST_TEST_CSTR_EQ(u.port(), "8080");
            }
            // Valid IPv6 address syntax is normalized
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("[::0:01]:2");
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::0:01]"));
                BOOST_TEST_NO_THROW(u.set_port("2"));
                // Boost.URL also normalizes the IPv6 when parsing
                // The value is recoded from the binary address representation
                BOOST_TEST_CSTR_EQ(u, "http://[::1]:2");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "[::1]:2");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "[::1]");
                BOOST_TEST_CSTR_EQ(u.port(), "2");
            }
            // IPv6 literal address with port, crbug.com/1012416
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("[2001:db8::2]:4002");
                BOOST_TEST_NO_THROW(u.set_encoded_host("[2001:db8::2]"));
                BOOST_TEST_NO_THROW(u.set_port("4002"));
                // BOOST_TEST_CSTR_EQ(u, "http://[2001:db8::2]:4002/");
                BOOST_TEST_CSTR_EQ(u, "http://[2001:db8::2]:4002");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "[2001:db8::2]:4002");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "[2001:db8::2]");
                BOOST_TEST_CSTR_EQ(u.port(), "4002");
            }
            // Default port number is removed in Ada
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:80");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_NO_THROW(u.set_port("80"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/");
                BOOST_TEST_CSTR_EQ(u, "http://example.com:80");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u.port(), "80");
            }
            // Default port number is removed
            {
                system::result<url> r = parse_uri_reference("https://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:443");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_NO_THROW(u.set_port("443"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "https://example.com:443");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:443");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "443");
            }
            // Default port number is only removed for the relevant scheme
            {
                system::result<url> r = parse_uri_reference("https://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:80");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_NO_THROW(u.set_port("80"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.com:80/");
                BOOST_TEST_CSTR_EQ(u, "https://example.com:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "80");
            }
            // Port number is removed if new port is scheme default and existing URL has a non-default port
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:80");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_NO_THROW(u.set_port("80"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "http://example.com:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "80");
            }
            // Stuff after a / delimiter is ignored in Ada
            // Boost.URL encodes invalid chars
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com/stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com/stuff"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%2Fstuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%2Fstuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%2Fstuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a / delimiter is ignored in Ada
            // Boost.URL does not accept invalid port numbers
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:8080/stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_THROWS(u.set_port("8080/stuff"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a ? delimiter is ignored in Ada
            // Boost.URL encodes invalid host chars
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com?stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com?stuff"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%3Fstuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%3Fstuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%3Fstuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a ? delimiter is ignored
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:8080?stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_THROWS(u.set_port("8080?stuff"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a # delimiter is ignored in Ada
            // Boost.URL encodes invalid host chars
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com#stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com#stuff"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%23stuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%23stuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%23stuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a # delimiter is ignored in Ada
            // Boost.URL does not accept invalid port chars
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:8080#stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_THROWS(u.set_port("8080#stuff"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a \ delimiter is ignored for special schemes in Ada
            // Boost.URL encodes invalid host chars
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com\\stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com\\stuff"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%5Cstuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%5Cstuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%5Cstuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a \ delimiter is ignored for special schemes in Ada
            // Boost.URL does not accept invalid port chars
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:8080\\stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_THROWS(u.set_port("8080\\stuff"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // \ is not a delimiter for non-special schemes, but still forbidden in hosts
            // Boost.URL encodes invalid host chars
            {
                system::result<url> r = parse_uri_reference("view-source+http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com\\stuff");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com\\stuff"));
                // BOOST_TEST_CSTR_EQ(u, "view-source+http://example.net/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "view-source+http://example.com%5Cstuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%5Cstuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%5Cstuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // In Ada, anything other than ASCII digit stops the port parser in a setter but is not an error
            // Boost.URL does nto accept invalid port chars
            {
                system::result<url> r = parse_uri_reference("view-source+http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:8080stuff2");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_THROWS(u.set_port("8080stuff2"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "view-source+http://example.com:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "view-source+http://example.com/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // In Ada, anything other than ASCII digit stops the port parser in a setter but is not an error
            // Boost.URL does not accept invalid port chars
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:8080stuff2");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_THROWS(u.set_port("8080stuff2"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Anything other than ASCII digit stops the port parser in a setter but is not an error
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:8080+2");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_THROWS(u.set_port("8080+2"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Port numbers are 16 bit integers
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:65535");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_NO_THROW(u.set_port("65535"));
                BOOST_TEST_CSTR_EQ(u, "http://example.com:65535/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:65535");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "65535");
            }
            // Ada fails at port overflow
            // Boost.URL parses any valid port number as a string
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("example.com:65536");
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                BOOST_TEST_NO_THROW(u.set_port("65536"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "http://example.com:65536/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:65536");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "65536");
            }
            // Broken IPv6 is invalid in Ada
            // Boost.URL identifies a ref-name that needs encoding
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("[google.com]");
                BOOST_TEST_NO_THROW(u.set_encoded_host("[google.com]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5Bgoogle.com%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5Bgoogle.com%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5Bgoogle.com%5D");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("[::1.2.3.4x]");
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::1.2.3.4x]"));
                // Ada rejects invalid Ip Address
                // Boost.URL identifies a valid reg-name
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5B%3A%3A1.2.3.4x%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5B%3A%3A1.2.3.4x%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5B%3A%3A1.2.3.4x%5D");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("[::1.2.3.]");
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::1.2.3.]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5B%3A%3A1.2.3.%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5B%3A%3A1.2.3.%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5B%3A%3A1.2.3.%5D");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("[::1.2.]");
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::1.2.]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5B%3A%3A1.2.%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5B%3A%3A1.2.%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5B%3A%3A1.2.%5D");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("[::1.]");
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::1.]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5B%3A%3A1.%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5B%3A%3A1.%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5B%3A%3A1.%5D");
            }
            // Ada treats file as a special scheme
            {
                system::result<url> r = parse_uri_reference("file://y/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("x:123");
                BOOST_TEST_NO_THROW(u.set_encoded_host("x"));
                BOOST_TEST_NO_THROW(u.set_port("123"));
                // BOOST_TEST_CSTR_EQ(u, "file://y/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "y");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "y");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "file://x:123/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "x:123");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "x");
                BOOST_TEST_CSTR_EQ(u.port(), "123");
            }
            {
                system::result<url> r = parse_uri_reference("file://y/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("loc%41lhost");
                BOOST_TEST_NO_THROW(u.set_encoded_host("loc%41lhost"));
                // BOOST_TEST_CSTR_EQ(u, "file:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "file://loc%41lhost/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "loc%41lhost");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "loc%41lhost");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            {
                system::result<url> r = parse_uri_reference("file://hi/x");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("");
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                BOOST_TEST_CSTR_EQ(u, "file:///x");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            {
                system::result<url> r = parse_uri_reference("sc://test@test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("");
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                // BOOST_TEST_CSTR_EQ(u, "sc://test@test/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "test");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "test");
                BOOST_TEST_CSTR_EQ(u, "sc://test@/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "test");
            }
            {
                system::result<url> r = parse_uri_reference("sc://test:12/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("");
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                // BOOST_TEST_CSTR_EQ(u, "sc://test:12/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "test:12");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
                // BOOST_TEST_CSTR_EQ(u.port(), "12");
                BOOST_TEST_CSTR_EQ(u, "sc://:12/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), ":12");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u.port(), "12");
            }
            // Leading / is not stripped
            {
                system::result<url> r = parse_uri_reference("http://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("///bad.com");
                BOOST_TEST_NO_THROW(u.set_encoded_host("///bad.com"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "http://%2F%2F%2Fbad.com/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%2F%2F%2Fbad.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%2F%2F%2Fbad.com");
            }
            // Leading / is not stripped
            {
                system::result<url> r = parse_uri_reference("sc://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("///bad.com");
                BOOST_TEST_NO_THROW(u.set_encoded_host("///bad.com"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%2F%2F%2Fbad.com/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%2F%2F%2Fbad.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%2F%2F%2Fbad.com");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("a%C2%ADb");
                BOOST_TEST_NO_THROW(u.set_encoded_host("a%C2%ADb"));
                // BOOST_TEST_CSTR_EQ(u, "https://ab/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "ab");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "ab");
                BOOST_TEST_CSTR_EQ(u, "https://a%C2%ADb/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "a%C2%ADb");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "a%C2%ADb");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("\xad" "");
                BOOST_TEST_NO_THROW(u.set_encoded_host("\xad" ""));
                // BOOST_TEST_CSTR_EQ(u, "https://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "https://%AD/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%AD");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%AD");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("%C2%AD");
                BOOST_TEST_NO_THROW(u.set_encoded_host("%C2%AD"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "https://%C2%AD/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%C2%AD");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%C2%AD");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // u.set_encoded_host_and_port("xn--");
                BOOST_TEST_NO_THROW(u.set_encoded_host("xn--"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "https://xn--/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "xn--");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "xn--");
            }
        }
        // hostname
        {
            // Non-special scheme
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("\x00" ""));
                // BOOST_TEST_CSTR_EQ(u, "sc://x/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "x");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "x");
                BOOST_TEST_CSTR_EQ(u, "sc:///");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("\t"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%09/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%09");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%09");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("\n"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%0A/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%0A");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%0A");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("\r"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%0D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%0D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%0D");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host(" "));
                // BOOST_TEST_CSTR_EQ(u, "sc://x/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "x");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "x");
                BOOST_TEST_CSTR_EQ(u, "sc://%20/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%20");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%20");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("#"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%23/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%23");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%23");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("/"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%2F/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%2F");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%2F");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("?"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%3F/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%3F");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%3F");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("@"));
                // BOOST_TEST_CSTR_EQ(u, "sc://x/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "x");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "x");
                BOOST_TEST_CSTR_EQ(u, "sc://%40/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%40");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%40");
            }
            // Cannot-be-a-base means no host in Ada
            {
                system::result<url> r = parse_uri_reference("mailto:me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                // BOOST_TEST_CSTR_EQ(u, "mailto:me@example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u, "mailto://example.com/me@example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
            }
            // Cannot-be-a-base means no host in Ada
            {
                system::result<url> r = parse_uri_reference("data:text/plain,Stuff");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.net"));
                // BOOST_TEST_CSTR_EQ(u, "data:text/plain,Stuff");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u, "data://example.net/text/plain,Stuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com:8080/");
                BOOST_TEST_CSTR_EQ(u, "http://example.com:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u.port(), "8080");
            }
            // The empty host is not valid for special schemes
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
            }
            // The empty host is OK for non-special schemes
            {
                system::result<url> r = parse_uri_reference("view-source+http://example.net/foo");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                BOOST_TEST_CSTR_EQ(u, "view-source+http:///foo");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
            }
            // Path-only URLs can gain a host
            {
                system::result<url> r = parse_uri_reference("a:/foo");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.net"));
                BOOST_TEST_CSTR_EQ(u, "a://example.net/foo");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
            }
            // IPv4 address syntax is normalized
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("0x7F000001"));
                // BOOST_TEST_CSTR_EQ(u, "http://127.0.0.1:8080/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "127.0.0.1:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "127.0.0.1");
                BOOST_TEST_CSTR_EQ(u, "http://0x7F000001:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "0x7F000001:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "0x7F000001");
                BOOST_TEST_CSTR_EQ(u.port(), "8080");
            }
            // IPv6 address syntax is normalized
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::0:01]"));
                // Boost.URL also normalizes the IPv6 when parsing
                // The value is recoded from the binary address representation
                BOOST_TEST_CSTR_EQ(u, "http://[::1]");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "[::1]");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "[::1]");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // : delimiter invalidates entire value
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com:8080"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%3A8080/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%3A8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%3A8080");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // : delimiter invalidates entire value
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com:"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%3A:8080/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%3A:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%3A");
                BOOST_TEST_CSTR_EQ(u.port(), "8080");
            }
            // Stuff after a / delimiter is ignored
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com/stuff"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%2Fstuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%2Fstuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%2Fstuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a ? delimiter is ignored
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com?stuff"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%3Fstuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%3Fstuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%3Fstuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a # delimiter is ignored
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com#stuff"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%23stuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%23stuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%23stuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a \ delimiter is ignored for special schemes
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com\\stuff"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "http://example.com%5Cstuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%5Cstuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%5Cstuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // \ is not a delimiter for non-special schemes, but still forbidden in hosts
            {
                system::result<url> r = parse_uri_reference("view-source+http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("example.com\\stuff"));
                // BOOST_TEST_CSTR_EQ(u, "view-source+http://example.net/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "view-source+http://example.com%5Cstuff/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com%5Cstuff");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com%5Cstuff");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Broken IPv6
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("[google.com]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5Bgoogle.com%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5Bgoogle.com%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5Bgoogle.com%5D");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::1.2.3.4x]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5B%3A%3A1.2.3.4x%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5B%3A%3A1.2.3.4x%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5B%3A%3A1.2.3.4x%5D");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::1.2.3.]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5B%3A%3A1.2.3.%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5B%3A%3A1.2.3.%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5B%3A%3A1.2.3.%5D");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::1.2.]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5B%3A%3A1.2.%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5B%3A%3A1.2.%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5B%3A%3A1.2.%5D");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("[::1.]"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u, "http://%5B%3A%3A1.%5D/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%5B%3A%3A1.%5D");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%5B%3A%3A1.%5D");
            }
            {
                system::result<url> r = parse_uri_reference("file://y/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("x:123"));
                // BOOST_TEST_CSTR_EQ(u, "file://y/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "y");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "y");
                BOOST_TEST_CSTR_EQ(u, "file://x%3A123/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "x%3A123");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "x%3A123");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            {
                system::result<url> r = parse_uri_reference("file://y/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("loc%41lhost"));
                // BOOST_TEST_CSTR_EQ(u, "file:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "file://loc%41lhost/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "loc%41lhost");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "loc%41lhost");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            {
                system::result<url> r = parse_uri_reference("file://hi/x");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                BOOST_TEST_CSTR_EQ(u, "file:///x");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            {
                system::result<url> r = parse_uri_reference("sc://test@test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                // BOOST_TEST_CSTR_EQ(u, "sc://test@test/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "test");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
                // BOOST_TEST_CSTR_EQ(u.encoded_user(), "test");
                BOOST_TEST_CSTR_EQ(u, "sc://test@/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_user(), "test");
            }
            {
                system::result<url> r = parse_uri_reference("sc://test:12/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                // BOOST_TEST_CSTR_EQ(u, "sc://test:12/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "test:12");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "test");
                BOOST_TEST_CSTR_EQ(u, "sc://:12/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), ":12");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u.port(), "12");
            }
            // Drop /. from path
            {
                system::result<url> r = parse_uri_reference("non-spec:/.//p");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("h"));
                // BOOST_TEST_CSTR_EQ(u, "non-spec://h//p");
                BOOST_TEST_CSTR_EQ(u, "non-spec://h/.//p");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "h");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "h");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//p");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/.//p");
            }
            {
                system::result<url> r = parse_uri_reference("non-spec:/.//p");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host(""));
                // BOOST_TEST_CSTR_EQ(u, "non-spec:////p");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//p");
                BOOST_TEST_CSTR_EQ(u, "non-spec:///.//p");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/.//p");
            }
            // Leading / is not stripped
            {
                system::result<url> r = parse_uri_reference("http://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("///bad.com"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "http://%2F%2F%2Fbad.com/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%2F%2F%2Fbad.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%2F%2F%2Fbad.com");
            }
            // Leading / is not stripped
            {
                system::result<url> r = parse_uri_reference("sc://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("///bad.com"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://%2F%2F%2Fbad.com/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%2F%2F%2Fbad.com");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%2F%2F%2Fbad.com");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("a%C2%ADb"));
                // BOOST_TEST_CSTR_EQ(u, "https://ab/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "ab");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "ab");
                BOOST_TEST_CSTR_EQ(u, "https://a%C2%ADb/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "a%C2%ADb");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "a%C2%ADb");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("\xad"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "https://%AD/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%AD");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%AD");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("%C2%AD"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "https://%C2%AD/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "%C2%AD");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "%C2%AD");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.com/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_host("xn--"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.com/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.com");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.com");
                BOOST_TEST_CSTR_EQ(u, "https://xn--/");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "xn--");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "xn--");
            }
        }
        // port
        {
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("8080"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/");
                BOOST_TEST_CSTR_EQ(u, "http://example.net:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "8080");
            }
            // Port number is removed if empty is the new value
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port(""));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                BOOST_TEST_CSTR_EQ(u, "http://example.net:");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Default port number is removed
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("80"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "http://example.net:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "80");
            }
            // Default port number is removed
            {
                system::result<url> r = parse_uri_reference("https://example.net:4433");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("443"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "https://example.net:443");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:443");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "443");
            }
            // Default port number is only removed for the relevant scheme
            {
                system::result<url> r = parse_uri_reference("https://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("80"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net:80/");
                BOOST_TEST_CSTR_EQ(u, "https://example.net:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:80");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "80");
            }
            // Stuff after a / delimiter is ignored
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("8080/stuff"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.net/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a ? delimiter is ignored
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("8080?stuff"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.net/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a # delimiter is ignored
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("8080#stuff"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.net/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Stuff after a \ delimiter is ignored for special schemes
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("8080\\stuff"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.net/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Anything other than ASCII digit stops the port parser in a setter but is not an error
            {
                system::result<url> r = parse_uri_reference("view-source+http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("8080stuff2"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "view-source+http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "view-source+http://example.net/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Anything other than ASCII digit stops the port parser in a setter but is not an error
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("8080stuff2"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.net/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Anything other than ASCII digit stops the port parser in a setter but is not an error
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("8080+2"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.net/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "");
            }
            // Port numbers are 16 bit integers
            {
                system::result<url> r = parse_uri_reference("http://example.net/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("65535"));
                BOOST_TEST_CSTR_EQ(u, "http://example.net:65535/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:65535");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "65535");
            }
            // Port numbers are 16 bit integers, overflowing is an error in Ada
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("65536"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.net:65536/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:65536");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "65536");
            }
            // Setting port to a string that doesn't parse as a number
            {
                system::result<url> r = parse_uri_reference("http://example.net:8080/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("randomstring"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "http://example.net:8080/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "8080");
            }
            // Port numbers are 16 bit integers, overflowing is an error
            {
                system::result<url> r = parse_uri_reference("non-special://example.net:8080/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("65536"));
                // BOOST_TEST_CSTR_EQ(u, "non-special://example.net:8080/path");
                // BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:8080");
                // BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u, "non-special://example.net:65536/path");
                BOOST_TEST_CSTR_EQ(u.encoded_host_and_port(), "example.net:65536");
                BOOST_TEST_CSTR_EQ(u.encoded_host(), "example.net");
                BOOST_TEST_CSTR_EQ(u.port(), "65536");
            }
            {
                system::result<url> r = parse_uri_reference("file://test/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("12"));
                // BOOST_TEST_CSTR_EQ(u, "file://test/");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "file://test:12/");
                BOOST_TEST_CSTR_EQ(u.port(), "12");
            }
            {
                system::result<url> r = parse_uri_reference("file://localhost/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("12"));
                // BOOST_TEST_CSTR_EQ(u, "file:///");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "file://localhost:12/");
                BOOST_TEST_CSTR_EQ(u.port(), "12");
            }
            {
                system::result<url> r = parse_uri_reference("non-base:value");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("12"));
                // BOOST_TEST_CSTR_EQ(u, "non-base:value");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "non-base://:12/value");
                BOOST_TEST_CSTR_EQ(u.port(), "12");
            }
            {
                system::result<url> r = parse_uri_reference("sc:///");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("12"));
                // BOOST_TEST_CSTR_EQ(u, "sc:///");
                // BOOST_TEST_CSTR_EQ(u.port(), "");
                BOOST_TEST_CSTR_EQ(u, "sc://:12/");
                BOOST_TEST_CSTR_EQ(u.port(), "12");
            }
            {
                system::result<url> r = parse_uri_reference("sc://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("12"));
                BOOST_TEST_CSTR_EQ(u, "sc://x:12/");
                BOOST_TEST_CSTR_EQ(u.port(), "12");
            }
            {
                system::result<url> r = parse_uri_reference("javascript://x/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_port("12"));
                BOOST_TEST_CSTR_EQ(u, "javascript://x:12/");
                BOOST_TEST_CSTR_EQ(u.port(), "12");
            }
            // Leading u0009 on special scheme
            {
                system::result<url> r = parse_uri_reference("https://domain.com:443");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("\t8080"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u.port(), "443");
            }
            // Leading u0009 on non-special scheme
            {
                system::result<url> r = parse_uri_reference("wpt++://domain.com:443");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("\t8080"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u.port(), "8080");
                BOOST_TEST_CSTR_EQ(u.port(), "443");
            }
            // Should use all ascii prefixed characters as port
            {
                system::result<url> r = parse_uri_reference("https://www.google.com:4343");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_THROWS(u.set_port("4wpt"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u.port(), "4");
                BOOST_TEST_CSTR_EQ(u.port(), "4343");
            }
        }
        // pathname
        {
            // Opaque paths cannot be set
            {
                system::result<url> r = parse_uri_reference("mailto:me@example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("/foo"));
                // BOOST_TEST_CSTR_EQ(u, "mailto:me@example.net");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "me@example.net");
                BOOST_TEST_CSTR_EQ(u, "mailto:/foo");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/foo");
            }
            {
                system::result<url> r = parse_uri_reference("data:original");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("new value"));
                // BOOST_TEST_CSTR_EQ(u, "data:original");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "original");
                BOOST_TEST_CSTR_EQ(u, "data:new%20value");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "new%20value");
            }
            {
                system::result<url> r = parse_uri_reference("sc:original");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("new value"));
                // BOOST_TEST_CSTR_EQ(u, "sc:original");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "original");
                BOOST_TEST_CSTR_EQ(u, "sc:new%20value");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "new%20value");
            }
            // Special URLs cannot have their paths erased
            {
                system::result<url> r = parse_uri_reference("file:///some/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path(""));
                // BOOST_TEST_CSTR_EQ(u, "file:///");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
                BOOST_TEST_CSTR_EQ(u, "file://");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "");
            }
            // Non-special URLs can have their paths erased
            {
                system::result<url> r = parse_uri_reference("foo://somehost/some/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path(""));
                BOOST_TEST_CSTR_EQ(u, "foo://somehost");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "");
            }
            // Non-special URLs with an empty host can have their paths erased
            {
                system::result<url> r = parse_uri_reference("foo:///some/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path(""));
                BOOST_TEST_CSTR_EQ(u, "foo://");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "");
            }
            // Path-only URLs cannot have their paths erased
            {
                system::result<url> r = parse_uri_reference("foo:/some/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path(""));
                // BOOST_TEST_CSTR_EQ(u, "foo:/");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/");
                BOOST_TEST_CSTR_EQ(u, "foo:");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "");
            }
            // Path-only URLs always have an initial slash
            {
                system::result<url> r = parse_uri_reference("foo:/some/path");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("test"));
                // BOOST_TEST_CSTR_EQ(u, "foo:/test");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/test");
                BOOST_TEST_CSTR_EQ(u, "foo:test");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "test");
            }
            {
                system::result<url> r = parse_uri_reference("unix:/run/foo.socket?timeout=10");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("/var/log/../run/bar.socket"));
                // BOOST_TEST_CSTR_EQ(u, "unix:/var/run/bar.socket?timeout=10");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/var/run/bar.socket");
                BOOST_TEST_CSTR_EQ(u, "unix:/var/log/../run/bar.socket?timeout=10");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/var/log/../run/bar.socket");
                u.normalize_path();
                BOOST_TEST_CSTR_EQ(u, "unix:/var/run/bar.socket?timeout=10");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/var/run/bar.socket");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("home"));
                BOOST_TEST_CSTR_EQ(u, "https://example.net/home#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/home");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("../home"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/home#nav");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/home");
                BOOST_TEST_CSTR_EQ(u, "https://example.net/../home#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/../home");
                u.normalize_path();
                BOOST_TEST_CSTR_EQ(u, "https://example.net/../home#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/../home");
            }
            // \ is a segment delimiter for 'special' URLs
            {
                system::result<url> r = parse_uri_reference("http://example.net/home?lang=fr#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("\\a\\%2E\\b\\%2e.\\c"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/a/c?lang=fr#nav");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/a/c");
                BOOST_TEST_CSTR_EQ(u, "http://example.net/%5Ca%5C%2E%5Cb%5C%2e.%5Cc?lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%5Ca%5C%2E%5Cb%5C%2e.%5Cc");
            }
            // \ is *not* a segment delimiter for non-'special' URLs
            {
                system::result<url> r = parse_uri_reference("view-source+http://example.net/home?lang=fr#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("\\a\\%2E\\b\\%2e.\\c"));
                // BOOST_TEST_CSTR_EQ(u, "view-source+http://example.net/\\a\\%2E\\b\\%2e.\\c?lang=fr#nav");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/\\a\\%2E\\b\\%2e.\\c");
                BOOST_TEST_CSTR_EQ(u, "view-source+http://example.net/%5Ca%5C%2E%5Cb%5C%2e.%5Cc?lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%5Ca%5C%2E%5Cb%5C%2e.%5Cc");
            }
            // UTF-8 percent encoding with the default encode set. Tabs and newlines are removed.
            {
                system::result<url> r = parse_uri_reference("a:/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // Boost.URL does not accept invalid hexdig
                BOOST_TEST_THROWS(u.set_encoded_path("%00" "\x01" "\t\n\r\x1f" " !\"#$%&\'()*+,-./09:;<=>?@AZ[\\]^_`az{|}~\x7f" "\x80" "\x81" "Éé"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "a:/%00%01%1F%20!%22%23$%&'()*+,-./09:;%3C=%3E%3F@AZ[\\]^_%60az%7B|%7D~%7F%C2%80%C2%81%C3%89%C3%A9");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%00%01%1F%20!%22%23$%&'()*+,-./09:;%3C=%3E%3F@AZ[\\]^_%60az%7B|%7D~%7F%C2%80%C2%81%C3%89%C3%A9");
            }
            // Bytes already percent-encoded are left as-is, including %2E outside dotted segments.
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("%2e%2E%c3%89té"));
                BOOST_TEST_CSTR_EQ(u, "http://example.net/%2e%2E%c3%89t%C3%A9");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%2e%2E%c3%89t%C3%A9");
            }
            // ? needs to be encoded
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("?"));
                BOOST_TEST_CSTR_EQ(u, "http://example.net/%3F");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%3F");
            }
            // # needs to be encoded
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("#"));
                BOOST_TEST_CSTR_EQ(u, "http://example.net/%23");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%23");
            }
            // ? needs to be encoded, non-special scheme
            {
                system::result<url> r = parse_uri_reference("sc://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("?"));
                BOOST_TEST_CSTR_EQ(u, "sc://example.net/%3F");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%3F");
            }
            // # needs to be encoded, non-special scheme
            {
                system::result<url> r = parse_uri_reference("sc://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("#"));
                BOOST_TEST_CSTR_EQ(u, "sc://example.net/%23");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%23");
            }
            // ? doesn't mess up encoding
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("/?é"));
                BOOST_TEST_CSTR_EQ(u, "http://example.net/%3F%C3%A9");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%3F%C3%A9");
            }
            // # doesn't mess up encoding
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("/#é"));
                BOOST_TEST_CSTR_EQ(u, "http://example.net/%23%C3%A9");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%23%C3%A9");
            }
            // File URLs and (back)slashes
            {
                system::result<url> r = parse_uri_reference("file://monkey/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("\\\\"));
                // BOOST_TEST_CSTR_EQ(u, "file://monkey//");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//");
                BOOST_TEST_CSTR_EQ(u, "file://monkey/%5C%5C");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/%5C%5C");
            }
            // File URLs and (back)slashes
            {
                system::result<url> r = parse_uri_reference("file:///unicorn");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("//\\/"));
                // BOOST_TEST_CSTR_EQ(u, "file://////");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "////");
                BOOST_TEST_CSTR_EQ(u, "file:////%5C/");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "//%5C/");
            }
            // File URLs and (back)slashes
            {
                system::result<url> r = parse_uri_reference("file:///unicorn");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("//monkey/..//"));
                // BOOST_TEST_CSTR_EQ(u, "file://///");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "///");
                BOOST_TEST_CSTR_EQ(u, "file:////monkey/..//");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "//monkey/..//");
            }
            // Serialize /. in path
            {
                system::result<url> r = parse_uri_reference("non-spec:/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("/.//p"));
                BOOST_TEST_CSTR_EQ(u, "non-spec:/.//p");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//p");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/.//p");
            }
            {
                system::result<url> r = parse_uri_reference("non-spec:/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("/..//p"));
                // BOOST_TEST_CSTR_EQ(u, "non-spec:/.//p");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//p");
                BOOST_TEST_CSTR_EQ(u, "non-spec:/..//p");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/..//p");
            }
            {
                system::result<url> r = parse_uri_reference("non-spec:/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("//p"));
                BOOST_TEST_CSTR_EQ(u, "non-spec:/.//p");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "//p");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "/.//p");
            }
            // Drop /. from path
            {
                system::result<url> r = parse_uri_reference("non-spec:/.//");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("p"));
                // BOOST_TEST_CSTR_EQ(u, "non-spec:/p");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/p");
                BOOST_TEST_CSTR_EQ(u, "non-spec:p");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "p");
            }
            // Non-special URLs with non-opaque paths percent-encode U+0020
            {
                system::result<url> r = parse_uri_reference("data:/nospace");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("space "));
                // BOOST_TEST_CSTR_EQ(u, "data:/space%20");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/space%20");
                BOOST_TEST_CSTR_EQ(u, "data:space%20");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "space%20");
            }
            {
                system::result<url> r = parse_uri_reference("sc:/nospace");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_path("space "));
                // BOOST_TEST_CSTR_EQ(u, "sc:/space%20");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "/space%20");
                BOOST_TEST_CSTR_EQ(u, "sc:space%20");
                BOOST_TEST_CSTR_EQ(u.encoded_path(), "space%20");
            }
        }
        // search
        {
            {
                system::result<url> r = parse_uri_reference("https://example.net#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query("lang=fr"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/?lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "lang=fr");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query("lang=fr"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/?lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "lang=fr");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query("?lang=fr"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net?lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u, "https://example.net??lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "?lang=fr");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query("??lang=fr"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/??lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u, "https://example.net???lang=fr#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "??lang=fr");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query("?"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/?#nav");
                BOOST_TEST_CSTR_EQ(u, "https://example.net??#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "?");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query(""));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/#nav");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?#nav");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query(""));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query(""));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "");
            }
            // UTF-8 percent encoding with the query encode set. Tabs and newlines are removed.
            {
                system::result<url> r = parse_uri_reference("a:/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // Boost.URL does not accept invalid octets
                BOOST_TEST_THROWS(u.set_encoded_query("%00" "\x01" "\t\n\r\x1f" " !\"#$%&\'()*+,-./09:;<=>?@AZ[\\]^_`az{|}~\x7f" "\x80" "\x81" "Éé"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "a:/?%00%01%1F%20!%22%23$%&'()*+,-./09:;%3C=%3E?@AZ[\\]^_`az{|}~%7F%C2%80%C2%81%C3%89%C3%A9");
                // BOOST_TEST_CSTR_EQ(u.encoded_query(), "%00%01%1F%20!%22%23$%&'()*+,-./09:;%3C=%3E?@AZ[\\]^_`az{|}~%7F%C2%80%C2%81%C3%89%C3%A9");
            }
            // Bytes already percent-encoded are left as-is
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_query("%c3%89té"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/?%c3%89t%C3%A9");
                // BOOST_TEST_CSTR_EQ(u.encoded_query(), "%c3%89t%C3%A9");
                BOOST_TEST_CSTR_EQ(u, "http://example.net?%c3%89t%C3%A9");
                BOOST_TEST_CSTR_EQ(u.encoded_query(), "%c3%89t%C3%A9");
            }
            // Drop trailing spaces from trailing opaque paths
            {
                system::result<url> r = parse_uri_reference("data:space ?query");
                // Boost.URL does not accept parsing invalid path chars
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_encoded_query(""));
                // BOOST_TEST_CSTR_EQ(u, "data:space");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "space");
                // BOOST_TEST_CSTR_EQ(u.encoded_query(), "");
            }
            {
                system::result<url> r = parse_uri_reference("sc:space ?query");
                // Boost.URL does not accept parsing invalid path chars
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_encoded_query(""));
                // BOOST_TEST_CSTR_EQ(u, "sc:space");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "space");
                // BOOST_TEST_CSTR_EQ(u.encoded_query(), "");
            }
            // Do not drop trailing spaces from non-trailing opaque paths
            {
                system::result<url> r = parse_uri_reference("data:space  ?query#fragment");
                // Boost.URL does not accept parsing invalid path chars
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_encoded_query(""));
                // BOOST_TEST_CSTR_EQ(u, "data:space  #fragment");
                // BOOST_TEST_CSTR_EQ(u.encoded_query(), "");
            }
            {
                system::result<url> r = parse_uri_reference("sc:space  ?query#fragment");
                // Boost.URL does not accept parsing invalid path chars
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_encoded_query(""));
                // BOOST_TEST_CSTR_EQ(u, "sc:space  #fragment");
                // BOOST_TEST_CSTR_EQ(u.encoded_query(), "");
            }
        }
        // hash
        {
            {
                system::result<url> r = parse_uri_reference("https://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("main"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/#main");
                BOOST_TEST_CSTR_EQ(u, "https://example.net#main");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "main");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("main"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/#main");
                BOOST_TEST_CSTR_EQ(u, "https://example.net#main");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "main");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("#nav"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/?lang=en-US##nav");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?lang=en-US##nav");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "#nav");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("main"));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/?lang=en-US#main");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?lang=en-US#main");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "main");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment(""));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/?lang=en-US#");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?lang=en-US#");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "");
            }
            {
                system::result<url> r = parse_uri_reference("https://example.net?lang=en-US#nav");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment(""));
                // BOOST_TEST_CSTR_EQ(u, "https://example.net/?lang=en-US");
                BOOST_TEST_CSTR_EQ(u, "https://example.net?lang=en-US#");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("foo bar"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/#foo%20bar");
                BOOST_TEST_CSTR_EQ(u, "http://example.net#foo%20bar");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "foo%20bar");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("foo\"bar"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/#foo%22bar");
                BOOST_TEST_CSTR_EQ(u, "http://example.net#foo%22bar");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "foo%22bar");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("foo<bar"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/#foo%3Cbar");
                BOOST_TEST_CSTR_EQ(u, "http://example.net#foo%3Cbar");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "foo%3Cbar");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("foo>bar"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/#foo%3Ebar");
                BOOST_TEST_CSTR_EQ(u, "http://example.net#foo%3Ebar");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "foo%3Ebar");
            }
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("foo`bar"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/#foo%60bar");
                BOOST_TEST_CSTR_EQ(u, "http://example.net#foo%60bar");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "foo%60bar");
            }
            // Simple percent-encoding; tabs and newlines are removed
            {
                system::result<url> r = parse_uri_reference("a:/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                // Boost.URL does not accept invalid octets
                BOOST_TEST_THROWS(u.set_encoded_fragment("%00" "\x01" "\t\n\r\x1f" " !\"#$%&\'()*+,-./09:;<=>?@AZ[\\]^_`az{|}~\x7f" "\x80" "\x81" "Éé"), system::system_error);
                // BOOST_TEST_CSTR_EQ(u, "a:/#%00%01%1F%20!%22#$%&'()*+,-./09:;%3C=%3E?@AZ[\\]^_%60az{|}~%7F%C2%80%C2%81%C3%89%C3%A9");
                // BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "%00%01%1F%20!%22#$%&'()*+,-./09:;%3C=%3E?@AZ[\\]^_%60az{|}~%7F%C2%80%C2%81%C3%89%C3%A9");
            }
            // Percent-encode NULLs in fragment
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("a%00" "b"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/#a%00b");
                BOOST_TEST_CSTR_EQ(u, "http://example.net#a%00b");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "a%00b");
            }
            // Percent-encode NULLs in fragment
            {
                system::result<url> r = parse_uri_reference("non-spec:/");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("a%00" "b"));
                BOOST_TEST_CSTR_EQ(u, "non-spec:/#a%00b");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "a%00b");
            }
            // Bytes already percent-encoded are left as-is
            {
                system::result<url> r = parse_uri_reference("http://example.net");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("%c3%89té"));
                // BOOST_TEST_CSTR_EQ(u, "http://example.net/#%c3%89t%C3%A9");
                BOOST_TEST_CSTR_EQ(u, "http://example.net#%c3%89t%C3%A9");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "%c3%89t%C3%A9");
            }
            {
                system::result<url> r = parse_uri_reference("javascript:alert(1)");
                if (!BOOST_TEST(r)) return;
                url u = *r;
                BOOST_TEST_NO_THROW(u.set_encoded_fragment("castle"));
                BOOST_TEST_CSTR_EQ(u, "javascript:alert(1)#castle");
                BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "castle");
            }
            // Drop trailing spaces from trailing opaque paths
            {
                system::result<url> r = parse_uri_reference("data:space                   #fragment");
                // Boost.URL does not accept parsing invalid path chars
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_encoded_fragment(""));
                // BOOST_TEST_CSTR_EQ(u, "data:space");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "space");
                // BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "");
            }
            {
                system::result<url> r = parse_uri_reference("sc:space    #fragment");
                // Boost.URL does not accept parsing invalid path chars
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_encoded_fragment(""));
                // BOOST_TEST_CSTR_EQ(u, "sc:space");
                // BOOST_TEST_CSTR_EQ(u.encoded_path(), "space");
                // BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "");
            }
            // Do not drop trailing spaces from non-trailing opaque paths
            {
                system::result<url> r = parse_uri_reference("data:space  ?query#fragment");
                // Boost.URL does not accept parsing invalid path chars
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_encoded_fragment(""));
                // BOOST_TEST_CSTR_EQ(u, "data:space  ?query");
                // BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "");
            }
            {
                system::result<url> r = parse_uri_reference("sc:space  ?query#fragment");
                // Boost.URL does not accept parsing invalid path chars
                BOOST_TEST_NOT(r);
                // if (!BOOST_TEST(r)) return;
                // url u = *r;
                // BOOST_TEST_NO_THROW(u.set_encoded_fragment(""));
                // BOOST_TEST_CSTR_EQ(u, "sc:space  ?query");
                // BOOST_TEST_CSTR_EQ(u.encoded_fragment(), "");
            }
        }
    }

    void
    run()
    {
        // An adaptation of ada unit tests to ensure boost.url
        // can handle all the same URLs or at least the differences
        // are explicit
        // https://github.com/ada-url/ada
        adaBasicTests();
        urlTestData();
        adaSettersTests();
    }
};

TEST_SUITE(ada_test, "boost.url.compat.ada");

} // urls
} // boost
