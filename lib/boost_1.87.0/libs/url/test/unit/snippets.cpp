//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#include "test_suite.hpp"

#include <boost/container/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/core/detail/string_view.hpp>
#include <boost/url/string_view.hpp>
#include <boost/url/grammar/string_token.hpp>

//[snippet_headers_1
#include <boost/url.hpp>
//]

#include <algorithm>
#include <iostream>
#include <vector>
#include <cctype>

//[snippet_headers_3
#include <boost/url.hpp>
using namespace boost::urls;
//]

#include <iostream>

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

void
using_url_views()
{
    {
        //[snippet_accessing_1
        url_view u( "https://user:pass@example.com:443/path/to/my%2dfile.txt?id=42&name=John%20Doe+Jingleheimer%2DSchmidt#page%20anchor" );
        assert(u.scheme() == "https");
        assert(u.authority().buffer() == "user:pass@example.com:443");
        assert(u.userinfo() == "user:pass");
        assert(u.user() == "user");
        assert(u.password() == "pass");
        assert(u.host() == "example.com");
        assert(u.port() == "443");
        assert(u.path() == "/path/to/my-file.txt");
        assert(u.query() == "id=42&name=John Doe Jingleheimer-Schmidt");
        assert(u.fragment() == "page anchor");
        //]
    }

    //[code_urls_parsing_1
    boost::core::string_view s = "https://user:pass@example.com:443/path/to/my%2dfile.txt?id=42&name=John%20Doe+Jingleheimer%2DSchmidt#page%20anchor";
    //]

    {
        //[code_urls_parsing_2
        boost::system::result<url_view> r = parse_uri( s );
        //]
        //[snippet_parsing_3
        url_view u = r.value();
        //]
        boost::ignore_unused(u);
    }

    {
        boost::system::result<url_view> r = parse_uri( s );
        //[snippet_parsing_4
        url_view u = *r;
        //]
        boost::ignore_unused(u);
    }

    url_view u( s );
    //[snippet_accessing_1b
    for (auto seg: u.segments())
        std::cout << seg << "\n";
    std::cout << "\n";

    for (auto param: u.params())
        std::cout << param.key << ": " << param.value << "\n";
    std::cout << "\n";
    //]

    {
        //[snippet_token_1
        std::string h = u.host();
        assert(h == "example.com");
        //]
        boost::ignore_unused(h);
    }

    {
        //[snippet_token_2
        std::string h = "host: ";
        u.host(string_token::append_to(h));
        assert(h == "host: example.com");
        //]
        boost::ignore_unused(h);
    }

    {
        //[snippet_accessing_2a
        url_view u1 = parse_uri( "http://www.example.com" ).value();
        assert(u1.fragment().empty());
        assert(!u1.has_fragment());
        //]
        boost::ignore_unused(u1);
    }

    {
        //[snippet_accessing_2b
        url_view u2 = parse_uri( "http://www.example.com/#" ).value();
        assert(u2.fragment().empty());
        assert(u2.has_fragment());
        //]
    }

    {
        //[snippet_accessing_3a
        url_view u1 = parse_uri( "http://www.example.com" ).value();
        std::cout << "has fragment 1 : " << u1.has_fragment() << "\n";
        std::cout << "fragment 1 : " << u1.fragment() << "\n\n";
        //]

        //[snippet_accessing_3b
        url_view u2 = parse_uri( "http://www.example.com/#" ).value();
        std::cout << "has fragment 2 : " << u2.has_fragment() << "\n";
        std::cout << "fragment 2 : " << u2.fragment() << "\n\n";
        //]
    }

    {
        //[snippet_accessing_4
        std::cout <<
            "url       : " << u                     << "\n"
            "scheme    : " << u.scheme()            << "\n"
            "authority : " << u.encoded_authority() << "\n"
            "userinfo  : " << u.encoded_userinfo()  << "\n"
            "user      : " << u.encoded_user()      << "\n"
            "password  : " << u.encoded_password()  << "\n"
            "host      : " << u.encoded_host()      << "\n"
            "port      : " << u.port()              << "\n"
            "path      : " << u.encoded_path()      << "\n"
            "query     : " << u.encoded_query()     << "\n"
            "fragment  : " << u.encoded_fragment()  << "\n";
        //]
    }

    {
        //[snippet_decoding_1
        decode_view dv("id=42&name=John%20Doe%20Jingleheimer%2DSchmidt");
        std::cout << dv << "\n";
        //]
    }
    {
        url u1 = u;
        url u2 = u;
        //[snippet_decoding_2
        u1.set_host(u2.host());
        //]
        std::cout << u1 << "\n";
    }
    {
        //[snippet_decoding_3
        boost::filesystem::path p;
        for (auto seg: u.segments())
            p.append(seg.begin(), seg.end());
        std::cout << "path: " << p << "\n";
        //]
    }
// transparent std::equal_to<> required
#if BOOST_CXX_VERSION >= 201402L && !defined(BOOST_CLANG)
    {
        auto handle_route = [](
            std::vector<std::string> const&,
            url_view)
        {};

        //[snippet_decoding_4a
        auto match = [](
            std::vector<std::string> const& route,
            url_view u)
        {
            auto segs = u.segments();
            if (route.size() != segs.size())
                return false;
            return std::equal(
                route.begin(),
                route.end(),
                segs.begin());
        };
        //]
        //[snippet_decoding_4b
        std::vector<std::string> route =
            {"community", "reviews.html"};
        if (match(route, u))
        {
            handle_route(route, u);
        }
        //]
    }
#endif
    {
        //[snippet_compound_elements_1
        segments_encoded_view segs = u.encoded_segments();
        for( auto v : segs )
        {
            std::cout << v << "\n";
        }
        //]
    }

    {
        //[snippet_encoded_compound_elements_1
        segments_encoded_view segs = u.encoded_segments();

        for( pct_string_view v : segs )
        {
            decode_view dv = *v;
            std::cout << dv << "\n";
        }
        //]
    }

    {
        //[snippet_encoded_compound_elements_2
        params_encoded_view params_ref = u.encoded_params();

        for( auto v : params_ref )
        {
            decode_view dk(v.key);
            decode_view dv(v.value);

            std::cout <<
                "key = " << dk <<
                ", value = " << dv << "\n";
        }
        //]
    }
}

void
using_urls()
{
    boost::core::string_view s = "https://user:pass@www.example.com:443/path/to/my%2dfile.txt?id=42&name=John%20Doe#page%20anchor";

    //[snippet_quicklook_modifying_1
    url u = parse_uri( s ).value();
    //]

    //[snippet_quicklook_modifying_1b
    static_url<1024> su = parse_uri( s ).value();
    //]
    (void)su;

    //[snippet_quicklook_modifying_2
    u.set_scheme( "https" );
    //]

    //[snippet_quicklook_modifying_3
    u.set_scheme_id( scheme::https ); // equivalent to u.set_scheme( "https" );
    //]

    //[snippet_quicklook_modifying_4
    u.set_host_ipv4( ipv4_address( "192.168.0.1" ) )
        .set_port_number( 8080 )
        .remove_userinfo();
    std::cout << u << "\n";
    //]

    //[snippet_quicklook_modifying_5
    params_ref p = u.params();
    p.replace(p.find("name"), {"name", "John Doe"});
    std::cout << u << "\n";
    //]
}

void
parsing_urls()
{
    {
        auto handle_error = [](boost::system::error_code e)
        {
            boost::ignore_unused(e);
        };
        //[snippet_parsing_url_1
        boost::system::result< url > ru = parse_uri_reference( "https://www.example.com/path/to/file.txt" );
        if ( ru )
        {
            url u = *ru;
            assert(u.encoded_path() == "/path/to/file.txt");
        }
        else
        {
            boost::system::error_code e = ru.error();
            handle_error(e);
        }
        //]
    }

    {
        auto handle_error = [](boost::system::system_error& e)
        {
            boost::ignore_unused(e);
        };
        //[snippet_parsing_url_1b
        try
        {
            url u = parse_uri_reference( "https://www.example.com/path/to/file.txt" ).value();
            assert(u.encoded_path() == "/path/to/file.txt");
        }
        catch (boost::system::system_error &e)
        {
            handle_error(e);
        }
        //]
    }

    {
        //[snippet_parsing_url_1b2
        boost::system::result< url > ru = parse_uri_reference( "https://www.example.com/path/to/file.txt" );
        if ( ru.has_value() )
        {
            url u = *ru;
            assert(u.encoded_path() == "/path/to/file.txt");
        }
        //]
        boost::ignore_unused(ru);
    }

    {
        //[snippet_parsing_url_1bb
        url u = parse_uri_reference( "https://www.example.com/path/to/file.txt" ).value();

        assert(u.encoded_path() == "/path/to/file.txt");
        //]
        boost::ignore_unused(u);
    }

    {
        //[snippet_parsing_url_1bc
        boost::system::result< url > rv = parse_uri_reference( "https://www.example.com/path/to/file.txt" );

        static_assert( std::is_convertible< boost::system::result< url_view >, boost::system::result< url > >::value, "" );
        //]
        boost::ignore_unused(rv);
    }

    {
        //[snippet_parsing_url_1bd
        boost::system::result< static_url<1024> > rv = parse_uri_reference( "https://www.example.com/path/to/file.txt" );

        static_assert( std::is_convertible< boost::system::result< static_url<1024> >, boost::system::result< url > >::value, "" );
        //]
        boost::ignore_unused(rv);
    }

    {
        //[snippet_parsing_url_1c
        boost::system::result< url_view > r0 = parse_relative_ref( "/path/to/file.txt" );
        assert( r0.has_value() );
        //]
        //[snippet_parsing_url_1d
        boost::system::result< url_view > r1 = parse_uri( "https://www.example.com" );
        assert( r1.has_value() );
        url dest;
        resolve(r1.value(), r0.value(), dest);
        assert(dest.buffer() == "https://www.example.com/path/to/file.txt");
        //]
        boost::ignore_unused(dest);
    }

    //[snippet_parsing_url_2
    // This will hold our copy
    std::shared_ptr<url_view const> sp;
    {
        std::string s = "/path/to/file.txt";

        // result::value() will throw an exception if an error occurs
        url_view u = parse_relative_ref( s ).value();

        // create a copy with ownership and string lifetime extension
        sp = u.persist();

        // At this point the string goes out of scope
    }

    // but `*sp` remains valid since it has its own copy
    std::cout << *sp << "\n";
    //]

    {
        //[snippet_parsing_url_3
        // This will hold our mutable copy
        url v;
        {
            std::string s = "/path/to/file.txt";

            // result::value() will throw an exception if an error occurs
            v = parse_relative_ref(s).value();

            // At this point the string goes out of scope
        }

        // but `v` remains valid since it has its own copy
        std::cout << v << "\n";

        // and it's mutable
        v.set_fragment("anchor");

        // path/to/file.txt#anchor
        std::cout << v << "\n";
        //]
    }
}

void
parsing_components()
{
    {
        //[snippet_components_2a
        url_view u( "https://www.ietf.org/rfc/rfc2396.txt" );
        assert(u.scheme() == "https");
        assert(u.host() == "www.ietf.org");
        assert(u.path() == "/rfc/rfc2396.txt");
        //]
    }
    {
        //[snippet_components_2b
        url_view u( "ftp://ftp.is.co.za/rfc/rfc1808.txt" );
        assert(u.scheme() == "ftp");
        assert(u.host() == "ftp.is.co.za");
        assert(u.path() == "/rfc/rfc1808.txt");
        //]
    }
    {
        //[snippet_components_2c
        url_view u( "mailto:John.Doe@example.com" );
        assert(u.scheme() == "mailto");
        assert(u.path() == "John.Doe@example.com");
        //]
    }
    {
        //[snippet_components_2d
        url_view u( "urn:isbn:096139210x" );
        assert(u.scheme() == "urn");
        assert(u.path() == "isbn:096139210x");
        //]
    }
    {
        //[snippet_components_2e
        url_view u( "magnet:?xt=urn:btih:d2474e86c95b19b8bcfdb92bc12c9d44667cfa36" );
        assert(u.scheme() == "magnet");
        assert(u.path() == "");
        assert(u.query() == "xt=urn:btih:d2474e86c95b19b8bcfdb92bc12c9d44667cfa36");
        //]
    }
}

void
formatting_components()
{
    // I have spent a lot of time on this and have no
    // idea how to fix this bug in GCC 4.8 and GCC 5.0
    // without help from the pros.
#if !BOOST_WORKAROUND( BOOST_GCC_VERSION, < 60000 )
    {
        //[snippet_format_1
        url u = format("{}://{}:{}/rfc/{}", "https", "www.ietf.org", 80, "rfc2396.txt");
        assert(u.buffer() == "https://www.ietf.org:80/rfc/rfc2396.txt");
        //]
    }

    {
        //[snippet_format_2
        url u = format("https://{}/{}", "www.boost.org", "Hello world!");
        assert(u.buffer() == "https://www.boost.org/Hello%20world!");
        //]
    }

    {
        //[snippet_format_3a
        url u = format("{}:{}", "mailto", "someone@example.com");
        assert(u.buffer() == "mailto:someone@example.com");
        assert(u.scheme() == "mailto");
        assert(u.path() == "someone@example.com");
        //]
    }
    {
        //[snippet_format_3b
        url u = format("{}{}", "mailto:", "someone@example.com");
        assert(u.buffer() == "mailto%3Asomeone@example.com");
        assert(!u.has_scheme());
        assert(u.path() == "mailto:someone@example.com");
        assert(u.encoded_path() == "mailto%3Asomeone@example.com");
        //]
    }

    {
        //[snippet_format_4
        static_url<50> u;
        format_to(u, "{}://{}:{}/rfc/{}", "https", "www.ietf.org", 80, "rfc2396.txt");
        assert(u.buffer() == "https://www.ietf.org:80/rfc/rfc2396.txt");
        //]
    }

    {
        //[snippet_format_5a
        url u = format("{0}://{2}:{1}/{3}{4}{3}", "https", 80, "www.ietf.org", "abra", "cad");
        assert(u.buffer() == "https://www.ietf.org:80/abracadabra");
        //]
    }
    {
        //[snippet_format_5b
        url u = format("https://example.com/~{username}", arg("username", "mark"));
        assert(u.buffer() == "https://example.com/~mark");
        //]
    }

    {
        //[snippet_format_5c
        boost::core::string_view fmt = "{scheme}://{host}:{port}/{dir}/{file}";
        url u = format(fmt, {{"scheme", "https"}, {"port", 80}, {"host", "example.com"}, {"dir", "path/to"}, {"file", "file.txt"}});
        assert(u.buffer() == "https://example.com:80/path/to/file.txt");
        //]
    }
#endif
}

void
parsing_scheme()
{
    {
        //[snippet_parsing_scheme_1
        url_view u("mailto:name@email.com" );
        assert( u.has_scheme() );
        assert( u.scheme() == "mailto" );
        //]
        boost::ignore_unused(u);
    }
    {
        //[snippet_parsing_scheme_3
        url_view u("file://host/path/to/file" );
        assert( u.scheme_id() == scheme::file );
        //]
        boost::ignore_unused(u);
    }
}

void
parsing_authority()
{
    {
        //[snippet_parsing_authority_1
        url_view u( "https:///path/to_resource" );
        assert( u.authority().buffer() == "");
        assert( u.has_authority() );
        assert( u.path() == "/path/to_resource" );
        //]
    }
    {
        //[snippet_parsing_authority_2
        url_view u("https://www.boost.org" );
        std::cout << "scheme:        " << u.scheme()            << "\n"
                     "has authority: " << u.has_authority()     << "\n"
                     "authority:     " << u.authority()         << "\n"
                     "path:          " << u.path()              << "\n";
        //]
    }
    {
        //[snippet_parsing_authority_3a
        url_view u( "https://www.boost.org/users/download/" );
        assert(u.has_authority());
        authority_view a = u.authority();
        //]
        //[snippet_parsing_authority_3b
        assert(a.host() == "www.boost.org");
        //]
    }
    {
        //[snippet_parsing_authority_4
        url_view u( "https://www.boost.org/" );
        std::cout << "scheme:        " << u.scheme()            << "\n"
                     "has authority: " << u.has_authority()     << "\n"
                     "authority:     " << u.authority()         << "\n"
                     "path:          " << u.path()              << "\n";
        //]
    }
    {
        //[snippet_parsing_authority_5
        url_view u( "mailto:John.Doe@example.com" );
        std::cout << "scheme:        " << u.scheme()            << "\n"
                     "has authority: " << u.has_authority()     << "\n"
                     "authority:     " << u.authority()         << "\n"
                     "path:          " << u.path()              << "\n";
        //]
    }
    {
        //[snippet_parsing_authority_6
        url_view u( "mailto://John.Doe@example.com" );
        std::cout << u << "\n"
            "scheme:        " << u.scheme()   << "\n"
            "has authority: " << u.has_authority()     << "\n"
            "authority:     " << u.authority()         << "\n"
            "path:          " << u.path()              << "\n";
        //]
    }
    {
        //[snippet_parsing_authority_7
        url_view u( "https://john.doe@www.example.com:123/forum/questions/" );
        std::cout << "scheme:        " << u.scheme()   << "\n"
            "has authority: " << u.has_authority()     << "\n"
            "authority:     " << u.authority()         << "\n"
            "host:          " << u.host()              << "\n"
            "userinfo:      " << u.userinfo()          << "\n"
            "port:          " << u.port()              << "\n"
            "path:          " << u.path()              << "\n";
        //]
    }
    {
        //[snippet_parsing_authority_8
        url_view u( "https://john.doe@www.example.com:123/forum/questions/" );
        assert(u.host() == "www.example.com");
        assert(u.port() == "123");
        //]
    }
    {
        //[snippet_parsing_authority_9
        url_view u( "https://john.doe@192.168.2.1:123/forum/questions/" );
        assert(u.host() == "192.168.2.1");
        assert(u.port() == "123");
        //]
    }
    {
        //[snippet_parsing_authority_9b
        url_view u( "https://www.example.com" );
        assert(u.host() == "www.example.com");
        assert(u.host() == u.encoded_host());
        //]
    }
    {
        struct resolve_f {
            boost::core::string_view
            operator()(boost::core::string_view h)
            {
                return h;
            }
        } resolve;
        struct write_request_f {
            void operator()(boost::core::string_view) {}
            void operator()(ipv4_address) {}
            void operator()(ipv6_address) {}
        } write_request;

        //[snippet_parsing_authority_10
        url_view u( "https://www.boost.org/users/download/" );
        switch (u.host_type())
        {
        case host_type::name:
            write_request(resolve(u.host()));
            break;
        case host_type::ipv4:
            write_request(u.host_ipv4_address());
            break;
        case host_type::ipv6:
            write_request(u.host_ipv6_address());
            break;
        default:
            break;
        }
        //]
    }
    {
        //[snippet_parsing_authority_10a
        url_view u( "https:///path/to_resource" );
        assert( u.has_authority() );
        assert( u.authority().buffer().empty() );
        assert( u.path() == "/path/to_resource" );
        //]
    }
    {
        //[snippet_parsing_authority_10b
        url_view u( "https://www.boost.org" );
        assert( u.host() == "www.boost.org" );
        assert( u.path().empty() );
        //]
    }
    {
        //[snippet_parsing_authority_10c
        url_view u( "https://www.boost.org/users/download/" );
        assert( u.host() == "www.boost.org" );
        assert( u.path() == "/users/download/" );
        //]
    }
    {
        //[snippet_parsing_authority_10d
        url_view u( "https://www.boost.org/" );
        assert( u.host() == "www.boost.org" );
        assert( u.path() == "/" );
        //]
    }
    {
        //[snippet_parsing_authority_10e
        url_view u( "mailto:John.Doe@example.com" );
        assert( !u.has_authority() );
        assert( u.path() == "John.Doe@example.com" );
        //]
    }
    {
        //[snippet_parsing_authority_10f
        url_view u( "mailto://John.Doe@example.com" );
        assert( u.authority().buffer() == "John.Doe@example.com" );
        assert( u.path().empty() );
        //]
    }
    {
        //[snippet_parsing_authority_11a
        url_view u( "https://john.doe@www.example.com:123/forum/questions/" );
        assert(u.userinfo() == "john.doe");
        assert(u.port() == "123");
        //]
    }
    {
        //[snippet_parsing_authority_11b
        url_view u( "https://john.doe:123456@www.somehost.com/forum/questions/" );
        assert(u.userinfo() == "john:doe");
        assert(u.user() == "john");
        assert(u.password() == "doe");
        //]
    }
    {
        //[snippet_parsing_authority_12
        authority_view a = parse_authority( "www.example.com:80" ).value();
        assert(!a.has_userinfo());
        assert(a.host() == "www.example.com");
        assert(a.port() == "80");
        //]
    }
    {
        //[snippet_parsing_authority_13
        authority_view a = parse_authority( "user:pass@www.example.com:443" ).value();
        assert(a.userinfo() == "user:pass");
        assert(a.user() == "user");
        assert(a.password() == "pass");
        assert(a.host() == "www.example.com");
        assert(a.port() == "443");
        //]
    }
}

void
parsing_path()
{
    {
        //[snippet_parsing_path_0
        url_view u("http://www.example.com/path/to/file.txt");
        assert(u.path() == "/path/to/file.txt");
        segments_view segs = u.segments();
        auto it = segs.begin();
        assert( *it++ == "path" );
        assert( *it++ == "to" );
        assert( *it   == "file.txt" );
        //]
        boost::ignore_unused(it);
    }
    {
        //[snippet_parsing_path_1
        url_view u("https://www.boost.org/doc/the%20libs/");
        assert(u.path() == "/doc/the libs/");
        assert(u.encoded_path() == "/doc/the%20libs/");
        //]

        //[snippet_parsing_path_1_b
        std::cout << u.encoded_segments().size() << " segments\n";
        for (auto seg: u.encoded_segments())
            std::cout << "segment: " << seg << "\n";
        //]
    }
    {
        //[snippet_parsing_path_2
        url_view u("https://www.boost.org/doc/libs");
        std::cout << u.segments().size() << " segments\n";
        for (auto seg: u.segments())
            std::cout << "segment: " << seg << "\n";
        //]
    }
    {
        //[snippet_parsing_path_3
        url_view u("https://www.boost.org");
        assert(u.path().empty());
        assert(u.encoded_path().empty());
        //]
    }
    {
        //[snippet_parsing_path_3a
        assert( url_view("urn:isbn:096139210x").path() == "isbn:096139210x" );
        //]
    }
    {
        //[snippet_parsing_path_4
        /*
        url_view u("https://www.boost.org//doc///libs");
        std::cout << u << "\n"
                  "path:             " << u.encoded_path()     << "\n"
                  "encoded segments: " << u.encoded_segments() << "\n"
                  "segments:         " << u.segments()         << "\n";
        std::cout << u.encoded_segments().size() << " segments\n";
        for (auto seg: u.encoded_segments())
            std::cout << "segment: " << seg << "\n";
        */
        //]
    }

    {
        {
            //[snippet_parsing_path_5_a
            boost::core::string_view s = "https://www.boost.org";
            url_view u(s);
            std::cout << u << "\n"
                      << "path:     " << u.encoded_host()            << "\n"
                      << "path:     " << u.encoded_path()            << "\n"
                      << "segments: " << u.encoded_segments().size() << "\n";
            //]
        }
        {
            //[snippet_parsing_path_5_b
            boost::core::string_view s = "https://www.boost.org/";
            url_view u(s);
            std::cout << u << "\n"
                      << "host:     " << u.encoded_host()            << "\n"
                      << "path:     " << u.encoded_path()            << "\n"
                      << "segments: " << u.encoded_segments().size() << "\n";
            //]
        }
        {
            //[snippet_parsing_path_5_c
            boost::core::string_view s = "https://www.boost.org//";
            url_view u(s);
            std::cout << u << "\n"
                      << "host:     " << u.encoded_host()            << "\n"
                      << "path:     " << u.encoded_path()            << "\n"
                      << "segments: " << u.encoded_segments().size() << "\n";
            //]
        }
    }

    {
        //[snippet_parsing_path_6
        url_view u("https://www.boost.org//doc/libs/");
        std::cout << u << "\n"
                  "authority: " << u.encoded_authority() << "\n"
                  "path:      " << u.encoded_path()      << "\n";
        std::cout << u.encoded_segments().size() << " segments\n";
        for (auto seg: u.encoded_segments())
            std::cout << "segment: " << seg << "\n";
        //]
    }

    {
        //[snippet_parsing_path_7
        url_view u("https://doc/libs/");
        std::cout << u << "\n"
                  "authority: " << u.encoded_authority() << "\n"
                  "path:      " << u.encoded_path()      << "\n";
        std::cout << u.encoded_segments().size() << " segments\n";
        for (auto seg: u.encoded_segments())
            std::cout << "segment: " << seg << "\n";
        //]
    }

    {
        //[snippet_parsing_path_8
        url_view u("https://www.boost.org/doc@folder/libs:boost");
        std::cout << u << "\n"
                  "authority: " << u.encoded_authority() << "\n"
                  "path:      " << u.encoded_path()      << "\n";
        std::cout << u.encoded_segments().size() << " segments\n";
        for (auto seg: u.encoded_segments())
            std::cout << "segment: " << seg << "\n";
        //]
    }

    {
        //[snippet_parsing_path_9
        /*
        segments_view segs = parse_path("/doc/libs").value();
        assert( segs.size() == 2 );
        */
        //]
        //boost::ignore_unused(segs);
    }

    {
        //[snippet_parsing_path_use_case_1
        url_view u("https://www.boost.org/doc/libs/");
        assert(u.host() == "www.boost.org");
        assert(u.path() == "/doc/libs/");
        //]
        boost::ignore_unused(u);
    }
    {
        //[snippet_parsing_path_use_case_2
        assert( parse_uri("https://www.boost.org").has_value() );
        assert( parse_uri("https://www.boost.org/").has_value() );
        assert( parse_uri("https://www.boost.org//").has_value() );
        //]
    }
    {
        //[snippet_parsing_path_use_case_3
        assert( url_view("https://www.boost.org").path().empty() );
        //]
    }
    {
        //[snippet_parsing_path_use_case_4
        url_view u("https://www.boost.org/doc@folder/libs:boost");
        assert(u.path() == "/doc@folder/libs:boost");
        //]
        boost::ignore_unused(u);
    }
    {
        //[snippet_parsing_path_use_case_5
        url_view u("https://www.boost.org/doc/libs/");
        segments_view segs = u.segments();
        auto it = segs.begin();
        assert(*it++ == "doc");
        assert(*it++ == "libs");
        assert(*it == "");
        //]
        boost::ignore_unused(it);
    }
    {
        //[snippet_parsing_path_use_case_6
        url_view u("https://www.boost.org//doc///libs");
        segments_view segs = u.segments();
        auto it = segs.begin();
        assert(*it++ == "");
        assert(*it++ == "doc");
        assert(*it++ == "");
        assert(*it++ == "");
        assert(*it == "libs");
        //]
        boost::ignore_unused(it);
    }
    {
        //[snippet_parsing_path_use_case_7
        url_view u("https://www.boost.org//doc/libs/");
        segments_view segs = u.segments();
        auto it = segs.begin();
        assert(*it++ == "");
        assert(*it++ == "doc");
        assert(*it++ == "libs");
        assert(*it == "");
        //]
        boost::ignore_unused(it);
    }
    {
        //[snippet_parsing_path_use_case_8
        url_view u("https://doc/libs/");
        assert(u.host() == "doc");
        segments_view segs = u.segments();
        auto it = segs.begin();
        assert(*it++ == "libs");
        assert(*it == "");
        //]
        boost::ignore_unused(it);
    }
}

void
parsing_query()
{
    {
        //[snippet_parsing_query_0
        url_view u("https://www.example.com/get-customer.php?id=409&name=Joe&individual");
        assert( u.query() == "id=409&name=Joe&individual" );
        //]
    }
    {
        //[snippet_parsing_query_1
        url_view u("https://www.example.com/get-customer.php?id=409&name=Joe&individual");
        params_view ps = u.params();
        assert(ps.size() == 3);
        //]
        //[snippet_parsing_query_1a
        auto it = ps.begin();
        assert((*it).key == "id");
        assert((*it).value == "409");
        ++it;
        assert((*it).key == "name");
        assert((*it).value == "Joe");
        ++it;
        assert((*it).key == "individual");
        assert(!(*it).has_value);
        //]
        boost::ignore_unused(it);
    }
    {
        //[snippet_parsing_query_2
        url_view u("https://www.example.com/get-customer.php?key-1=value-1&key-2=&key-3&&=value-2");
        std::cout << u << "\n"
                  "has query:     " << u.has_query()     << "\n"
                  "encoded query: " << u.encoded_query() << "\n"
                  "query:         " << u.query()         << "\n";
        std::cout << u.encoded_params().size() << " parameters\n";
        for (auto p: u.encoded_params())
        {
            if (p.has_value)
            {
                std::cout <<
                    "parameter: <" << p.key <<
                    ", " << p.value << ">\n";
            } else {
                std::cout << "parameter: " << p.key << "\n";
            }
        }
        //]
    }
    {
        //[snippet_parsing_query_3
        url_view u("https://www.example.com/get-customer.php?email=joe@email.com&code=a:2@/!");
        std::cout << u << "\n"
                  "has query:     " << u.has_query()     << "\n"
                  "encoded query: " << u.encoded_query() << "\n"
                  "query:         " << u.query()         << "\n";
        std::cout << u.encoded_params().size() << " parameters\n";
        for (auto p: u.encoded_params())
        {
            if (p.has_value)
            {
                std::cout <<
                    "parameter: <" << p.key <<
                    ", " << p.value << ">\n";
            } else {
                std::cout << "parameter: " << p.key << "\n";
            }
        }
        //]
    }
    {
        //[snippet_parsing_query_4
        url_view u("https://www.example.com/get-customer.php?name=joe");
        std::cout << u << "\n"
                  "query: " << u.query() << "\n";
        //]
    }
    {
        //[snippet_parsing_query_5
        assert(url_view("https://www.example.com/get-customer.php?").has_query());
        assert(!url_view("https://www.example.com/get-customer.php").has_query());
        //]
    }
    {
        //[snippet_parsing_query_6
        url_view u("https://www.example.com/get-customer.php?name=John%20Doe");
        std::cout << u << "\n"
                  "has query:     " << u.has_query()     << "\n"
                  "encoded query: " << u.encoded_query() << "\n"
                  "query:         " << u.query()         << "\n";
        //]
    }
    {
        //[snippet_parsing_query_7
        url_view u("https://www.example.com/get-customer.php?name=John%26Doe");
        assert(u.query() == "name=John&Doe");
        //]
    }
    {
        //[snippet_parsing_query_8
        url_view u("https://www.example.com/get-customer.php?key-1=value-1&key-2=&key-3&&=value-4");
        params_view ps = u.params();
        assert(ps.size() == 5);
        //]
        //[snippet_parsing_query_8a
        auto it = ps.begin();
        assert((*it).key == "key-1");
        assert((*it).value == "value-1");
        //]
        //[snippet_parsing_query_8b
        ++it;
        assert((*it).key == "key-2");
        assert((*it).value == "");
        assert((*it).has_value);
        //]
        //[snippet_parsing_query_8c
        ++it;
        assert((*it).key == "key-3");
        assert(!(*it).has_value);
        //]
        //[snippet_parsing_query_8d
        ++it;
        assert((*it).key == "");
        assert((*it).value == "value-4");
        //]
    }
    {
        //[snippet_parsing_query_9
        url_view u("https://www.example.com/get-customer.php?email=joe@email.com&code=a:2@/!");
        assert(u.has_query());
        //]
    }
}

void
parsing_fragment()
{
    {
        //[snippet_parsing_fragment_1
        url_view u("https://www.example.com/index.html#section%202");
        std::cout << u << "\n"
                  "has fragment:     " << u.has_fragment()     << "\n"
                  "fragment:         " << u.fragment()         << "\n"
                  "encoded fragment: " << u.encoded_fragment() << "\n";
        //]
    }
    {
        //[snippet_parsing_fragment_2_a
        url_view u("https://www.example.com/index.html#");
        std::cout << u << "\n"
                  "has fragment:     " << u.has_fragment()     << "\n"
                  "fragment:         " << u.fragment()         << "\n";
        //]
    }
    {
        //[snippet_parsing_fragment_2_b
        url_view u("https://www.example.com/index.html");
        std::cout << u << "\n"
                  "has fragment:     " << u.has_fragment()     << "\n"
                  "fragment:         " << u.fragment()         << "\n";
        //]
    }
    {
        //[snippet_parsing_fragment_3
        url_view u("https://www.example.com/index.html#code%20:a@b?c/d");
        std::cout << u << "\n"
                  "has fragment:     " << u.has_fragment()     << "\n"
                  "fragment:         " << u.fragment()         << "\n";
        //]
    }
    {
        //[snippet_parsing_fragment_4
        url_view u("https://www.example.com/index.html#section2");
        assert(u.fragment() == "section2");
        //]
    }
    {
        //[snippet_parsing_fragment_5
        assert(url_view("https://www.example.com/index.html#").has_fragment());
        assert(!url_view("https://www.example.com/index.html").has_fragment());
        //]
    }
    {
        //[snippet_parsing_fragment_6
        url_view u("https://www.example.com/index.html#code%20:a@b?c/d");
        assert(u.fragment() == "code :a@b?c/d");
        //]
    }
}

void
using_modifying()
{
    {
        //[snippet_modifying_1
        url_view u("https://www.example.com");
        url v(u);
        //]

        //[snippet_modifying_2
        assert(v.scheme() == "https");
        assert(v.has_authority());
        assert(v.encoded_authority() == "www.example.com");
        assert(v.encoded_path() == "");
        //]

        //[snippet_modifying_3
        v.set_host("my website.com");
        v.set_path("my file.txt");
        v.set_query("id=42&name=John Doe");
        assert(v.buffer() == "https://my%20website.com/my%20file.txt?id=42&name=John%20Doe");
        //]

        //[snippet_modifying_4
        v.set_scheme("http");
        assert(v.buffer() == "http://my%20website.com/my%20file.txt?id=42&name=John%20Doe");
        v.set_encoded_host("www.my%20example.com");
        assert(v.buffer() == "http://my%20example.com/my%20file.txt?id=42&name=John%20Doe");
        //]


    }
}

void
grammar_parse()
{
    {
        //[snippet_parse_1
        // VFALCO we should not show this example
        /*
        boost::core::string_view s = "http:after_scheme";
        const char* it = s.begin();
        auto rv = grammar::parse(it, s.end(), scheme_rule() );
        if( ! rv )
        {
            std::cout << "scheme: " << rv->scheme << '\n';
            std::cout << "suffix: " << it << '\n';
        }
        */
        //]
    }

    {
        //[snippet_parse_2
        // VFALCO This needs refactoring
        /*
        boost::core::string_view s = "?key=value#anchor";
        const char* it = s.begin();
        boost::system::error_code ec;
        if (grammar::parse(it, s.end(), ec, r1))
        {
            auto r2 = grammar::parse( it, s.end(), fragment_part_rule );
            if( r2 )
            {
                std::cout << "query: " << r1.query_part << '\n';
                std::cout << "fragment: " << std::get<1>(*r2.value()).encoded() << '\n';
            }
        }
        */
        //]
    }

    {
        //[snippet_parse_3
        // VFALCO This needs refactoring
        /*
        boost::core::string_view s = "?key=value#anchor";
        query_part_rule r1;
        const char* it = s.begin();
        boost::system::error_code ec;
        auto r2 = grammar::parse( it, s.end(), ec, fragment_part_rule );
        if( ! ec.failed() )
        {
            std::cout << "query: " << r1.query_part << '\n';
            std::cout << "fragment: " << r2.fragment.encoded() << '\n';
        }
        */
        //]
    }

    {
        //[snippet_parse_4
        /* VFALCO This will be removed
        boost::core::string_view s = "http://www.boost.org";
        uri_rule r;
        boost::system::error_code ec;
        if (grammar::parse_string(s, ec, r))
        {
            std::cout << "scheme: " << r.scheme_part.scheme << '\n';
            std::cout << "host: " << r.hier_part.authority.host.host_part << '\n';
        }
        */
        //]
    }
}

//[snippet_customization_1
/* VFALCO This needs rewriting
struct lowercase_rule
{
    boost::core::string_view str;

    friend
    void
    tag_invoke(
        grammar::parse_tag const&,
        char const*& it,
        char const* const end,
        boost::system::error_code& ec,
        lowercase_rule& t) noexcept
    {
        ec = {};
        char const* begin = it;
        while (it != end && std::islower(*it))
        {
            ++it;
        }
        t.str = boost::core::string_view(begin, it);
    }
};
*/
//]

void
grammar_customization()
{
    {
        //[snippet_customization_2
        // VFALCO THIS NEEDS TO BE PORTED
        /*
        boost::core::string_view s = "http:somelowercase";
        scheme_rule r1;
        lowercase_rule r2;
        boost::system::error_code ec;
        if (grammar::parse_string(s, ec, r1, ':', r2))
        {
            std::cout << "scheme: " << r1.scheme << '\n';
            std::cout << "lower:  " << r2.str << '\n';
        }
        */
        //]
    }
}

//[code_charset_1
struct CharSet
{
    bool operator()( char c ) const noexcept;

    // These are both optional. If either or both are left
    // unspecified, a default implementation will be used.
    //
    char const* find_if( char const* first, char const* last ) const noexcept;
    char const* find_if_not( char const* first, char const* last ) const noexcept;
};
//]

void
modifying_path()
{
    {
        //[snippet_modifying_path_1
        url_view u("https://www.boost.org");

        //]
        BOOST_TEST_NOT(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 0u);
    }
    {
        //[snippet_modifying_path_2
        url_view u("https://www.boost.org/");
        //]
        BOOST_TEST(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 0u);
    }
    {
        //[snippet_modifying_path_1_2
        assert( !url_view("https://www.boost.org").is_path_absolute() );
        assert( url_view("https://www.boost.org/").is_path_absolute() );
        //]
    }

    {
        //[snippet_modifying_path_3
        url u("https://www.boost.org/./a/../b");
        u.normalize();
        assert(u.buffer() == "https://www.boost.org/b");
        //]
        BOOST_TEST(u.is_path_absolute());
        BOOST_TEST_EQ(u.buffer(), "https://www.boost.org/b");
        BOOST_TEST_EQ(u.encoded_segments().size(), 1u);
    }
    {
        //[snippet_modifying_path_4
        // scheme and a relative path
        url_view u("https:path/to/file.txt");
        //]
        BOOST_TEST_EQ(u.scheme(), "https");
        BOOST_TEST_NOT(u.has_authority());
        BOOST_TEST_NOT(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 3u);
    }

    {
        //[snippet_modifying_path_5
        // scheme and an absolute path
        url_view u("https:/path/to/file.txt");
        //]
        BOOST_TEST_EQ(u.scheme(), "https");
        BOOST_TEST_NOT(u.has_authority());
        BOOST_TEST(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 3u);
    }

    {
        //[snippet_modifying_path_6
        // "//path" will be considered the authority component
        url_view u("https://path/to/file.txt");
        //]
        BOOST_TEST_EQ(u.scheme(), "https");
        BOOST_TEST(u.has_authority());
        BOOST_TEST(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 2u);
    }

    {
        //[snippet_modifying_path_4_5_6
        // scheme and a relative path
        url_view u1("https:path/to/file.txt");
        assert(!u1.has_authority());
        assert(!u1.is_path_absolute());
        assert(u1.path() == "path/to/file.txt");

        // scheme and an absolute path
        url_view u2("https:/path/to/file.txt");
        assert(!u2.has_authority());
        assert(u2.is_path_absolute());
        assert(u2.path() == "/path/to/file.txt");

        // "//path" will be considered the authority component
        url_view u3("https://path/to/file.txt");
        assert(u3.has_authority());
        assert(u3.authority().buffer() == "path");
        assert(u3.is_path_absolute());
        assert(u3.path() == "/to/file.txt");
        //]
    }

    {
        //[snippet_modifying_path_7
        // only a relative path
        url_view u = parse_uri_reference("path-to/file.txt").value();
        //]
        BOOST_TEST_NOT(u.has_scheme());
        BOOST_TEST_NOT(u.has_authority());
        BOOST_TEST_NOT(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 2u);
    }

    {
        //[snippet_modifying_path_8
        // "path:" will be considered the scheme component
        // instead of a substring of the first segment
        url_view u = parse_uri_reference("path:to/file.txt").value();
        //]
        BOOST_TEST(u.has_scheme());
        BOOST_TEST_NOT(u.has_authority());
        BOOST_TEST_NOT(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 2u);
    }

    {
        //[snippet_modifying_path_7_8
        // only a relative path
        url_view u1 = parse_uri_reference("path-to/file.txt").value();
        assert(!u1.has_scheme());
        assert(u1.path() == "path-to/file.txt");

        // "path:" will be considered the scheme component
        // instead of a substring of the first segment
        url_view u2 = parse_uri_reference("path:to/file.txt").value();
        assert(u2.has_scheme());
        assert(u2.path() == "to/file.txt");
        //]
    }

    {
        //[snippet_modifying_path_9
        // "path" should not become the authority component
        url u = parse_uri("https:path/to/file.txt").value();
        u.set_encoded_path("//path/to/file.txt");
        //]
        BOOST_TEST_EQ(u.scheme(), "https");
        BOOST_TEST_NOT(u.has_authority());
        BOOST_TEST(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 4u);
    }

    {
        //[snippet_modifying_path_10
        // "path:to" should not make the scheme become "path:"
        url u = parse_uri_reference("path-to/file.txt").value();
        u.set_encoded_path("path:to/file.txt");
        //]
        BOOST_TEST_NOT(u.has_scheme());
        BOOST_TEST_NOT(u.has_authority());
        BOOST_TEST_NOT(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 2u);
    }
    {
        //[snippet_modifying_path_9_10
        url u1("https:path/to/file.txt" );
        u1.set_encoded_path("//path/to/file.txt");
        assert(u1.buffer() == "https:/.//path/to/file.txt");
        //]
    }
    {
        //[snippet_modifying_path_11
        // should not insert as "pathto/file.txt"
        url u = parse_uri_reference("to/file.txt").value();
        segments_ref segs = u.segments();
        segs.insert(segs.begin(), "path");
        assert(u.buffer() == "path/to/file.txt");
        //]
        BOOST_TEST_NOT(u.has_scheme());
        BOOST_TEST_NOT(u.has_authority());
        BOOST_TEST_NOT(u.is_path_absolute());
        BOOST_TEST_EQ(u.encoded_segments().size(), 3u);
    }
}

void
normalizing()
{
    {
        //[snippet_normalizing_1
        url_view u1("https://www.boost.org/index.html");
        url_view u2("https://www.boost.org/doc/../index.html");
        assert(u1.buffer() != u2.buffer());
        //]
    }

    {
        //[snippet_normalizing_2
        url_view u1("https://www.boost.org/index.html");
        url_view u2("https://www.boost.org/doc/../index.html");
        assert(u1 == u2);
        //]
    }

    {
        //[snippet_normalizing_3
        url_view u1("https://www.boost.org/index.html");
        url u2("https://www.boost.org/doc/../index.html");
        assert(u1.buffer() != u2.buffer());
        assert(u1 == u2);
        u2.normalize();
        assert(u1.buffer() == u2.buffer());
        assert(u1 == u2);
        //]
    }

    {
        //[snippet_normalizing_4
        url u("https://www.boost.org/doc/../%69%6e%64%65%78%20file.html");
        u.normalize();
        assert(u.buffer() == "https://www.boost.org/index%20file.html");
        //]
    }

    {
        //[snippet_normalizing_5
        auto normalize_http_url =
            [](url& u)
        {
            u.normalize();
            if (u.port() == "80" ||
                u.port().empty())
                u.remove_port();
            if (u.has_authority() &&
                u.encoded_path().empty())
                u.set_path_absolute(true);
        };

        url u1("https://www.boost.org");
        normalize_http_url(u1);
        url u2("https://www.boost.org/");
        normalize_http_url(u2);
        url u3("https://www.boost.org:/");
        normalize_http_url(u3);
        url u4("https://www.boost.org:80/");
        normalize_http_url(u4);

        assert(u1.buffer() == "https://www.boost.org/");
        assert(u2.buffer() == "https://www.boost.org/");
        assert(u3.buffer() == "https://www.boost.org/");
        assert(u4.buffer() == "https://www.boost.org/");
        //]
    }
}

void
decode_with_token()
{
    {
        //[snippet_string_token_1
        url_view u("http://www.example.com/my%20file.txt");
        pct_string_view sv = u.encoded_path();
        assert(sv == "/my%20file.txt");
        std::string s = u.path();
        assert(s == "/my file.txt");
        //]
    }

    {
        //[snippet_string_token_2
        url_view u("http://www.example.com/my%20file.txt");
        std::string s = "existing string";
        u.path(string_token::assign_to(s));
        assert(s == "/my file.txt");
        //]
    }

    {
        //[snippet_string_token_3
        url_view u("http://www.example.com/my%20file.txt");
        std::string s = "existing string";
        u.path(string_token::append_to(s));
        assert(s == "existing string/my file.txt");
        //]
    }

    {
        //[snippet_string_token_4
        url_view u("http://www.example.com/my%20file.txt");
        std::string s = "existing string";
        boost::core::string_view sv = u.path(string_token::preserve_size(s));
        assert(sv == "/my file.txt");
        //]
    }

    {
#if defined(__cpp_static_assert) && __cpp_static_assert >= 201411L
        //[snippet_string_token_5
        static_assert(
            string_token::is_token<string_token::return_string>::value);
        //]
#endif
    }
}

void
encoding()
{
    {
        //[snippet_encoding_1
        std::string s = encode("hello world!", unreserved_chars);
        assert(s == "hello%20world%21");
        //]
    }

    {
        //[snippet_encoding_2
        encoding_opts opt;
        opt.space_as_plus = true;
        std::string s = encode("msg=hello world", pchars, opt);
        assert(s == "msg=hello+world");
        //]
    }

    {
        //[snippet_encoding_3
        std::string s;
        encode("hello ", pchars, {}, string_token::assign_to(s));
        encode("world", pchars, {}, string_token::append_to(s));
        assert(s == "hello%20world");
        //]
    }

    {
        //[snippet_encoding_4
        boost::core::string_view e = "hello world";
        std::string s;
        s.reserve(encoded_size(e, pchars));
        encode(e, pchars, {}, string_token::assign_to(s));
        assert(s == "hello%20world");
        //]
    }

    {
        //[snippet_encoding_5
        boost::core::string_view e = "hello world";
        std::string s;
        s.resize(encoded_size(e, pchars));
        encode(&s[0], s.size(), e, pchars);
        assert(s == "hello%20world");
        //]
    }

    {
        //[snippet_encoding_6
        pct_string_view sv = "hello%20world";
        assert(sv == "hello%20world");
        //]
    }

    {
        //[snippet_encoding_7
        boost::system::result<pct_string_view> rs =
            make_pct_string_view("hello%20world");
        assert(rs.has_value());
        pct_string_view sv = rs.value();
        assert(sv == "hello%20world");
        //]
    }

    {
        //[snippet_encoding_8
        pct_string_view s = "path/to/file";
        url u;
        u.set_encoded_path(s);
        assert(u.buffer() == "path/to/file");
        //]
    }

    {
        //[snippet_encoding_9
        url u;
        u.set_encoded_path("path/to/file");
        assert(u.buffer() == "path/to/file");
        //]
    }

    {
        //[snippet_encoding_10
        url_view uv("path/to/file");
        url u;
        u.set_encoded_path(uv.encoded_path());
        assert(u.buffer() == "path/to/file");
        //]
    }

    {
        //[snippet_encoding_11
        pct_string_view es("hello%20world");
        assert(es == "hello%20world");

        decode_view dv("hello%20world");
        assert(dv == "hello world");
        //]
    }

    {
        //[snippet_encoding_12
        boost::system::result<pct_string_view> rs =
            make_pct_string_view("hello%20world");
        assert(rs.has_value());
        pct_string_view s = rs.value();
        decode_view dv = *s;
        assert(dv == "hello world");
        //]
    }

    {
        //[snippet_encoding_13
        url_view u =
            parse_relative_ref("user/john%20doe/profile%20photo.jpg").value();
        std::vector<std::string> route =
            {"user", "john doe", "profile photo.jpg"};
        auto segs = u.encoded_segments();
        auto it0 = segs.begin();
        auto end0 = segs.end();
        auto it1 = route.begin();
        auto end1 = route.end();
        while (
            it0 != end0 &&
            it1 != end1)
        {
            pct_string_view seg0 = *it0;
            decode_view dseg0 = *seg0;
            boost::core::string_view seg1 = *it1;
            if (dseg0 == seg1)
            {
                ++it0;
                ++it1;
            }
            else
            {
                break;
            }
        }
        bool route_match = it0 == end0 && it1 == end1;
        assert(route_match);
        //]
    }

    {
        //[snippet_encoding_14
        pct_string_view s = "user/john%20doe/profile%20photo.jpg";
        std::string buf;
        buf.resize(s.decoded_size());
        s.decode({}, string_token::assign_to(buf));
        assert(buf == "user/john doe/profile photo.jpg");
        //]
    }
}

void
readme_snippets()
{
    // Parse a URL. This allocates no memory. The view
    // references the character buffer without taking ownership.
    //
    url_view uv( "https://www.example.com/path/to/file.txt?id=1001&name=John%20Doe&results=full" );

    // Print the query parameters with percent-decoding applied
    //
    for( auto v : uv.params() )
    {
        std::cout << v.key << "=" << v.value << " ";
    }

    // Prints: id=1001 name=John Doe results=full

    // Create a modifiable copy of `uv`, with ownership of the buffer
    //
    url u = uv;

    // Change some elements in the URL
    //
    u.set_scheme( "http" )
        .set_encoded_host( "boost.org" )
        .set_encoded_path( "/index.htm" )
        .remove_query()
        .remove_fragment()
        .params().append( {"key", "value"} );

    std::cout << u;
}

//[snippet_using_static_pool_1
// VFALCO NOPE
//]

namespace boost {
namespace urls {

class snippets_test
{
public:
    void
    run()
    {
        ignore_unused(&using_url_views);
        ignore_unused(&using_urls);
        // parsing_urls();
        parsing_components();
        formatting_components();
        ignore_unused(&parsing_scheme);
        ignore_unused(&parsing_authority);
        ignore_unused(&parsing_path);
        ignore_unused(&parsing_query);
        ignore_unused(&parsing_fragment);
        ignore_unused(&using_modifying);
        ignore_unused(&grammar_parse);
        ignore_unused(&grammar_customization);
        ignore_unused(&modifying_path);
        normalizing();
        decode_with_token();
        encoding();
        ignore_unused(&readme_snippets);

        BOOST_TEST_PASS();
    }
};

TEST_SUITE(snippets_test, "boost.url.snippets");

} // urls
} // boost
