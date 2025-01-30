//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#include <boost/url.hpp>
#include "test_suite.hpp"

#include <boost/core/ignore_unused.hpp>
#include <iostream>

namespace boost {
namespace urls {
namespace grammar {

namespace impl {

// tag::code_grammar_1_1[]
template< class Rule >
auto parse( core::string_view s, Rule const& r) -> system::result< typename Rule::value_type >;

template< class Rule >
auto parse( char const *& it, char const* end, Rule const& r) -> system::result< typename Rule::value_type >;
// end::code_grammar_1_1[]

} // impl

//------------------------------------------------

namespace {

// tag::code_grammar_2_2[]
struct ws_chars_t
{
    constexpr bool operator()( char c ) const noexcept
    {
        return c == '\t' || c == ' ' || c == '\r' || c == '\n';
    }
};
// end::code_grammar_2_2[]

// tag::code_grammar_2_3[]
static_assert( is_charset< ws_chars_t >::value, "CharSet requirements not met" );
// end::code_grammar_2_3[]

// tag::code_grammar_2_4[]
constexpr ws_chars_t ws_chars{};
// end::code_grammar_2_4[]

// tag::code_grammar_2_5[]
core::string_view get_token( core::string_view s ) noexcept
{
    auto it0 = s.data();
    auto const end = it0 + s.size();

    // find the first non-whitespace character
    it0 = find_if_not( it0, end, ws_chars );

    if( it0 == end )
    {
        // all whitespace or empty string
        return {};
    }

    // find the next whitespace character
    auto it1 = find_if( it0, end, ws_chars );

    // [it0, it1) is the part we want
    return core::string_view( it0, it1 - it0 );
}
// end::code_grammar_2_5[]

// tag::code_grammar_2_9[]
struct is_visible
{
    constexpr bool operator()( char ch ) const noexcept
    {
        return ch >= 33 && ch <= 126;
    }
};
constexpr lut_chars visible_chars( is_visible{} ); // (since C++11)
// end::code_grammar_2_9[]

#if __cplusplus >= 201703L
namespace cpp17 {
// tag::code_grammar_2_10[]
constexpr lut_chars visible_chars( [](char ch) { return ch >= 33 && ch <= 126; } ); // (since C++17)
// end::code_grammar_2_10[]
} // cpp17
#endif

} // (anon)

struct doc_grammar_test
{
    void snippets1()
    {
        // tag::code_grammar_1_2[]
        struct comma_rule_t
        {
            // The type of value returned upon success
            using value_type = core::string_view;

            // The algorithm which checks for a match
            system::result< value_type >
            parse( char const*& it, char const* end ) const
            {
                if( it != end && *it == ',')
                    return core::string_view( it++, 1 );

                return error::mismatch;
            }
        };
        // end::code_grammar_1_2[]

        // tag::code_grammar_1_3[]
        constexpr comma_rule_t comma_rule{};
        // end::code_grammar_1_3[]

        {
        // tag::code_grammar_1_4[]
        system::result< core::string_view > rv = parse( ",", comma_rule );

        assert( rv.has_value() && rv.value() == "," );
        // end::code_grammar_1_4[]
        }
        {
        // tag::code_grammar_1_5[]
        system::result< unsigned short > rv =
            parse( "16384", unsigned_rule< unsigned short >{} );
        // end::code_grammar_1_5[]
        ignore_unused(rv);
        }
        {
        // tag::code_grammar_1_6[]
        system::result< core::string_view > rv = parse( ",", delim_rule(',') );
        // end::code_grammar_1_6[]
        ignore_unused(rv);
        }
    }
    //--------------------------------------------

    void snippets2()
    {
        {
        // tag::code_grammar_2_6[]
        assert( get_token( " \t john-doe\r\n \t jane-doe\r\n") == "john-doe" );
        // end::code_grammar_2_6[]
        }
        static
        // tag::code_grammar_2_7[]
        constexpr lut_chars vowels = "AEIOU" "aeiou";
        // end::code_grammar_2_7[]
        {
        // tag::code_grammar_2_8[]
        constexpr auto vowels_and_y = vowels + 'y' + 'Y';
        // end::code_grammar_2_8[]
        ignore_unused(vowels_and_y);
        }
        {
        // tag::code_grammar_2_11[]
        constexpr auto visible_non_vowels = visible_chars - vowels;
        // end::code_grammar_2_11[]
        ignore_unused(visible_non_vowels);
        }
        {
        // tag::code_grammar_2_12[]
        constexpr auto visible_non_vowels_or_y = visible_chars - vowels - 'y';
        // end::code_grammar_2_12[]
        ignore_unused(visible_non_vowels_or_y);
        }
    }

    //--------------------------------------------

    void snippets3()
    {
        {
        // tag::code_grammar_3_1[]
        constexpr auto version_rule = tuple_rule( delim_rule( 'v' ), dec_octet_rule, delim_rule( '.' ), dec_octet_rule );
        // end::code_grammar_3_1[]

        // tag::code_grammar_3_2[]
        system::result< std::tuple< core::string_view, unsigned char, core::string_view, unsigned char > > rv = parse( "v42.44800", version_rule );
        // end::code_grammar_3_2[]
        ignore_unused(rv);
        }
        {
        // tag::code_grammar_3_3[]
        constexpr auto version_rule = tuple_rule( squelch( delim_rule( 'v' ) ), dec_octet_rule, squelch( delim_rule( '.' ) ), dec_octet_rule );

        system::result< std::tuple< unsigned char, unsigned char > > rv = parse( "v42.44800", version_rule );
        // end::code_grammar_3_3[]
        ignore_unused(rv);
        }
        {
        // tag::code_grammar_3_4[]
        // port     = ":" unsigned-short

        constexpr auto port_rule = tuple_rule( squelch( delim_rule( ':' ) ), unsigned_rule< unsigned short >{} );

        system::result< unsigned short > rv = parse( ":443", port_rule );
        // end::code_grammar_3_4[]
        ignore_unused(rv);
        }
        {
        // tag::code_grammar_3_5[]
        // port     = [ ":" unsigned-short ]

        constexpr auto port_rule = optional_rule( tuple_rule( squelch( delim_rule( ':' ) ), unsigned_rule< unsigned short >{} ) );

        system::result< boost::optional< unsigned short > > rv = parse( ":8080", port_rule );

        assert( rv->has_value() && rv->value() == 8080 );
        // end::code_grammar_3_5[]
        }
        {
        // tag::code_grammar_3_6[]
        // ipv4_address = dec-octet "." dec-octet "." dec-octet "." dec-octet
        //
        // port         = ":" unsigned-short
        //
        // endpoint     = ipv4_address [ port ]

        constexpr auto endpoint_rule = tuple_rule(
            tuple_rule(
                dec_octet_rule, squelch( delim_rule( '.' ) ),
                dec_octet_rule, squelch( delim_rule( '.' ) ),
                dec_octet_rule, squelch( delim_rule( '.' ) ),
                dec_octet_rule ),
            optional_rule(
                tuple_rule(
                    squelch( delim_rule( ':' ) ),
                    unsigned_rule< unsigned short >{} ) ) );
        // end::code_grammar_3_6[]
        ignore_unused(endpoint_rule);
        }
        {
        // tag::code_grammar_3_7[]
        constexpr auto endpoint_rule = tuple_rule(
            ipv4_address_rule,
            optional_rule(
                tuple_rule(
                    squelch( delim_rule( ':' ) ),
                    unsigned_rule< unsigned short >{} ) ) );

        system::result< std::tuple< ipv4_address, boost::optional< unsigned short > > > rv = parse( "192.168.0.1:443", endpoint_rule );
        // end::code_grammar_3_7[]
        ignore_unused(rv);
        }
        {
        // tag::code_grammar_3_8[]
        constexpr auto request_target_rule = variant_rule(
            origin_form_rule,
            absolute_uri_rule,
            authority_rule,
            delim_rule('*') );

        system::result< variant2::variant< url_view, url_view, authority_view, core::string_view > > rv = parse( "/results.htm?page=4", request_target_rule );
        // end::code_grammar_3_8[]
        }
    }

    void snippets4()
    {
        {
        // tag::code_grammar_4_1[]
        constexpr auto chunk_ext_rule = range_rule(
            tuple_rule( squelch( delim_rule( ';' ) ), token_rule( alnum_chars ) ) );
        // end::code_grammar_4_1[]
        // tag::code_grammar_4_2[]
        system::result< range< core::string_view > > rv = parse( ";johndoe;janedoe;end", chunk_ext_rule );

        for( auto s : rv.value() )
            std::cout << s << "\n";
        // end::code_grammar_4_2[]
        }
        {
        // tag::code_grammar_4_3[]
        constexpr auto token_list_rule = range_rule(
            token_rule( alnum_chars ),
            tuple_rule( squelch( delim_rule( ',' ) ), token_rule( alnum_chars ) ),
            1 );
        // end::code_grammar_4_3[]
        // tag::code_grammar_4_4[]
        system::result< range< core::string_view > > rv = parse( "johndoe,janedoe,end", token_list_rule );

        for( auto s : rv.value() )
            std::cout << s << "\n";
        // end::code_grammar_4_4[]
        }
    }

    void
    snippets5()
    {
        {
        core::string_view s;
        // tag::code_grammar_5_1[]
        system::result< pct_string_view > rv = parse( s, pct_encoded_rule( pchars ) );
        // end::code_grammar_5_1[]
        ignore_unused(rv);
        }
        {
        // tag::code_grammar_5_2[]
        // request-line   = method SP request-target SP HTTP-version CRLF

        constexpr auto request_line_rule = tuple_rule(
            not_empty_rule( token_rule( alpha_chars ) ),    // method
            squelch( delim_rule( ' ' ) ),                   // SP
            variant_rule(
                absolute_uri_rule,                          // absolute-uri or
                relative_ref_rule),                         // relative-ref
            squelch( delim_rule( ' ' ) ),
            squelch( literal_rule( "HTTP/" ) ),             // "HTTP/"
            delim_rule( digit_chars ),                      // DIGIT
            squelch( delim_rule( '.' ) ),                   // "."
            delim_rule( digit_chars ),                      // DIGIT
            squelch( literal_rule( "\r\n" ) ) );            // CRLF
        // end::code_grammar_5_2[]
        ignore_unused(request_line_rule);
        }
    }

    void
    run()
    {
        snippets1();
        snippets2();
        snippets3();
        snippets4();
    }
};

TEST_SUITE(
    doc_grammar_test,
    "boost.url.doc.grammar");

} // grammar
} // urls
} // boost