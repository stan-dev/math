//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/encode.hpp>

#include <boost/url/rfc/pchars.hpp>
#include <boost/core/ignore_unused.hpp>

#include "test_suite.hpp"

#include <memory>

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

namespace boost {
namespace urls {

class encode_test
{
public:
    struct test_chars
    {
        constexpr
        bool
        operator()(char c) const noexcept
        {
            return c == 'A' || c == '+';
        }
    };

    void
    check(
        core::string_view s,
        core::string_view m0,
        bool space_as_plus = false)
    {
        // encoded_size
        {
            encoding_opts opt;
            opt.space_as_plus =
                space_as_plus;
            BOOST_TEST(encoded_size(
                s, test_chars{}, opt) ==
                    m0.size());
        }
        // encode
        {
            encoding_opts opt;
            opt.space_as_plus =
                space_as_plus;
            std::string t;
            t.resize(
                encoded_size(s, test_chars{}, opt));
            encode(
                &t[0], t.size(), s, test_chars{}, opt);
            BOOST_TEST(t == m0);
        }
        encoding_opts opt;
        opt.space_as_plus =
            space_as_plus;
        auto const m = encode(
            s, test_chars{}, opt, {});
        if(! BOOST_TEST(m == m0))
            return;
        char buf[64];
        BOOST_ASSERT(
            m.size() < sizeof(buf));
        for(std::size_t i = 0;
            i <= sizeof(buf); ++i)
        {
            char* dest = buf;
            std::size_t n = encode(
                dest, i, s, test_chars{}, opt);
            core::string_view r(buf, n);
            if(n == m.size())
            {
                BOOST_TEST_EQ(i, m.size());
                BOOST_TEST_EQ(r, m);
                break;
            }
            BOOST_TEST(
                core::string_view(buf, n) ==
                m.substr(0, n));
        }
    };

    void
    testEncode()
    {
        check("", "");
        check(" ", "%20");
        check("A", "A");
        check("B", "%42");
        check("AB", "A%42");
        check("A B", "A%20%42");

        check("", "", true);
        check(" ", "+", true);
        check("A", "A", true);
        check("B", "%42", true);
        check("AB", "A%42", true);
        check("A B", "A+%42", true);
    }

    void
    testEncodeExtras()
    {
        // space_as_plus
        {
            BOOST_TEST(encode(
                " ", test_chars{}, {}, {}) == "%20");
            encoding_opts opt;
            BOOST_TEST_EQ(opt.space_as_plus, false);
            BOOST_TEST(encode(
                " ", test_chars{}, opt, {}) == "%20");
            BOOST_TEST(encode(
                "A", test_chars{}, opt, {}) == "A");
            BOOST_TEST(encode(
                " A+", test_chars{}, opt, {}) == "%20A+");
            opt.space_as_plus = true;
            BOOST_TEST(encode(
                " ", test_chars{}, opt, {}) == "+");
            BOOST_TEST(encode(
                "A", test_chars{}, opt, {}) == "A");
            BOOST_TEST(encode(
                " A+", test_chars{}, opt, {}) == "+A+");
        }
    }

    void
    testJavadocs()
    {
        // encoded_size()
        {
    assert( encoded_size( "My Stuff", pchars ) == 10 );
        }

        // encode()
        {
    char buf[100];
    assert( encode( buf, sizeof(buf), "Program Files", pchars ) == 15 );

    ignore_unused(buf);
        }
    }

    void
    run()
    {
        testEncode();
        testEncodeExtras();
        testJavadocs();
    }
};

TEST_SUITE(
    encode_test,
    "boost.url.encode");

} // urls
} // boost
