//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/unsigned_rule.hpp>

#include <boost/url/grammar/alpha_chars.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/token_rule.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct unsigned_rule_test
{
    void
    run()
    {
        // constexpr
        {
            constexpr auto r =
                unsigned_rule<unsigned short>{};
            (void)r;
        }

        // javadoc
        {
            auto rv = parse( "32767", unsigned_rule< unsigned short >{} );

            (void)rv;
        }

        {
            using T = std::uint8_t;
            constexpr auto t =
                unsigned_rule<T>{};

            ok(t, "0", T(0));
            ok(t, "1", T(1));
            ok(t, "9", T(9));
            ok(t, "255", T(255));

            bad(t, "00", error::invalid);
            bad(t, "01", error::invalid);
            bad(t, "256", error::invalid);
            bad(t, "300", error::invalid);
            bad(t, "2555", error::invalid);
            bad(t, "25555", error::invalid);
        }
        {
            using T = std::uint16_t;
            constexpr auto t =
                unsigned_rule<T>{};

            ok(t, "0", T(0));
            ok(t, "1", T(1));
            ok(t, "99", T(99));
            ok(t, "65535", T(65535));

            bad(t, "", error::mismatch);
            bad(t, "a", error::mismatch);
            bad(t, "00", error::invalid);
            bad(t, "01", error::invalid);
            bad(t, "65536", error::invalid);
            bad(t, "70000", error::invalid);
        }
        {
            using T = std::uint32_t;
            constexpr auto t =
                unsigned_rule<T>{};

            ok(t, "0", T(0));
            ok(t, "1", T(1));
            ok(t, "999", T(999));
            ok(t, "4294967295", T(4294967295));

            bad(t, "00", error::invalid);
            bad(t, "01", error::invalid);
            bad(t, "4294967296", error::invalid);
            bad(t, "5000000000", error::invalid);
        }
    }
};

TEST_SUITE(
    unsigned_rule_test,
    "boost.url.grammar.unsigned_rule");

} // grammar
} // urls
} // boost
