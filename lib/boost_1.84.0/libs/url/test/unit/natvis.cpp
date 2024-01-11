//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#include <boost/url/error.hpp>
#include <boost/url/string_view.hpp>
#include <boost/assert/source_location.hpp>
#include <boost/core/ignore_unused.hpp>
#include <system_error>
#include "test_suite.hpp"

namespace boost {
namespace urls {

struct natvis_test
{
    struct yesexcept
    {
        int id;
        yesexcept()
            : id([]
                {
                    static int id_ = 0;
                    return ++id_;
                }())
        {
        }
        yesexcept(yesexcept&& u) { id = u.id; }
        yesexcept(yesexcept const& u) { id = u.id; }
        yesexcept& operator=(yesexcept&& u) { id = u.id; return *this; }
        yesexcept& operator=(yesexcept const& u) { id = u.id; return *this; }
    };

    struct my_category : boost::system::error_category
    {
        my_category() noexcept
            : boost::system::error_category(0xabadfadeadeadfad)
        {
        }

        const char* name() const noexcept override
        {
            return "boost.url.natvis";
        }

        std::string message(int) const override
        {
            return {};
        }

        boost::system::error_condition default_error_condition(
            int ev) const noexcept override
        {
            return {ev, *this};
        }
    };

    // these are here to view the results of
    // .natvis definitions in the debugger.
    void
    run()
    {
        // boost::assert::source_location
        {
            static auto loc = BOOST_CURRENT_LOCATION;
            ignore_unused(loc);
        }

        // boost::variant2::variant
        {
        }

        // boost::system::error_category
        {
            auto const& c1 = boost::system::generic_category();
            auto const& c2 = boost::system::system_category();
            auto const& c3 = system::error_code(std::error_code()).category();
            auto const& c4 = my_category();
            auto const& c5 = system::error_code(error::not_a_base).category();
            ignore_unused(c1, c2, c3, c4, c5);
        }

        // boost::system::error_code
        {
            static auto loc = BOOST_CURRENT_LOCATION;
            auto const e0 = system::error_code();
            auto const e1 = system::error_code(std::make_error_code(std::errc::address_in_use));
            auto const e2 = system::error_code(error::success);
            auto const e3a = system::error_code(error::not_a_base);
            auto const e3b = system::error_code(make_error_code(boost::system::errc::bad_address));
            auto const e4 = system::error_code(system::error_code(error::success), &loc);
            auto const e5 = system::error_code(system::error_code(error::not_a_base), &loc);
            ignore_unused(e0, e1, e2, e3a, e3b, e4, e5);
        }

        // boost::system::result
        {
            system::result<double> rv1;
            system::result<double> rv2 = 3.14;
            system::result<double> rv3 = error::not_a_base;
            system::result<yesexcept> rv4;
            system::result<yesexcept> rv5 = yesexcept{};
            system::result<yesexcept> rv6 = error::not_a_base;
            ignore_unused(rv1, rv2, rv3, rv4, rv5, rv6);
        }

        // boost::core::core::string_view
        {
            core::string_view s1;
            core::string_view s2 = "This is how we do it";
            core::string_view s3 = s2.substr(8, 3);
            ignore_unused(s1, s2, s3);
        }
    }
};

TEST_SUITE(
    natvis_test,
    "boost.url.natvis");

} // urls
} // boost

