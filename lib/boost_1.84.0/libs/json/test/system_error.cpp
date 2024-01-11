//
// Copyright (c) 2022 Dmitry Arkhipov (grisumbras@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/system_error.hpp>
// Test that header file is properly guarder
#include <boost/json/system_error.hpp>

#include <boost/json/value_to.hpp>

#include "test_suite.hpp"

namespace boost {
namespace json {

class system_error_test
{
public:
    void
    run()
    {
        BOOST_STATIC_CONSTEXPR source_location loc = BOOST_CURRENT_LOCATION;
        auto const res = result_from_errno<int>(EDOM, &loc);
        BOOST_TEST( res.has_error() );

        auto const ec = res.error();
        BOOST_TEST( ec.has_location() );
        BOOST_TEST( ec.location() == loc );
        BOOST_TEST( ec.category() == system::generic_category() );
        BOOST_TEST( ec == std::errc::argument_out_of_domain );
    }
};

TEST_SUITE(system_error_test, "boost.json.system_error");

} // namespace json
} // namespace boost
