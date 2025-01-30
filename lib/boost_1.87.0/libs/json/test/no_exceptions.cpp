//
// Copyright (c) 2023 Dmitry Arkhipov (grisumbras@yandex.ru)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

#include <boost/json.hpp>

#include <cstdio>

#include "test_suite.hpp"

namespace boost {

void
throw_exception( std::exception const&, boost::source_location const&)
{
    std::printf("Exceptions are not supported!");
    std::abort();
}

namespace json {

class no_exceptions
{
public:
    void
    run()
    {
        system::error_code ec;
        value jv = parse(" 1 ", ec);
        BOOST_TEST( !ec.failed() );
        BOOST_TEST( jv == 1 );

        jv = parse(" x ", ec);
        BOOST_TEST( ec == error::syntax );
    }
};

TEST_SUITE(no_exceptions, "boost.json.no_exceptions");

} // namespace json
} // namespace boost
