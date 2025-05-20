// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <boost/leaf/config.hpp>

#ifdef BOOST_LEAF_NO_EXCEPTIONS

#include <iostream>

int main()
{
    std::cout << "Unit test not applicable." << std::endl;
    return 0;
}

#else

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/handle_errors.hpp>
#   include <boost/leaf/exception.hpp>
#   include <boost/leaf/pred.hpp>
#endif

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

enum class ErrorCode
{
    E_GENERIC_UNEXPECTED,
    E_GENERIC_PARSE
};

int test1( int a )
{
    if( a == 1 )
        BOOST_LEAF_THROW_EXCEPTION(ErrorCode::E_GENERIC_UNEXPECTED);
    return 32;
}

float test2(int a)
{
    return leaf::try_catch(
        [&]
        {
            auto val = test1(a);
            (void) val;
            return 4.5;
        },
        [](leaf::match<ErrorCode,ErrorCode::E_GENERIC_UNEXPECTED>)
        {
            BOOST_LEAF_THROW_EXCEPTION(ErrorCode::E_GENERIC_PARSE);
            return 0;
        }
    );
}

void test3(int a)
{
    int x = 0;
    leaf::try_catch(
        [&]
        {
            auto val = test2(a);
            (void) val;
        },
        [&](leaf::match<ErrorCode, ErrorCode::E_GENERIC_PARSE>)
        {
            x = 1;
        },
        [&](ErrorCode e)
        {
            x = 2;
        },
        [&]()
        {
            x = 3;
        }
    );
    BOOST_TEST_EQ(x, 1);
}

int main()
{
    test3(1);
    return boost::report_errors();
}

#endif
