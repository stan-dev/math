// Copyright 2018-2023 Emil Dotchevski and Reverge Studios, Inc.

// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/handle_errors.hpp>
#   include <boost/leaf/result.hpp>
#   include <boost/leaf/pred.hpp>
#endif

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

enum class ErrorCode
{
    E_GENERIC_UNEXPECTED,
    E_GENERIC_PARSE
};

leaf::result<int> test1( int a )
{
    if( a == 1 )
        return leaf::new_error(ErrorCode::E_GENERIC_UNEXPECTED);
    return 32;
}

leaf::result<float> test2(int a)
{
    return leaf::try_handle_some(
        [&]() -> leaf::result<float>
        {
            BOOST_LEAF_AUTO(val, test1(a));
            (void) val;
            return 4.5;
        },
        [](leaf::match<ErrorCode,ErrorCode::E_GENERIC_UNEXPECTED>) -> leaf::result<float>
        {
            return leaf::new_error(ErrorCode::E_GENERIC_PARSE);
        }
    );
}

void test3(int a)
{
    int x = 0;
    leaf::try_handle_all(
        [&]() -> leaf::result<void>
        {
            BOOST_LEAF_AUTO(val, test2(a));
            (void) val;
            return {};
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
