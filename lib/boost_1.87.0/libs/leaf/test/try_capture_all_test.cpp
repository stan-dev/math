// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/leaf/config.hpp>

#if BOOST_LEAF_CFG_STD_STRING
#   include <sstream>
#   include <iostream>
#endif

#if !BOOST_LEAF_CFG_CAPTURE

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
#   include <boost/leaf/result.hpp>
#   include <boost/leaf/handle_errors.hpp>
#endif

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

namespace
{
    template <int>
    struct err
    {
        static int count;
        static int move_count;

        err()
        {
            ++count;
        }

        err( err const & )
        {
            ++count;
        }

        err( err && )
        {
            ++count;
            ++move_count;
        }

        ~err()
        {
            --count;
        }
    };
    template <int N> int err<N>::count = 0;
    template <int N> int err<N>::move_count = 0;
}

int main()
{
#ifndef BOOST_LEAF_NO_EXCEPTIONS
    {
        leaf::result<void> r = leaf::try_capture_all([]() { throw std::runtime_error("x"); });
        BOOST_TEST(!r);
        int r1 = leaf::try_handle_all(
            [&]() -> leaf::result<int>
            {
                BOOST_LEAF_CHECK(r);
                return 0;
            },
            [](std::runtime_error const &)
            {
                return 1;
            },
            []
            {
                return 2;
            } );
        BOOST_TEST_EQ(r1, 1);
    }
#endif
    {
        leaf::result<void> r = leaf::try_capture_all(
            []() -> leaf::result<void>
            {
                return leaf::new_error(err<1>{}, err<2>{});
            });
        BOOST_TEST_EQ(err<1>::count, 1);
        BOOST_TEST_EQ(err<2>::count, 1);
        BOOST_TEST(!r);
        int r1 = leaf::try_handle_all(
            [&]() -> leaf::result<int>
            {
                BOOST_LEAF_CHECK(r);
                return 0;
            },
            [](err<1>, err<2>)
            {
                return 1;
            },
            []
            {
                return 2;
            } );
        BOOST_TEST_EQ(r1, 1);
    }
    BOOST_TEST_EQ(err<1>::count, 0);
    BOOST_TEST_EQ(err<2>::count, 0);

    {
        leaf::result<int> r = leaf::try_capture_all(
            []() -> leaf::result<int>
            {
                return leaf::new_error(err<1>{}, err<2>{});
            });
        BOOST_TEST(!r);
        BOOST_TEST_EQ(err<1>::count, 1);
        BOOST_TEST_EQ(err<2>::count, 1);
        int r1 = leaf::try_handle_all(
            [&]() -> leaf::result<int>
            {
                BOOST_LEAF_CHECK(r);
                return 0;
            },
            [](err<1>, err<2>)
            {
                return 1;
            },
            []
            {
                return 2;
            } );
        BOOST_TEST_EQ(r1, 1);
    }
    BOOST_TEST_EQ(err<1>::count, 0);
    BOOST_TEST_EQ(err<2>::count, 0);

    return boost::report_errors();
}

#endif
