// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/leaf/config.hpp>

#if defined(BOOST_LEAF_NO_EXCEPTIONS) || !BOOST_LEAF_CFG_CAPTURE

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
#   include <boost/leaf/on_error.hpp>
#   include <boost/leaf/exception.hpp>
#endif

#include "_test_ec.hpp"
#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

template <int> struct info { int value; };

template <class F>
void test( F f )
{
    {
        int c=0;
        auto r = f();
        leaf::try_catch(
            [&r]
            {
                (void) r.value();
            },
            [&c]( info<1> const & x )
            {
                BOOST_TEST_EQ(x.value, 1);
                BOOST_TEST_EQ(c, 0);
                c = 1;
            },
            [&c]
            {
                BOOST_TEST_EQ(c, 0);
                c = 2;
            } );
        BOOST_TEST_EQ(c, 1);
    }

    {
        int c=0;
        auto r = f();
        leaf::try_catch(
            [&r]
            {
                (void) r.value();
            },
            [&c]( info<2> const & x )
            {
                BOOST_TEST_EQ(x.value, 2);
                BOOST_TEST_EQ(c, 0);
                c = 1;
            },
            [&c]
            {
                BOOST_TEST_EQ(c, 0);
                c = 2;
            } );
        BOOST_TEST_EQ(c, 2);
    }

    {
        auto r = f();
        int what = leaf::try_catch(
            [&r]
            {
                (void) r.value();
                return 0;
            },
            []( info<1> const & x )
            {
                BOOST_TEST_EQ(x.value, 1);
                return 1;
            },
            []
            {
                return 2;
            } );
        BOOST_TEST_EQ(what, 1);
    }

    {
        auto r = f();
        int what = leaf::try_catch(
            [&r]
            {
                (void) r.value();
                return 0;
            },
            []( info<2> const & x )
            {
                BOOST_TEST_EQ(x.value, 2);
                return 1;
            },
            []
            {
                return 2;
            } );
        BOOST_TEST_EQ(what, 2);
    }
}

int main()
{
    test( []
    {
        return leaf::try_capture_all(
            []() -> int
            {
                leaf::throw_exception(errc_a::a0, info<1>{1}, info<3>{3});
            } );
    } );

    test( []
    {
        return leaf::try_capture_all(
            []() -> void
            {
                leaf::throw_exception(errc_a::a0, info<1>{1}, info<3>{3});
            } );
    } );

    test( []
    {
        return leaf::try_capture_all(
            []() -> int
            {
                auto load = leaf::on_error(errc_a::a0, info<1>{1}, info<3>{3});
                leaf::throw_exception();
            } );
    } );

    test( []
    {
        return leaf::try_capture_all(
            []() -> void
            {
                auto load = leaf::on_error(errc_a::a0, info<1>{1}, info<3>{3});
                leaf::throw_exception();
            } );
    } );

    test( []
    {
        return leaf::try_capture_all(
            []() -> int
            {
                auto load = leaf::on_error(errc_a::a0, info<1>{1}, info<3>{3});
                throw std::exception();
            } );
    } );

    test( []
    {
        return leaf::try_capture_all(
            []() -> void
            {
                auto load = leaf::on_error(errc_a::a0, info<1>{1}, info<3>{3});
                throw std::exception();
            } );
    } );

    return boost::report_errors();
}

#endif
