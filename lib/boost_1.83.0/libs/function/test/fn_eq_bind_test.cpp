// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/function.hpp>
#include <boost/bind/bind.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::placeholders;

int f1() { return 1; }
int f2() { return 2; }

int main()
{
    {
        boost::function<int()> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int, int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int, int, int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int, int, int, int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int, int, int, int, int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int, int, int, int, int, int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int, int, int, int, int, int, int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int, int, int, int, int, int, int, int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    {
        boost::function<int(int, int, int, int, int, int, int, int, int)> fn( boost::bind( f1 ) );

        BOOST_TEST( fn == boost::bind( f1 ) );
        BOOST_TEST( fn != boost::bind( f2 ) );
    }

    return boost::report_errors();
}
