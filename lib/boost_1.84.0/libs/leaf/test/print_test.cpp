// Copyright 2018-2023 Emil Dotchevski and Reverge Studios, Inc.

// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/leaf/config.hpp>

#if !BOOST_LEAF_CFG_DIAGNOSTICS && BOOST_LEAF_CFG_STD_STRING

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
#   include <boost/leaf/detail/print.hpp>
#endif

#include <sstream>
#include <iostream>

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

struct c0
{
    friend std::ostream & operator<<( std::ostream & os, c0 const & )
    {
        return os << "c0";
    }
};

struct c1
{
    int value;

    friend std::ostream & operator<<( std::ostream & os, c1 const & )
    {
        return os << "c1";
    }
};

struct c2
{
    int value;
};

std::ostream & operator<<( std::ostream & os, c2 const & )
{
    return os << "c2";
}

struct c3
{
    int value;
};

struct c4
{
    struct unprintable { };
    unprintable value;;
};

template <int Line, class T>
bool check( T const & x, char const * sub )
{
    using namespace leaf::leaf_detail;
    std::ostringstream s;
    diagnostic<T>::print(s,x);
    std::string q = s.str();
    std::cout << "LINE " << Line << ": " << q << std::endl;
    return q.find(sub)!=q.npos;
}

struct my_exception: std::exception
{
    char const * what() const noexcept override { return "my_exception_what"; }
};

int main()
{
    BOOST_TEST(check<__LINE__>(c0{ },"c0"));
    BOOST_TEST(check<__LINE__>(c1{42},"c1"));
    {
        c1 x;
        c1 & y = x;
        BOOST_TEST(check<__LINE__>(x,"c1"));
        BOOST_TEST(check<__LINE__>(y,"c1"));
    }
    BOOST_TEST(check<__LINE__>(c2{42},"c2"));
    {
        c2 x = {42};
        c2 & y = x;
        BOOST_TEST(check<__LINE__>(x,"c2"));
        BOOST_TEST(check<__LINE__>(y,"c2"));
    }
    BOOST_TEST(check<__LINE__>(c3{42},"c3"));
    BOOST_TEST(check<__LINE__>(c3{42},"42"));
    {
        c3 x = {42};
        c3 & y = x;
        BOOST_TEST(check<__LINE__>(x,"c3"));
        BOOST_TEST(check<__LINE__>(x,"42"));
        BOOST_TEST(check<__LINE__>(y,"c3"));
        BOOST_TEST(check<__LINE__>(y,"42"));
    }
    BOOST_TEST(check<__LINE__>(c4(),"c4"));
    BOOST_TEST(check<__LINE__>(c4(),"{not printable}"));
    {
        c4 x;
        c4 & y = x;
        BOOST_TEST(check<__LINE__>(x,"c4"));
        BOOST_TEST(check<__LINE__>(x,"{not printable}"));
        BOOST_TEST(check<__LINE__>(y,"c4"));
        BOOST_TEST(check<__LINE__>(y,"{not printable}"));
    }
    BOOST_TEST(check<__LINE__>(my_exception{}, "std::exception::what(): my_exception_what"));
    return boost::report_errors();
}

#endif
