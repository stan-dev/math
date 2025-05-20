// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/invoke.hpp>
#include <boost/core/lightweight_test_trait.hpp>

struct F
{
    void operator()()
    {
    }

    char operator()( char x1 ) noexcept
    {
        return x1;
    }

    int operator()( int x1, int x2 ) const
    {
        return 10*x1+x2;
    }

    double operator()( float x1, float x2, float x3 ) const noexcept
    {
        return 100*x1 + 10*x2 + x3;
    }
};

struct X
{
};

int main()
{
    using boost::compat::is_invocable;

    // nonfunction

    BOOST_TEST_TRAIT_FALSE(( is_invocable<int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<int, int, int> ));

    // function reference

    BOOST_TEST_TRAIT_TRUE(( is_invocable<void(&)()> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<char(&)(int), char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<int(&)(int, int), int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<double(&)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable<void(&)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<char(&)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<char(&)(int), int, int> ));

    // function pointer

    BOOST_TEST_TRAIT_TRUE(( is_invocable<void(*)()> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<char(*)(int), char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<int(*)(int, int), int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<double(*)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable<void(*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<char(*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<char(*)(int), int, int> ));

    // object

    BOOST_TEST_TRAIT_TRUE(( is_invocable<F> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<F, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<F, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<F, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable<F, int, int, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<F const> ));

    // member function pointer

    BOOST_TEST_TRAIT_TRUE(( is_invocable<void(X::*)(), X> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<char(X::*)(int), X, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<int(X::*)(int, int), X, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<double(X::*)(double, double, double), X, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable<void(X::*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<void(X::*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<void(X::*)(), X, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<char(X::*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<char(X::*)(int), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<char(X::*)(int), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<char(X::*)(int), X, int, int> ));

    // member data pointer

    BOOST_TEST_TRAIT_TRUE(( is_invocable<int X::*, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<int X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<int X::*, X&> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<int X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<int X::*, X*> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable<int X::*, X const*> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable<int X::*> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<int X::*, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable<int X::*, X, int> ));

    return boost::report_errors();
}
