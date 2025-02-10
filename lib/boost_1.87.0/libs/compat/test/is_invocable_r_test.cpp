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
    using boost::compat::is_invocable_r;

    // nonfunction

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, int, int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, int, int, int> ));

    // function reference

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, void(&)()> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, char(&)(int), char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int(&)(int, int), int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, double(&)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, void(&)()> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, char(&)(int), char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int(&)(int, int), int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, double(&)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, void(&)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, char(&)(int), char> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int(&)(int, int), int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, double(&)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, void(&)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, char(&)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, char(&)(int), int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, void(&)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, char(&)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, char(&)(int), int, int> ));

    // function pointer

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, void(*)()> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, char(*)(int), char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int(*)(int, int), int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, double(*)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, void(*)()> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, char(*)(int), char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int(*)(int, int), int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, double(*)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, void(*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, char(*)(int), char> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int(*)(int, int), int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, double(*)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, void(*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, char(*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, char(*)(int), int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, void(*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, char(*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, char(*)(int), int, int> ));

    // object

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, F> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, F, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, F, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, F, float, float, float> ));

    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, F> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, F, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, F, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, F, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, F> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, F, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, F, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, F, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, F, int, int, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, F const> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, F, int, int, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, F const> ));

    // member function pointer

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, void(X::*)(), X> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, char(X::*)(int), X, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int(X::*)(int, int), X, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, double(X::*)(double, double, double), X, float, float, float> ));

    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, void(X::*)(), X> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, char(X::*)(int), X, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int(X::*)(int, int), X, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, double(X::*)(double, double, double), X, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, void(X::*)(), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, char(X::*)(int), X, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int(X::*)(int, int), X, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, double(X::*)(double, double, double), X, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, void(X::*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, void(X::*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, void(X::*)(), X, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, char(X::*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, char(X::*)(int), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, char(X::*)(int), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, char(X::*)(int), X, int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, void(X::*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, void(X::*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, void(X::*)(), X, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, char(X::*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, char(X::*)(int), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, char(X::*)(int), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, char(X::*)(int), X, int, int> ));

    // member data pointer

    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int X::*, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int X::*, X&> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int X::*, X*> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<long, int X::*, X const*> ));

    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int X::*, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int X::*, X&> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int X::*, X*> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<void, int X::*, X const*> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int X::*, X> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int X::*, X const> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int X::*, X&> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int X::*, X const&> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int X::*, X*> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<X, int X::*, X const*> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<int&, int X::*, X> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<int&, int X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<int&, int X::*, X&> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<int&, int X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<int&, int X::*, X*> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<int&, int X::*, X const*> ));

    BOOST_TEST_TRAIT_TRUE(( is_invocable_r<int&&, int X::*, X> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, int X::*> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, int X::*, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<long, int X::*, X, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, int X::*> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, int X::*, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_invocable_r<void, int X::*, X, int> ));

    return boost::report_errors();
}
