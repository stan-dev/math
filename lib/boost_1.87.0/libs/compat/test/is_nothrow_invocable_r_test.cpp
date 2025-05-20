// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/invoke.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>

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

struct Y
{
    Y( int );
    Y( double ) noexcept;
};

int main()
{
    using boost::compat::is_nothrow_invocable_r;

    // nonfunction

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int, int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int, int, int> ));

    // function reference

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(&)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(&)(int), char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int(&)(int, int), int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, double(&)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, void(&)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(&)(int), char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int(&)(int, int), int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, double(&)(double, double, double), float, float, float> ));

#if defined(__cpp_noexcept_function_type)

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(&)() noexcept> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, char(&)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int(&)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, double(&)(double, double, double) noexcept, float, float, float> ));

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, void(&)() noexcept> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, char(&)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int(&)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, double(&)(double, double, double) noexcept, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, void(&)() noexcept> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, char(&)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int(&)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, double(&)(double, double, double) noexcept, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, void(&)() noexcept> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, char(&)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int(&)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double(&)(double, double, double) noexcept, float, float, float> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(&)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(&)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(&)(int), int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, void(&)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(&)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(&)(int), int, int> ));

    // function pointer

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(*)(int), char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int(*)(int, int), int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, double(*)(double, double, double), float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, void(*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(*)(int), char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int(*)(int, int), int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, double(*)(double, double, double), float, float, float> ));

#if defined(__cpp_noexcept_function_type)

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(*)() noexcept> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, char(*)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int(*)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, double(*)(double, double, double) noexcept, float, float, float> ));

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, void(*)() noexcept> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, char(*)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int(*)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, double(*)(double, double, double) noexcept, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, void(*)() noexcept> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, char(*)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int(*)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, double(*)(double, double, double) noexcept, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, void(*)() noexcept> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, char(*)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int(*)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double(*)(double, double, double) noexcept, float, float, float> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(*)(int), int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, void(*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(*)(int), int, int> ));

    // object

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, F> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, F, int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, F> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, F, int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, F> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, F, int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, F> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, F, int, int> ));

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1910)

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, F, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, F, float, float, float> ));

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, F, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, F, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, F, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, F, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, F, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, F, float, float, float> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, F, int, int, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, F const> ));

    // member function pointer

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(X::*)(), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(X::*)(int), X, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int(X::*)(int, int), X, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, double(X::*)(double, double, double), X, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, void(X::*)(), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(X::*)(int), X, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int(X::*)(int, int), X, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, double(X::*)(double, double, double), X, float, float, float> ));

#if defined(__cpp_noexcept_function_type)

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(X::*)() noexcept, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, char(X::*)(int) noexcept, X, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int(X::*)(int, int) noexcept, X, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, double(X::*)(double, double, double) noexcept, X, float, float, float> ));

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, void(X::*)() noexcept, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, char(X::*)(int) noexcept, X, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int(X::*)(int, int) noexcept, X, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, double(X::*)(double, double, double) noexcept, X, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, void(X::*)() noexcept, X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, char(X::*)(int) noexcept, X, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int(X::*)(int, int) noexcept, X, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, double(X::*)(double, double, double) noexcept, X, float, float, float> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, void(X::*)() noexcept, X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, char(X::*)(int) noexcept, X, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int(X::*)(int, int) noexcept, X, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double(X::*)(double, double, double) noexcept, X, float, float, float> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(X::*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(X::*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, void(X::*)(), X, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(X::*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(X::*)(int), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(X::*)(int), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, char(X::*)(int), X, int, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, void(X::*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, void(X::*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, void(X::*)(), X, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(X::*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(X::*)(int), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(X::*)(int), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, char(X::*)(int), X, int, int> ));

    // member data pointer

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1910)

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int X::*, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int X::*, X&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int X::*, X*> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<long, int X::*, X const*> ));

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int X::*, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int X::*, X&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int X::*, X*> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<void, int X::*, X const*> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int X::*, X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int X::*, X const> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int X::*, X&> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int X::*, X const&> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int X::*, X*> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<X, int X::*, X const*> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int X::*, X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int X::*, X const> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int X::*, X&> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int X::*, X const&> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int X::*, X*> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<Y, int X::*, X const*> ));

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double X::*, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double X::*, X&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double X::*, X*> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<Y, double X::*, X const*> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<int&, int X::*, X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<int&, int X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<int&, int X::*, X&> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<int&, int X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<int&, int X::*, X*> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<int&, int X::*, X const*> ));

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable_r<int&&, int X::*, X> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int X::*> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int X::*, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<long, int X::*, X, int> ));

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int X::*> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int X::*, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable_r<void, int X::*, X, int> ));

    return boost::report_errors();
}
