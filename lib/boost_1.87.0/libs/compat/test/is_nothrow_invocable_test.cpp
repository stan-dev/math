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

int main()
{
    using boost::compat::is_nothrow_invocable;

    // nonfunction

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int, int, int> ));

    // function reference

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<void(&)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(&)(int), char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int(&)(int, int), int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<double(&)(double, double, double), float, float, float> ));

#if defined(__cpp_noexcept_function_type)

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<void(&)() noexcept> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<char(&)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int(&)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<double(&)(double, double, double) noexcept, float, float, float> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<void(&)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(&)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(&)(int), int, int> ));

    // function pointer

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<void(*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(*)(int), char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int(*)(int, int), int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<double(*)(double, double, double), float, float, float> ));

#if defined(__cpp_noexcept_function_type)

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<void(*)() noexcept> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<char(*)(int) noexcept, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int(*)(int, int) noexcept, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<double(*)(double, double, double) noexcept, float, float, float> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<void(*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(*)(int), int, int> ));

    // object

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<F> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<F, int, int> ));

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1910)

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<F, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<F, float, float, float> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<F, int, int, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<F const> ));

    // member function pointer

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<void(X::*)(), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(X::*)(int), X, char> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int(X::*)(int, int), X, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<double(X::*)(double, double, double), X, float, float, float> ));

#if defined(__cpp_noexcept_function_type)

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<void(X::*)() noexcept, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<char(X::*)(int) noexcept, X, char> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int(X::*)(int, int) noexcept, X, int, int> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<double(X::*)(double, double, double) noexcept, X, float, float, float> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<void(X::*)()> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<void(X::*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<void(X::*)(), X, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(X::*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(X::*)(int), int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(X::*)(int), X> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<char(X::*)(int), X, int, int> ));

    // member data pointer

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1910)

    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int X::*, X> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int X::*, X const> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int X::*, X&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int X::*, X const&> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int X::*, X*> ));
    BOOST_TEST_TRAIT_TRUE(( is_nothrow_invocable<int X::*, X const*> ));

#endif

    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int X::*> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int X::*, int> ));
    BOOST_TEST_TRAIT_FALSE(( is_nothrow_invocable<int X::*, X, int> ));

    return boost::report_errors();
}
