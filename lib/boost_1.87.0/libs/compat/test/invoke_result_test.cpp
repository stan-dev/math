// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/invoke.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/mp11/utility.hpp>

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
    using boost::compat::invoke_result_t;
    using boost::mp11::mp_valid;

    // nonfunction

    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, int, int, int> ));

    // function reference

    BOOST_TEST_TRAIT_SAME( invoke_result_t<void(&)()>, void );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<char(&)(int), char>, char );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<int(&)(int, int), int, int>, int );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<double(&)(double, double, double), float, float, float>, double );

    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, void(&)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, char(&)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, char(&)(int), int, int> ));

    // function pointer

    BOOST_TEST_TRAIT_SAME( invoke_result_t<void(*)()>, void );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<char(*)(int), char>, char );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<int(*)(int, int), int, int>, int );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<double(*)(double, double, double), float, float, float>, double );

    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, void(*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, char(*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, char(*)(int), int, int> ));

    // object

    BOOST_TEST_TRAIT_SAME( invoke_result_t<F>, void );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<F, char>, char );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<F, int, int>, int );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<F, float, float, float>, double );

    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, F, int, int, int, int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, F const> ));

    // member function pointer

    BOOST_TEST_TRAIT_SAME( invoke_result_t<void(X::*)(), X>, void );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<char(X::*)(int), X, char>, char );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<int(X::*)(int, int), X, int, int>, int );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<double(X::*)(double, double, double), X, float, float, float>, double );

    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, void(X::*)()> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, void(X::*)(), int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, void(X::*)(), X, int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, char(X::*)(int)> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, char(X::*)(int), int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, char(X::*)(int), X> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, char(X::*)(int), X, int, int> ));

    // member data pointer

    BOOST_TEST_TRAIT_SAME( invoke_result_t<int X::*, X>, int&& );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<int X::*, X const>, int const&& );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<int X::*, X&>, int& );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<int X::*, X const&>, int const& );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<int X::*, X*>, int& );
    BOOST_TEST_TRAIT_SAME( invoke_result_t<int X::*, X const*>, int const& );

    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, int X::*> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, int X::*, int> ));
    BOOST_TEST_TRAIT_FALSE(( mp_valid<invoke_result_t, int X::*, X, int> ));

    return boost::report_errors();
}
