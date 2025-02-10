//
// sp_constexpr_test2.cpp
//
// Copyright 2017 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//

#include <boost/config.hpp>
#include <boost/config/workaround.hpp>
#include <boost/config/pragma_message.hpp>
#include <boost/config/helper_macros.hpp>

#if BOOST_WORKAROUND( BOOST_MSVC, < 1930 )

// MSVC does not implement static initialization for constexpr constructors
BOOST_PRAGMA_MESSAGE("Skipping test due to BOOST_MSVC < 1930")
int main() {}

#elif defined(__clang__) && defined( BOOST_NO_CXX14_CONSTEXPR )

// Clang only implements static initialization for constexpr in C++14 mode
BOOST_PRAGMA_MESSAGE("Skipping test due to __clang__ and BOOST_NO_CXX14_CONSTEXPR")
int main() {}

#elif defined( _LIBCPP_VERSION ) && ( _LIBCPP_VERSION < 6000 )

// in libc++, atomic_flag has a non-constexpr constructor from bool
BOOST_PRAGMA_MESSAGE("Skipping test due to _LIBCPP_VERSION " BOOST_STRINGIZE(_LIBCPP_VERSION))
int main() {}

#else

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/core/lightweight_test.hpp>

struct X: public boost::enable_shared_from_this<X>
{
    int v_;

    constexpr X() noexcept: v_( 1 )
    {
    }
};

struct Z
{
    Z();
};

static Z z;
static X x;

Z::Z()
{
    BOOST_TEST_EQ( x.v_, 1 );
}

int main()
{
    return boost::report_errors();
}

#endif
