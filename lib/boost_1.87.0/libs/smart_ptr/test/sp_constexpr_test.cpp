//
// sp_constexpr_test.cpp
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
};

struct Z
{
    Z();
};

static Z z;

static boost::shared_ptr<X> p1;
static boost::weak_ptr<X> p2;

static boost::shared_ptr<X> p3( nullptr );

Z::Z()
{
    p1.reset( new X );
    p2 = p1;
    p3.reset( new X );
}

int main()
{
    BOOST_TEST( p1.get() != 0 );
    BOOST_TEST_EQ( p1.use_count(), 1 );

    BOOST_TEST_EQ( p2.use_count(), 1 );
    BOOST_TEST_EQ( p2.lock(), p1 );

    BOOST_TEST( p3.get() != 0 );
    BOOST_TEST_EQ( p3.use_count(), 1 );

    return boost::report_errors();
}

#endif
