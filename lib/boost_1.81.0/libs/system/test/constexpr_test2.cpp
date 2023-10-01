// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/error_code.hpp>
#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>

#if defined(BOOST_GCC) && BOOST_GCC >= 40700 && BOOST_GCC < 40800

BOOST_PRAGMA_MESSAGE("Skipping test, BOOST_GCC is 407xx")

#else

struct X
{
    boost::system::error_code ec;
};

X const& f()
{
    BOOST_STATIC_CONSTEXPR X x;
    return x;
}

#endif
