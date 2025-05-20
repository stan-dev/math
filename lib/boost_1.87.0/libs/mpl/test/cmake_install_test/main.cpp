// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

#include <boost/mpl/assert.hpp>

int main()
{
    BOOST_MPL_ASSERT_RELATION( 1, <, 5 );
}
