/* test_splitmix64.cpp
 *
 * Copyright Steven Watanabe 2011
 * Copyright Matt Borland 2022
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * $Id$
 *
 */

#include <boost/random/splitmix64.hpp>
#include <cstdint>
#include <cmath>

#define BOOST_RANDOM_URNG boost::random::splitmix64

// principal operation validated with CLHEP, values by experiment
#define BOOST_RANDOM_VALIDATION_VALUE UINT64_C(542758903869407163)
#define BOOST_RANDOM_SEED_SEQ_VALIDATION_VALUE UINT64_C(7044339178176046395)

#include "test_generator.ipp"
