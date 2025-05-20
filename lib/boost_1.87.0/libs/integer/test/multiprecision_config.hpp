// Copyright (c) 2018, 2024 Andrey Semashev
//
// Use, modification, and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// https://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_INTEGER_TEST_MULTIPRECISION_CONFIG_HPP_INCLUDED_
#define BOOST_INTEGER_TEST_MULTIPRECISION_CONFIG_HPP_INCLUDED_

#include <boost/config.hpp>

// Boost.Multiprecision requires a number of C++11 features, see boost/multiprecision/detail/check_cpp11_config.hpp.
// Also, Boost.Multiprecision internally uses Boost.Math, which requires C++14 and a recent enough MSVC, see boost/math/tools/config.hpp.
#if (BOOST_CXX_VERSION < 201402) || \
    (defined(_MSC_VER) && (_MSC_VER <= 1900))
#define DISABLE_MP_TESTS
#endif

#endif // BOOST_INTEGER_TEST_MULTIPRECISION_CONFIG_HPP_INCLUDED_
