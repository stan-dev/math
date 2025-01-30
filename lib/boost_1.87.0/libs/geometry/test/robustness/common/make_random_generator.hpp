// Boost.Geometry (aka GGL, Generic Geometry Library)
// Robustness Test

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_TEST_ROBUSTNESS_COMMON_MAKE_RANDOM_GENERATOR_HPP
#define BOOST_GEOMETRY_TEST_ROBUSTNESS_COMMON_MAKE_RANDOM_GENERATOR_HPP

#include <functional>
#include <random>

inline auto make_int_generator(int seed, int size)
{
    return std::bind(std::uniform_int_distribution<int>(0, size - 1),
              std::default_random_engine(seed == -1 ? std::random_device()() : seed));    
}

inline auto make_real_generator(int seed, double a, double b)
{
    return std::bind(std::uniform_real_distribution<double>(a, b),
              std::default_random_engine(seed == -1 ? std::random_device()() : seed));    
}

#endif // BOOST_GEOMETRY_TEST_ROBUSTNESS_COMMON_MAKE_RANDOM_GENERATOR_HPP
