#ifndef STAN_TEST_UNIT_MATH_MIX_UTIL_HPP
#define STAN_TEST_UNIT_MATH_MIX_UTIL_HPP

#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

struct mathMix : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
  void TearDown() {
    // make sure memory's clean after each test
    stan::math::recover_memory();
  }
};

#endif
