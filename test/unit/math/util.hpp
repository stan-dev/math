#ifndef TEST_UNIT_MATH_UTIL_HPP
#define TEST_UNIT_MATH_UTIL_HPP
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

struct TestUnitMathTestAd : public testing::Test {
  inline void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
};

#endif
