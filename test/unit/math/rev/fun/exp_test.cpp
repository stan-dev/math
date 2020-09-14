#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, expInt) {
  using stan::math::exp;
  EXPECT_FLOAT_EQ(std::exp(3), exp(3));
  EXPECT_FLOAT_EQ(std::exp(3.1), exp(3.1));
  EXPECT_FLOAT_EQ(std::exp(3.0), exp(3.0));
}

TEST(MathFunctions, exp_map) {
  using stan::math::exp;
  using stan::math::var;
  using stan::math::vector_v;
  vector_v in1 = vector_v::Random(100000);
  vector_v out = exp(in1);
}
