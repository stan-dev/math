#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, exp_map) {
  using stan::math::exp;
  using stan::math::var;
  using stan::math::vector_v;
  vector_v in1 = vector_v::Random(100000);
  vector_v out = exp(in1);
}
