#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, pow_map) {
  using stan::math::pow;
  using stan::math::var;
  using stan::math::matrix_v;
  matrix_v in1 = matrix_v::Random(1000, 1000);
  matrix_v out = pow(in1, in1);
  EXPECT_MATRIX_EQ(out.val(),pow(in1.val(), in1.val()));
}
