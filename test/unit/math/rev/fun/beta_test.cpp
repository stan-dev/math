#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, beta_map) {
  using stan::math::beta;
  using stan::math::matrix_v;
  matrix_v in1 = matrix_v::Random(100, 100);
  matrix_v out = beta(in1, in1);
  EXPECT_MATRIX_EQ(out.val(), beta(in1.val(), in1.val()));
}
