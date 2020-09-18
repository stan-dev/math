#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, exp_map) {
  using stan::math::exp;
  using stan::math::var;
  using stan::math::matrix_v;
  matrix_v in1 = matrix_v::Random(1000, 1000);
  matrix_v out = exp(in1);
  EXPECT_MATRIX_EQ(out.val(),in1.val().array().exp());
}
