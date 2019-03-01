#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/scal/meta/is_var_or_arithmetic.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMeta, is_var_or_arithmetic_simple) {
  using stan::is_var_or_arithmetic;
  EXPECT_TRUE(stan::is_var_or_arithmetic<stan::math::var>::value);
  EXPECT_FALSE(stan::is_var_or_arithmetic<stan::math::var&>::value);
  bool temp
      = is_var_or_arithmetic<stan::math::var, std::vector<stan::math::var>,
                             stan::math::var, stan::math::var, stan::math::var,
                             stan::math::var>::value;
  EXPECT_TRUE(temp);
  temp = is_var_or_arithmetic<std::vector<stan::math::var>, stan::math::var,
                              std::vector<stan::math::var&>>::value;
  EXPECT_FALSE(temp);
}
