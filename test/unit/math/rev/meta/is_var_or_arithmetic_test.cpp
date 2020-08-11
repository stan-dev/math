#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaRevScal, is_var_or_arithmetic_simple) {
  using stan::is_var_or_arithmetic;
  EXPECT_TRUE(stan::is_var_or_arithmetic<stan::math::var>::value);
  EXPECT_TRUE(stan::is_var_or_arithmetic<stan::math::var&>::value);
  bool temp
      = is_var_or_arithmetic<stan::math::var, std::vector<stan::math::var>,
                             stan::math::var, stan::math::var, stan::math::var,
                             stan::math::var>::value;
  EXPECT_TRUE(temp);
  temp = is_var_or_arithmetic<std::vector<stan::math::var>, stan::math::var,
                              std::vector<stan::math::var const*>>::value;
  EXPECT_FALSE(temp);
  temp = is_var_or_arithmetic<std::vector<stan::math::var>, stan::math::var,
                              std::vector<stan::math::var>>::value;
  EXPECT_TRUE(temp);
}
