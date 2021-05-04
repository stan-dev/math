#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, is_var_or_arithmetic_simple) {
  using stan::is_var_or_arithmetic;
  EXPECT_TRUE(stan::is_var_or_arithmetic<double>::value);
  EXPECT_TRUE(stan::is_var_or_arithmetic<double&>::value);
  EXPECT_FALSE(stan::is_var_or_arithmetic<std::vector<double const*>>::value);

  EXPECT_TRUE(stan::is_var_or_arithmetic<int&>::value);
  EXPECT_TRUE(stan::is_var_or_arithmetic<int>::value);

  bool temp = is_var_or_arithmetic<int&, double>::value;
  EXPECT_TRUE(temp);
  temp = is_var_or_arithmetic<double, std::vector<double>>::value;
  EXPECT_TRUE(temp);
  temp = is_var_or_arithmetic<double, double, std::vector<double>, double,
                              std::vector<double>, double&>::value;
  EXPECT_TRUE(temp);
}
