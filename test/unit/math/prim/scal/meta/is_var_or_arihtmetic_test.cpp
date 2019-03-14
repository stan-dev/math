#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/meta/is_var_or_arithmetic.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMeta, is_var_or_arithmetic_simple) {
  using stan::is_var_or_arithmetic;
  EXPECT_TRUE(stan::is_var_or_arithmetic<double>::value);
  EXPECT_FALSE(stan::is_var_or_arithmetic<double&>::value);
  EXPECT_FALSE(stan::is_var_or_arithmetic<std::vector<double>>::value);

  EXPECT_FALSE(stan::is_var_or_arithmetic<int&>::value);
  EXPECT_TRUE(stan::is_var_or_arithmetic<int>::value);

  bool temp = is_var_or_arithmetic<int&, double>::value;
  EXPECT_FALSE(temp);
  temp = is_var_or_arithmetic<double, std::vector<double>>::value;
  EXPECT_FALSE(temp);
  temp = is_var_or_arithmetic<double, double, std::vector<double>, double,
                              std::vector<double>, double&>::value;
  EXPECT_FALSE(temp);
}
