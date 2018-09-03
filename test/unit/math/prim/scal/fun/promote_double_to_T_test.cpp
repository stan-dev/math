#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <type_traits>

TEST(MathFunctions, promote_double_to_T_int) {
  int x = 3;
  auto y = stan::math::promote_double_to_T<double>(x);
  EXPECT_TRUE((std::is_same<const int&, decltype(y)>::value));
  EXPECT_FALSE((std::is_same<double, decltype(y)>::value));
  EXPECT_EQ(x, y);
}

TEST(MathFunctions, promote_double_to_T_double) {
  double x = 3.0;
  auto y = stan::math::promote_double_to_T<float>(x);
  EXPECT_TRUE((std::is_same<float, decltype(y)>::value));
  EXPECT_FALSE((std::is_same<const double&, decltype(y)>::value));
  EXPECT_FLOAT_EQ(static_cast<float>(x), y);
}
