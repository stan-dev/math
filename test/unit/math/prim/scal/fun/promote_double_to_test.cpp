#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <type_traits>

TEST(MathFunctions, promote_double_to_int) {
  int x = 3;
  auto y = stan::math::promote_double_to<double>(x);
  EXPECT_TRUE((std::is_same<const int&, decltype(y)>::value));
  EXPECT_FALSE((std::is_same<double, decltype(y)>::value));
  EXPECT_EQ(x, y);
}

TEST(MathFunctions, promote_double_to_double) {
  double x = 3.0;
  auto y = stan::math::promote_double_to<float>(x);
  EXPECT_TRUE((std::is_same<float, decltype(y)>::value));
  EXPECT_FALSE((std::is_same<const double&, decltype(y)>::value));
  EXPECT_FLOAT_EQ(static_cast<float>(x), y);
}
