#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <type_traits>

TEST(MathFunctions, promote_double_to_T_int) {
  int x = 3;
  EXPECT_TRUE((
      std::is_same<const int&, decltype(stan::math::promote_double_to_T<double>(
                                   x))>::value));
  EXPECT_FALSE(
      (std::is_same<double, decltype(stan::math::promote_double_to_T<double>(
                                x))>::value));
}

TEST(MathFunctions, promote_double_to_T_double) {
  double x = 3.0;
  EXPECT_TRUE(
      (std::is_same<float, decltype(stan::math::promote_double_to_T<float>(
                               x))>::value));
  EXPECT_FALSE((std::is_same<const double&,
                             decltype(stan::math::promote_double_to_T<float>(
                                 x))>::value));
}
