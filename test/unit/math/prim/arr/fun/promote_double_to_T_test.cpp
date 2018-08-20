#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

TEST(MathFunctions, promote_double_to_T_std_vector_int) {
  std::vector<int> x(3);
  EXPECT_TRUE((std::is_same<const std::vector<int>&,
                            decltype(stan::math::promote_double_to_T<double>(
                                x))>::value));
  EXPECT_FALSE((std::is_same<std::vector<double>,
                             decltype(stan::math::promote_double_to_T<double>(
                                 x))>::value));
}

TEST(MathFunctions, promote_double_to_T_std_vector_double) {
  std::vector<double> x(3);
  EXPECT_TRUE((std::is_same<std::vector<float>,
                            decltype(stan::math::promote_double_to_T<float>(
                                x))>::value));
  EXPECT_FALSE((std::is_same<const std::vector<double>&,
                             decltype(stan::math::promote_double_to_T<float>(
                                 x))>::value));
}
