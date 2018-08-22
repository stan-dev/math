#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

TEST(MathFunctions, promote_double_to_T_int) {
  int x;
  EXPECT_TRUE((std::is_same<std::tuple<>,
                            decltype(stan::math::promote_double_to_T<double>(
                                std::make_tuple(x)))>::value));
  EXPECT_FALSE((std::is_same<std::tuple<double>,
                             decltype(stan::math::promote_double_to_T<double>(
                                 std::make_tuple(x)))>::value));
}

TEST(MathFunctions, promote_double_to_T_double) {
  double x;
  EXPECT_TRUE((std::is_same<std::tuple<float>,
                            decltype(stan::math::promote_double_to_T<float>(
                                std::make_tuple(x)))>::value));
  EXPECT_FALSE((std::is_same<std::tuple<>,
                             decltype(stan::math::promote_double_to_T<float>(
                                 std::make_tuple(x)))>::value));
}

TEST(MathFunctions, promote_double_to_T_mix) {
  int xi;
  double xd;
  auto a = stan::math::promote_double_to_T<float>(std::make_tuple(xd, xi, xd));
  EXPECT_TRUE((std::is_same<std::tuple<float, float>,
                            decltype(stan::math::promote_double_to_T<float>(
                                std::make_tuple(xd, xi, xd)))>::value));
  EXPECT_TRUE((std::is_same<std::tuple<float>,
                            decltype(stan::math::promote_double_to_T<float>(
                                std::make_tuple(xd, xi)))>::value));
}

TEST(MathFunctions, promote_double_to_T_std_vector_int) {
  std::vector<int> x(3);
  EXPECT_TRUE((std::is_same<std::tuple<>,
                            decltype(stan::math::promote_double_to_T<double>(
                                std::make_tuple(x)))>::value));
  EXPECT_FALSE((std::is_same<std::tuple<std::vector<double> >,
                             decltype(stan::math::promote_double_to_T<double>(
                                 std::make_tuple(x)))>::value));
}

TEST(MathFunctions, promote_double_to_T_std_vector_double) {
  std::vector<double> x(3);
  EXPECT_TRUE((std::is_same<std::tuple<std::vector<float> >,
                            decltype(stan::math::promote_double_to_T<float>(
                                std::make_tuple(x)))>::value));
  EXPECT_FALSE((std::is_same<std::tuple<>,
                             decltype(stan::math::promote_double_to_T<float>(
                                 std::make_tuple(x)))>::value));
}

TEST(MathFunctions, promote_double_to_T_std_vector_mix) {
  std::vector<int> xi(3);
  std::vector<double> xd(3);
  EXPECT_TRUE((std::is_same<std::tuple<std::vector<float>, std::vector<float> >,
                            decltype(stan::math::promote_double_to_T<float>(
                                std::make_tuple(xd, xi, xd)))>::value));
  EXPECT_TRUE((std::is_same<std::tuple<std::vector<float> >,
                            decltype(stan::math::promote_double_to_T<float>(
                                std::make_tuple(xd, xi)))>::value));
}
