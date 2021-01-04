#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

TEST(MathFunctions, for_each_basic_unary) {
  int count = 0;

  auto f = [&count](const auto& x, const auto& i) { count++; };

  std::tuple<> x1;
  stan::math::for_each(f, x1);

  EXPECT_EQ(count, 0);

  auto x2 = std::make_tuple(1.0, 2.0);
  stan::math::for_each(f, x2);

  EXPECT_EQ(count, 2);
}

TEST(MathFunctions, for_each_basic_unary_index) {
  std::vector<int> v = {-5, 2};
  auto x = std::make_tuple(v[0], v[1]);

  auto f = [&v](const auto& y, const auto& i) { EXPECT_EQ(y, v[i]); };

  stan::math::for_each(f, x);
}

TEST(MathFunctions, for_each_basic_binary) {
  int count = 0;

  auto f = [&count](const auto& x, const auto& y, const auto& i) { count++; };

  std::tuple<> x1;
  std::tuple<> y1;
  stan::math::for_each(f, x1, y1);

  EXPECT_EQ(count, 0);

  auto x2 = std::make_tuple(1, 2.0);
  auto y2 = std::make_tuple(3.0, 5);
  stan::math::for_each(f, x2, y2);

  EXPECT_EQ(count, 2);

  EXPECT_THROW(stan::math::for_each(f, x1, x2), std::invalid_argument);
}

TEST(MathFunctions, for_each_basic_binary_index) {
  std::vector<int> v1 = {-5, 2};
  std::vector<int> v2 = {3, -4};
  auto x1 = std::make_tuple(v1[0], v1[1]);
  auto x2 = std::make_tuple(v2[0], v2[1]);

  auto f = [&v1, &v2](const auto& y, const auto& z, const auto& i) {
    EXPECT_EQ(y, v1[i]);
    EXPECT_EQ(z, v2[i]);
  };

  stan::math::for_each(f, x1, x2);
}
