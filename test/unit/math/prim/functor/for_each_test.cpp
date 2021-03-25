#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

TEST(MathFunctions, for_each_basic_unary) {
  int count = 0;

  auto f = [&count](const auto& x) { count++; };

  std::tuple<> x1;
  stan::math::for_each(f, x1);

  EXPECT_EQ(count, 0);

  auto x2 = std::make_tuple(1.0, 2.0);
  stan::math::for_each(f, x2);

  EXPECT_EQ(count, 2);
}

TEST(MathFunctions, for_each_basic_unary_index) {
  std::vector<int> v = {-5, 2};
  auto x = std::make_tuple(v, v);

  auto f = [](auto& y) {
    y[0] += 1;
    y[1] += 1;
    return;
  };

  stan::math::for_each(f, x);
  EXPECT_EQ(std::get<0>(x)[0], -4);
  EXPECT_EQ(std::get<0>(x)[1], 3);
  EXPECT_EQ(std::get<1>(x)[0], -4);
  EXPECT_EQ(std::get<1>(x)[1], 3);
}

TEST(MathFunctions, for_each_basic_binary) {
  int count = 0;

  auto f = [&count](const auto& x, const auto& y) { count++; };

  std::tuple<> x1;
  std::tuple<> y1;
  stan::math::for_each(f, x1, y1);

  EXPECT_EQ(count, 0);

  auto x2 = std::make_tuple(1, 2.0);
  auto y2 = std::make_tuple(3.0, 5);
  stan::math::for_each(f, x2, y2);

  EXPECT_EQ(count, 2);
}

TEST(MathFunctions, for_each_basic_binary_index) {
  std::vector<int> v1 = {-5, 2};
  auto x1 = std::make_tuple(v1[0], v1[1]);
  auto f = [](const auto& y, const auto& z) { EXPECT_EQ(y, z); };

  stan::math::for_each(f, x1, x1);
}

TEST(MathFunctions, for_each_basic_trinary) {
  int count = 0;

  auto f = [&count](const auto& x, const auto& y, const auto& z) { count++; };

  std::tuple<> x1;
  std::tuple<> y1;
  std::tuple<> z1;
  stan::math::for_each(f, x1, y1, z1);

  EXPECT_EQ(count, 0);

  auto x2 = std::make_tuple(1, 2.0);
  auto y2 = std::make_tuple(3.0, 5);
  auto z2 = std::make_tuple(8.0, 3);
  stan::math::for_each(f, x2, y2, z2);

  EXPECT_EQ(count, 2);
}

TEST(MathFunctions, for_each_basic_trinary_index) {
  std::vector<int> v1 = {-5, 2, 3};
  auto x1 = std::make_tuple(v1[0], v1[1], v1[2]);
  auto f = [](const auto& x, const auto& y, const auto& z) {
    EXPECT_EQ(x, z);
    EXPECT_EQ(x, y);
  };

  stan::math::for_each(f, x1, x1, x1);
}
