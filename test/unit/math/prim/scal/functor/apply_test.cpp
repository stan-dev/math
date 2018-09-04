#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

struct func {
  template <typename T>
  T operator()(T t) {
    return t;
  }
};

TEST(MathFunctions, apply_basic_double) {
  std::tuple<double> x = {1.0};

  auto y = stan::math::apply(func{}, x);

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_basic_int) {
  std::tuple<int> x = {1};

  auto y = stan::math::apply(func{}, x);

  EXPECT_EQ(1, y);
  EXPECT_TRUE((std::is_same<int, decltype(y)>::value));
}

TEST(MathFunctions, apply_const_double) {
  const std::tuple<const double> x = {1.0};

  auto y = stan::math::apply(func{}, x);

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_double) {
  auto y = stan::math::apply(func{}, std::make_tuple(1.0));

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function) {
  auto y = stan::math::apply([](auto x) { return x; }, std::make_tuple(1.0));

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function_reference) {
  auto y = stan::math::apply([](auto& x) { return x; }, std::make_tuple(1.0));

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function_const_reference) {
  auto y = stan::math::apply([](const auto& x) { return x; },
                             std::make_tuple(1.0));

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function_return_reference) {
  std::tuple<double> x = {1.0};

  decltype(auto) y = stan::math::apply([](auto& x) -> auto& { return x; }, x);

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double&, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function_return_const_reference) {
  std::tuple<double> x = {1.0};

  decltype(auto) y
      = stan::math::apply([](const auto& x) -> const auto& { return x; }, x);

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<const double&, decltype(y)>::value));
}
